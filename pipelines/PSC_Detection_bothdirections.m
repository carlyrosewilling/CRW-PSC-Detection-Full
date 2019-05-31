%% Pipeline for Bayesian Algorithm Without Sampler: Seul Ah iPSC and ePSC at same time%%

%Written by CRW, 12 Dec 2018
    %last updated: 13 Dec 2018
    
    %line 93 and 136 reflect local template paths -- make sure to change
    
    %lines 74, 77, and 85 are important if your ePSC events seem wrong
    %lines 116, 119, and 127 are important if your iPSC events seem wrong
        %***See comments and powerpoint for more detailed information on how
        %to adjust these accordingly***

%% Turn off dumb warning for loading "wave" struct and directory
warning('off', 'MATLAB:unknownObjectNowStruct');
warning('off', 'MATLAB:MKDIR:DirectoryExists');
warning('off', 'signal:findpeaks:largeMinPeakHeight');

%% Initialize %%

%User selects exact file where data is stored
    [~,prepath] = uiputfile('*.*','Select data folder', 'datapath.mat');

    if isnumeric(prepath)
        disp('User must select an input');
        return
    end
    
%Finds all files in input folder which are a .mat and creates a list
	Files=dir(fullfile(prepath,'*.mat'));
    
%User selects where to save output files
    [~, savePath]=uiputfile('output.mat', 'Select path for output');

    
%% Organize by file 
    %loads individual file of traces
    for i = 1:length(Files)
        load(strcat(Files(i).folder,'/', Files(i).name));
        baseline = eval(Files(i).name(1:end-4));
        [nACQ tracelength] = size(baseline.data);
        
        %Creates directory for output path according to original file loaded.
        savePath1 = fullfile(savePath, Files(i).name(1:end-4));
        mkdir(savePath1)
            
    %% Runs through every trace stored in each file
        for j = 1:nACQ
            %Creates string identifying which number trace this is -- used to create variable name
            string = strcat('trace_', num2str(j)); 
            params = [];
            params.savepath = savePath;
            params.traces_file = string;
            
            %Isolates data for only one trace
            trace.data = baseline.data(j:j,:);
            
            %removes QC Spike
            iQCSpikeBegin = find(trace.data == min(trace.data));
            iQCSpikeEnd = find(trace.data == max(trace.data));
            trace.data = trace.data(1:iQCSpikeBegin-5);
            
            %% Set major overall parameters
            params.dt = 1/10000;
            params.savename = [string '-proc.mat'];
            params.full_save_string = fullfile(savePath1, params.savename);
            params.nACQ = nACQ;
            
            %% Set parameters specific to ePSC detection
            %Event direction 
            paramse.event_sign = -1;
            
            %In std of filtered trace (uses the zscore, so any spike above
                %this number std gets marked as an event. This really depends on the noise of the data   
            paramse.init_method.threshold = 2.25; 
                                                   
            %Detects events only once in a this many frame span.
            paramse.init_method.min_interval = 250;
            
            %Time in seconds per sample (your sampling rate)
            paramse.dt = 1/10000; 
            
            paramse.init_method.ar_noise_params.sigma_sq = 3;
            %Last number here is your noise parameter-- important when data is really low or really high noise
                %See powerpoint for more information about how to adjust this parameter accordingly
            paramse.init_method.ar_noise_params.phi = [1.000000000000000, 1.0, -0.15]; 
            
            paramse.traces_file = string;
            params.paramse = paramse;
            
            %% Detect ePSCs
            
            %Local path to template file for ePSC detection
            paramse.init_method.template_file_e = fullfile('/Volumes', 'Carly Rose', '2 - Code', '1 - MATLAB', 'CRW-PSC-Detection-master', 'template', 'epsc-template.mat');
            
            %Loads template and prepares data for wiener filtering
            load_struct = load(paramse.init_method.template_file_e);
            template_e = load_struct.template;
            tracein_e = paramse.event_sign*trace.data;
            tracein_e = tracein_e - min(tracein_e);

            %Run wiener filter using template matching and deconvolution to
                %detect event times and calculate amplitude.
            nfft = length(tracein_e) + length(template_e);
            [trace.filtered_trace_e, trace.event_times_e, trace.event_amp_e] = wiener_filter_biphasic(tracein_e*paramse.event_sign, paramse, template_e, nfft);
            
            %Creates matrix of event times and amplitudes to combine with
                %iPSC results
            time_amp_e = [trace.event_times_e; trace.event_amp_e];
           
            %% Set parameters for iPSC Detection
            %Event direction
            paramsi.event_sign = 1;
            
            %In std of filtered trace (uses the zscore, so any spike above
               %this number std gets marked as an event. This really depends on the noise of the data 
            paramsi.init_method.threshold = 2.28;     
            
            %Detects events only once in a this many frame span.
            paramsi.init_method.min_interval = 175; 
            
            %Time in seconds per sample (your sampling rate)
            paramsi.dt = 1/10000; 
            
            paramsi.init_method.ar_noise_params.sigma_sq = 3;
            %Last number here is your noise parameter-- important when data is really low or really high noise
                %See powerpoint for more information about how to adjust this parameter accordingly
            paramsi.init_method.ar_noise_params.phi = [1.000000000000000, 1.0, -0.28];
            
            %update params to reflect current trace being analyzed 
            paramsi.traces_file = string;
            params.paramsi = paramsi;

            %% Detect iPSCs
           
            %Local path to template file for iPSC detection
            paramsi.init_method.template_file_i = fullfile('/Volumes', 'Carly Rose', '2 - Code', '1 - MATLAB', 'CRW-PSC-Detection-master', 'template', 'ipsc-template.mat');
            
            %Loads template and prepares data for wiener filtering
            load_struct = load(paramsi.init_method.template_file_i);
            template_i = load_struct.template;
            tracein_i = paramsi.event_sign*trace.data;
            tracein_i = tracein_i - min(tracein_i);
            
            %Run wiener filter using template matching and deconvolution to
                %detect event times and calculate amplitude.
            nfft = length(tracein_i) + length(template_i);
            [trace.filtered_trace_i, trace.event_times_i, trace.event_amp_i] = wiener_filter_biphasic(tracein_i*paramsi.event_sign, paramsi, template_i, nfft);
            
            %Creates matrix of event times and amplitudes to combine with
                %ePSC results
            time_amp_i = [trace.event_times_i; trace.event_amp_i];
            
            %% Create full matrix of event times and amplitudes.
            
            %Combines event times and amplitudes for iPSCs and ePSCs
            time_amp_final = [time_amp_e time_amp_i];
            
            %Sorts event times and amplitudes by time.
            time_amp_final = sortrows(time_amp_final', 1)';
          
            %store full parameters in output variable
            trace.params = params;
            
            %isolate event times/amplitudes and store to output variable
            trace.event_times = time_amp_final(1:1,:);
            trace.event_amp = time_amp_final(2:2,:);
            
            %create ISI structure 
            if isempty(trace.event_times) == 1
                trace.ISIs = [];
            else
                trace.ISIs = [trace.event_times(1) diff(trace.event_times)];
            end 
            
            %Create spike train for later analysis
            trace.SpikeTrain = zeros([1 length(trace.data)]);
            
            for times = trace.event_times;
                trace.SpikeTrain(times) = 1;
            end

            trace.name = string;
            
        %% Save time~ 
            %save output variable to unique variable name
            assignin('base', string, trace);
            
            %save to output folder
            save(params.full_save_string, string);
            
            %clear extraneous variables before loop restarts
            clear(string);
            clear params paramse paramsi 
        end
    end
    
%% Turn back on Warnings        
warning('on', 'signal:findpeaks:largeMinPeakHeight');
warning('on', 'MATLAB:unknownObjectNowStruct');
warning('on', 'MATLAB:MKDIR:DirectoryExists');