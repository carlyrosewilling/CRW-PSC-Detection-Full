%% Pipeline for Bayesian Algorithm Without Sampler: Includes Read Excel Functionality%%

%Written by CRW, 4 Oct 2018
    %last updated: 8 Jan 2019
    
    %line 130 and 132 is really the important thing to change if your init
        %times seem wrong. I suggest plotting a refline at several
        %different points on a zscore of the filtered trace to make a more
        %educated decision for threshold. Decrease init_window if you have
        %spikes that are closer together.
    %change lines 26, 29, 118 and 124 to reflect local path.

%% Turn off dumb warning for loading "wave" struct and directory
warning('off', 'MATLAB:unknownObjectNowStruct');
warning('off', 'MATLAB:MKDIR:DirectoryExists');

%% Initialize %%

%User selects exact file where data is stored
    %prepath = fullfile('/Users', 'carlyrosewilling', 'Desktop', 'Data', 'Minis', 'KM092618/');
    %User inputs date of acquisition
    date = input('Input date of recording (i.e. 01/06/2019): ', 's');
    
%Makes input path given date information
    datedfolder = strcat('KM', date(1:2), date(4:5), date(9:end));
    prepath = fullfile('//Volumes', 'Neurobio', 'MICROSCOPE', 'Kevin', '3-Experiments', '4-SliceEphys', '6-WT Normal Development', '1-Raw Data', datedfolder);

%Makes save path given date information
    savePath = fullfile('//Volumes', 'Neurobio', 'MICROSCOPE', 'Kevin', '3-Experiments', '4-SliceEphys', '6-WT Normal Development', '2-Output', strcat(datedfolder, '_output'));

%User enters date of aquisition for parsing Excel Sheet
    tabledate = str2num(strcat(date(7:end), date(1:2), date(4:5)));
    
%Loads overview Excel sheet 
    %User selects Excel file
    PSC_LoadExcel
    
    PSCTableDate = {};
    [row column] = size(PSCTableRaw);
    for j = 2:row
        if PSCTableRaw{j,1} == tabledate
            row = horzcat(PSCTableRaw(j,:));
            PSCTableDate = vertcat(PSCTableDate, row);
        else
            PSCTableDate = PSCTableDate;
        end
    end
    [nrows ncolumns] = size(PSCTableDate);
    
%User enters which Epochs to run
    Epochs = input(['You have ', num2str(nrows) ' epochs. Input epochs to run in matrix form ']);
    %Epochs = [1];
    disp(' ');
%% Organize by Epoch 
    %input for acquisitions and creation of filename directories
    for i = Epochs;
        disp(['Beginning epoch #' num2str(PSCTableDate{i,17})]);
        disp('------------------------------');
        acqsweeps = PSCTableDate{i, 15}:PSCTableDate{i,16};
        nACQ = length(acqsweeps);
        event_sign = PSCTableDate{i, 14};
        names = cell(1, nACQ); %preallocate size of names
        savePath1 = fullfile(savePath, strcat('cell_', num2str(PSCTableDate{i, 3})), strcat('epoch_', num2str(PSCTableDate{i,17})));
        mkdir(savePath1);
        
        %lists file names for every sweep in Epoch
        for nameindex = 1:nACQ
            names{nameindex} = strcat('AD0_', num2str(acqsweeps(nameindex)), '.mat');
        end

%% Set Parameters for Each Epoch and load files

        for file = 1:nACQ
            %reset base_params and params for each wave file
            params = [];
            
            %initialize for get_params function
            params.traces_filename = strcat(prepath, '/', names{file});
            params.savepath = savePath1;
            

            try 
                load(params.traces_filename);
            catch
                fprintf(strcat(names{file}, ' failed to load... skipping'));
                disp(' ');
                continue;
            end 
              
            name = eval(names{file}(1:end-4));
           
            %Use only data that does not contain quality control spike,
                %saves QC spike
                %also subtracting last 20 and first 25 frames due to non-event spikes in data here
                %also eliminates 1500 frames after QC spike to allow for
                %time for steady baseline to form
            %name.data = name.data(25:end);
            %iQCSpikeBegin = find(name.data == min(name.data));
            %iQCSpikeEnd = find(name.data == max(name.data));
            %if iQCSpikeEnd > 10000
                %traces = name.data(1:iQCSpikeBegin-200);
            %elseif iQCSpikeEnd < 10000
                %traces = name.data(iQCSpikeEnd+1500:end-20);
            %else 
                %traces = name.data(iQCSpikeEnd+200:end-20);
            %end
            %QC = name.data(iQCSpikeBegin-2:iQCSpikeEnd+250);
            traces = name.data;
            
            %Use event sign gathered from excel in line 52
            params.event_sign = event_sign;
            
            %Get a_min/max using baseline data and event sign. (these don't
                %really matter anymore since we aren't using sampler)
            if params.event_sign == -1
                %params.b_min = mean(traces) - std(traces);
                %params.b_max = max(traces);
                %params.a_min = params.b_min - 50;
                %params.a_max = mean(traces) + 1.5 * std(traces); 
                params.init_method.template_file = fullfile('//Volumes', 'Carly Rose', '2 - Code', '1 - MATLAB', 'CRW-PSC-Detection-master', 'template', 'epsc-template.mat');
            else
                %params.b_min = min(traces);
                %params.b_max = mean(traces) + std(traces);
                %params.a_min = mean(traces) - 1.5 * std(traces);
                %params.a_max = params.b_max + 50; 
                params.init_method.template_file = fullfile('//Volumes', 'Carly Rose', '2 - Code', '1 - MATLAB', 'CRW-PSC-Detection-master', 'template', 'ipsc-template.mat');
            end 
           
            %These are the wiener filter parameters I've found work best,
                %may be changed depending on how data is structured?
            
            params.init_method.threshold = 2.375; %in std of filtered trace (uses the zscore, so any spike above 2.35 std gets marked as an
                                                   %event. This really depends on the noise ofthe data
            params.init_method.min_interval = 100; %detects events only once in a 170 frame span.
            params.dt = 1/10000; %time in seconds per sample
            params.location = PSCTableDate{i, 4};
            params.init_method.ar_noise_params.sigma_sq = 3;
            params.init_method.ar_noise_params.phi = [1.000000000000000, 1.0, -0.3];
            %update params to reflect current trace being analyzed 
            params.traces_file = names{file};
            params.savename = [params.traces_file(1:end-4) '-proc.mat'];
            params.full_save_string = fullfile(params.savepath, params.savename);
            params.number_of_trace = file;
            params.nACQ = nACQ;

%% Run Wiener Filter For Trace
            disp(['Now running inference on ' params.traces_file(1:end-4) ' data...']);

            %load template for init method
            load_struct = load(params.init_method.template_file);
            template = load_struct.template;
            
            raw_trace = traces; %preserve raw trace
            trace = params.event_sign*traces;
            trace = trace - min(trace);

            %run wiener filter to find init event times
            nfft = length(trace) + length(template);
        
            [filtered_trace, event_times, event_amp] = wiener_filter(trace*params.event_sign, params, template, nfft);

            %store initialization results in structure and save
            params.epoch = num2str(PSCTableDate{i,17});
            params.cell = num2str(PSCTableDate{i, 3});
            params.date = date;
            params.mouseID = PSCTableDate{i, 2};
            name.name = params.traces_file(1:end-4);
            name.params = params;
            name.event_times = event_times;
            name.event_amp = event_amp;
            if isempty(name.event_times) == 1
                name.ISIs = [];
            else
                name.ISIs = [event_times(1) diff(event_times)];
            end 
            
            name.SpikeTrain = zeros([1 length(raw_trace)]);
            
            for times = event_times;
                name.SpikeTrain(times) = 1;
            end

            name.filtered_trace = filtered_trace;
            name.raw_trace = raw_trace;
            %name.QC = QC;
            
            assignin('base', params.traces_file(1:end-4), name);
            
            save(params.full_save_string, params.traces_file(1:end-4));

%% Plot things and Move On!

            
            %plot some shit! 
            %disp(' ');
            %disp('Plotting Traces and Events');
            
            %plot_init_events_over_trace(raw_trace, results, params)
            %Uncomment if you want a plot for every trace
            
            %saveas(gcf, strcat(params.savepath, params.savename(1:end-4),
            %'.fig')); uncomment this if you want to save a figure for
            %every sweep.
            
            % ...and we're out.
            %disp('~~This trace is done!~~');
            
            %displays to user based on current stage of analysis
            if file ~= nACQ
                if file == nACQ-1
                    disp('On to the next! Only one more to go!')
                    disp('------------------------------')
                else
                    disp(['**On to the next! ' num2str(nACQ-file) ' traces left for epoch ' num2str(PSCTableDate{i,17}) '***']);
                    disp('------------------------------');
                end
            else
                disp(['Done with epoch #' num2str(PSCTableDate{i,17})]);
            end
           
            %clear all variables from this trace to save memory
            clear(params.traces_file(1:end-4))
            clear name base_params params results traces trace times event_times event_amp filtered_trace raw_trace iQCSpikeBegin iQCSpikeEnd QC 
        
        end
    clear names acsweeps nACQ savePath1
    end

%% End display and turn back on warnings
disp('Aaaaand we''re out!');
disp('~~Analysis is complete!~~');            
warning('on', 'MATLAB:unknownObjectNowStruct');
warning('on', 'MATLAB:MKDIR:DirectoryExists');
