%% Paired Pulse Ratio Analysis Pipeline%%

%Written by CRW, 27 May 2019
    
    %line 130 and 132 is really the important thing to change if your init
        %times seem wrong. I suggest plotting a refline at several
        %different points on a zscore of the filtered trace to make a more
        %educated decision for threshold. Decrease init_window if you have
        %spikes that are closer together.
    %change lines 26, 29, 118 and 124 to reflect local path.
    
function PSC_Process_PPR(PSCTableDate, datedfolder);

%% Turn off dumb warning for loading "wave" struct and directory
warning('off', 'MATLAB:unknownObjectNowStruct');
warning('off', 'MATLAB:MKDIR:DirectoryExists');

%% Initialize %%
    [nrows ncolumns] = size(PSCTableDate);

%Makes input path given date information
    prepath = fullfile('//Volumes', 'Neurobio', 'MICROSCOPE', 'Kevin', '3-Experiments', '4-SliceEphys', '7-Paired Pulse', 'Preprocessed Data', datedfolder);

%Makes save path given date information
    savePath = fullfile('//Volumes', 'Neurobio', 'MICROSCOPE', 'Kevin', '3-Experiments', '4-SliceEphys', '7-Paired Pulse', '2-Output', strcat(datedfolder, '_output'));

%User enters which Epochs to run
    Epochs = input(['You have ', num2str(nrows) ' epochs. Input epochs to run in matrix form ']);
    disp(' ');
%% Organize by Epoch 
    %input for acquisitions and creation of filename directories
    for i = Epochs;
        disp(['Beginning epoch #' num2str(PSCTableDate{i,13})]);
        disp('------------------------------');
        acqsweeps = PSCTableDate{i, 11}:PSCTableDate{i,12};
        nACQ = length(acqsweeps);
        experiment = PSCTableDate{i, 5};
        mouseID = PSCTableDate{i, 2};
        celll = num2str(PSCTableDate{i, 3});
        epochh = num2str(PSCTableDate{i,13});
        names = cell(1, nACQ); %preallocate size of names
        savePath1 = fullfile(savePath, strcat('cell_', num2str(PSCTableDate{i, 3})), strcat('epoch_', num2str(PSCTableDate{i,13})));
        mkdir(savePath1);
        %prepath = fullfile(prepath, strcat('cell_', celll), strcat('epoch_', epochh));
        
        %lists file names for every sweep in Epoch
        for nameindex = 1:nACQ
            names{nameindex} = strcat('AD0_', num2str(acqsweeps(nameindex)), '.mat');
        end

%% Set Parameters for Each Epoch and load files
        figure
        raw_concatenated_traces = [];
        for file = 1:nACQ
            %reset base_params and params for each wave file
            params = [];
            
            %initialize for get_params function
            params.traces_filename = fullfile(prepath, strcat('cell_', celll), strcat('epoch_', epochh), names{file});
            params.savepath = savePath1;
            
            try 
                load(params.traces_filename);
            catch
                fprintf(strcat(names{file}, ' failed to load... skipping'));
                disp(' ');
                continue;
            end 
              
            name = eval(names{file}(1:end-4));

            
            %isolate peristimulus windows and assign event sign.
            name.experiment = experiment;

            if name.experiment == -20
                name.peristim1 = name.data(10010:10500);
                name.peristim2 = name.data(10510:11010);
                params.event_sign = -1;
                subplot(2,1,1)
                plot(name.peristim1);
                title('After First Stimulus')
                hold on;
                subplot(2,1,2)
                plot(name.peristim2);
                title('After Second Stimulus');
                xlim([0 500]);
                hold on;
            elseif name.experiment == -70
                name.peristim1 = name.data(10010:10700);
                params.event_sign = -1;
                plot(name.peristim1)
                title('After Stimulus');
                hold on;
            elseif name.experiment == 40
                name.peristim1 = name.data(10010:10700);
                params.event_sign = 1;
                plot(name.peristim1);
                title('After Stimulus');
                hold on;
            end
            
           
            %define approproate template file
            if params.event_sign == -1
                params.init_method.template_file = fullfile('//Volumes', 'Carly Rose', '2 - Code', '1 - MATLAB', 'CRW-PSC-Detection-master', 'template', 'epsc-template.mat');
            else
                params.init_method.template_file = fullfile('//Volumes', 'Carly Rose', '2 - Code', '1 - MATLAB', 'CRW-PSC-Detection-master', 'template', 'ipsc-template.mat');
            end 
           
            %Initialize Wiener filter parameters
            
            params.init_method.threshold = 2.1; %in std of filtered trace (uses the zscore, so any spike above 2.35 std gets marked as an
                                                   %event. This really depends on the noise ofthe data
            params.init_method.min_interval = 100; %detects events only once in a 100 frame span.
            params.dt = 1/10000; %time in seconds per sample
            params.animal_information = PSCTableDate{i, 4};
            params.init_method.ar_noise_params.sigma_sq = 3;
            params.init_method.ar_noise_params.phi = [1.000000000000000, 1.0, -0.15];
            
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
            
            baseline = mean(name.data(1:10000));
            
            if experiment == -20
                trace1 = name.peristim1(15:400) * params.event_sign;
                trace1 = trace1 - min(trace1);
                nfft1 = length(trace1) + length(template);
                [filtered_trace1, event_times1] = wiener_filter_PPR(trace1*params.event_sign, params, template, nfft1);
                event_times1 = event_times1+15;
                event_times1(event_times1>300)=[];
                name.event_times1 = event_times1;
                name.filtered_trace1 = filtered_trace1;
                name.tau1 = [];
 
                name.tau2 = [];
                trace2 = name.peristim2(15:330) * params.event_sign;
                trace2 = trace2 - min(trace2);
                nfft2 = length(trace2) + length(template);
                [filtered_trace2, event_times2] = wiener_filter_PPR(trace2*params.event_sign, params, template, nfft2);
                event_times2 = event_times2+15;
                event_times2(event_times2>300)=[];
                name.event_times2 = event_times2;
                name.filtered_trace2 = filtered_trace2;

                if ~isempty(name.event_times2) == 1;
                    name.event_times2 = name.event_times2(1);
                    name.event_amp2 = (min(name.peristim2(25:name.event_times2+50))-baseline)*params.event_sign;
                    [peak, loc] = min(name.peristim2(25:name.event_times2+50));
                    if loc < name.event_times2 
                        name.event_times2 = loc;
                    end 
                    fittingTracee = [name.peristim2(loc:loc+200)-baseline];
                    xx = linspace(0, length(fittingTracee)/10, length(fittingTracee));
                    ff = fit(xx', fittingTracee', 'exp2', 'Robust', 'LAR', 'Upper', [-10 0 -5 0]);
                    name.tau2 = [name.tau2; ff.b ff.d];
                else
                    name.tau2 = [name.tau2];
                    name.event_amp2 = [];
                end 

                clear xx ff peak loc fittingTracee 
                
                if ~isempty(name.event_times1) == 1;
                    name.event_times1 = name.event_times1(1);
                    name.event_amp1 = (min(name.peristim1(25:name.event_times1+50))-baseline)*params.event_sign;
                    [peak, loc] = min(name.peristim1(25:name.event_times1+50));
                    if loc < name.event_times1
                        name.event_times1 = loc;
                    end
                    fittingTracee = [name.peristim1(loc:loc+200)-baseline];
                    xx = linspace(0, length(fittingTracee)/10, length(fittingTracee));
                    ff = fit(xx', fittingTracee', 'exp2', 'Robust', 'LAR', 'Upper', [-10 0 -5 0]);
                    name.tau1 = [name.tau1; ff.b ff.d];
                else
                    name.tau1 = [name.tau1];
                    name.event_amp1 = [];
                end
                
                clear xx ff peak loc fittingTracee 
                
                if ~isempty(name.event_times1) == 1 && ~isempty(name.event_times2)==1
                    name.PPR = name.event_amp2/name.event_amp1;
                else
                    name.PPR = [];
                end
            end
            

            %store initialization results in structure and save
            params.epoch = num2str(PSCTableDate{i,13});
            params.cell = num2str(PSCTableDate{i, 3});
            params.date = date;
            params.mouseID = PSCTableDate{i, 2};
            name.name = params.traces_file(1:end-4);
            name.params = params;
        
            if experiment == -70
                trace1 = name.peristim1(15:end) * params.event_sign;
                trace1 = trace1 - min(trace1);
                nfft1 = length(trace1) + length(template);
                [filtered_trace1, event_times1] = wiener_filter_PPR(trace1*params.event_sign, params, template, nfft1);
                event_times1 = event_times1+15;
                event_times1(event_times1>300)=[];
                name.event_times1 = event_times1;
                name.filtered_trace1 = filtered_trace1;
                name.tau1 = [];
                
                if ~isempty(name.event_times1) == 1
                    name.event_times1 = name.event_times1(1);
                    name.event_amp1 = (min(name.peristim1(25:name.event_times1+50))-baseline)*params.event_sign;
                    [peak, loc] = min(name.peristim1(25:name.event_times1+50));
                    fittingTracee = [name.peristim1(loc:loc+300)-baseline];
                    xx = linspace(0, length(fittingTracee)/10, length(fittingTracee));
                    ff = fit(xx', fittingTracee', 'exp2', 'Robust', 'LAR', 'Upper', [-10 0 -5 0]);
                    name.tau1 = [name.tau1; ff.b ff.d];
                else
                    name.tau1 = [name.tau1];
                    name.event_amp1 = [];
                end 
            end 
            
            if experiment == 40
                name.params.init_method.ar_noise_params.phi = [1.000000000000000, 1.0, -0.10];
                trace1 = name.peristim1 * params.event_sign;
                trace1 = trace1 - min(trace1);
                nfft1 = length(trace1) + length(template);
                [filtered_trace1, event_times1] = wiener_filter_PPR(trace1*params.event_sign, params, template, nfft1);
                event_times1 = 400;
                name.event_times1 = event_times1;
                name.filtered_trace1 = filtered_trace1;
                name.tau1 = [];
                
                if ~isempty(name.event_times1) == 1
                    name.event_times1 = name.event_times1(1);
                    name.event_amp1 = (name.peristim1(400)-baseline)*params.event_sign;
                    %[peak, loc] = min(name.peristim1(25:name.event_times1+50));
                    %fittingTracee = [name.peristim1(loc:loc+300)-baseline];
                    %xx = linspace(0, length(fittingTracee)/10, length(fittingTracee));
                    %ff = fit(xx', fittingTracee', 'exp2', 'Robust', 'LAR', 'Upper', [-10 0 -5 0]);
                    %name.tau1 = [name.tau1; ff.b ff.d];
                else
                    name.event_amp1 = [];
                end 
            end 
                
            
            clear  peak loc fittingTracee xx ff 
            
            
            assignin('base', params.traces_file(1:end-4), name);
            
            raw_concatenated_traces = [raw_concatenated_traces name];

            save(params.full_save_string, params.traces_file(1:end-4));
            
%% Plot things and Move On!

            %displays to user based on current stage of analysis
            if file ~= nACQ
                if file == nACQ-1
                    disp('On to the next! Only one more to go!')
                    disp('------------------------------')
                else
                    disp(['**On to the next! ' num2str(nACQ-file) ' traces left for epoch ' num2str(PSCTableDate{i,13}) '***']);
                    disp('------------------------------');
                end
            else
                disp(['Done with epoch #' num2str(PSCTableDate{i,13})]);
            end
           
            %clear all variables from this trace to save memory
            clear(params.traces_file(1:end-4))
            epoch = num2str(params.epoch);
            clear name params results trace times event_times1 event_amp1 event_times2 event_amp2 filtered_trace1 filtered_trace2 QC 
        
        end
        
        %save overlaid traces
        print(fullfile(savePath1, strcat(mouseID, '_cell', celll, '_epoch', epoch, '_Overlay_Traces')), '-dpdf', '-fillpage', '-r1000');
        save(fullfile(savePath1, strcat(mouseID, '_cell', celll,  '_epoch', epoch,  'raw_concatenated_traces.mat')), 'raw_concatenated_traces');
      
        
        close all
        
        disp('Running Quality Control Check');
        QCs = zeros(length(raw_concatenated_traces), length(raw_concatenated_traces(1).QC));
        for QCCheck = 1:length(raw_concatenated_traces)
            QCs(QCCheck,:) = raw_concatenated_traces(QCCheck).QC;
        end
        
        for l = 1:size(QCs,1)
            RCtrace=QCs(l,1001:end-1003); %isolate only the beginning and downward peak of the QC check
            steadyRC = mean(RCtrace((end-40):end));   %find the mean steady state current of 4ms before upward peak
            dV = 5; % voltage step in mV
            [peakRC,peakloc] = min(RCtrace); %peak of downward spike 
            baselineRC(l) = mean(QCs(l,1:999));  %find the mean baseline current
        
            %Rs(l) = dV/abs(steadyRC(l)-maxpeak).*1000;
    
            fittingTrace = [RCtrace((peakloc):(end-60))-steadyRC];  %fit the trace from the peak
            x=linspace(0,length(fittingTrace)/10,length(fittingTrace)); %in ms
            f=fit(x',fittingTrace','exp2','Robust','LAR','Upper',[-10 0 -5 0]); %fit two exponentials
    
            % extrapolation back to onset of current change
            fittingTrace2=[RCtrace((peakloc):(end-65))-steadyRC];
            [~,locOnset]=max(abs(diff(fittingTrace2)));   
            x2=x-5.*0.1002;
            fittedVal =feval(f,x2);
            extrapolPkRC = feval(f,(x(1)-0.1002.*(5-locOnset)))+steadyRC;
            %plot(f,x2,fittingTrace2)
    
            tauParams(l,:)=[f.b, f.d];  % save both taus (in ms)   
            Raccess(l) = dV./abs(extrapolPkRC-baselineRC(l)).*1000;   % Access resistance (in MOhm)
            Rm(l) = dV./abs(steadyRC-baselineRC(l)).*1000-Raccess(l);   % membrane resistance (in MOhm)
            hold on;
        end

        tau = max(-tauParams'); % slowest component is the membrane time constant
        tauCap = min(-tauParams');
   
    
    %Calculate capacitance per trace
        Cm = tau./Rm*1000;   %membrane capacitance (in picoFarad)

    %Plot everything
    figure
        for m = 1:size(QCs,1)
            subplot(3,2,1);
            plot(QCs(m, 900:1200));
            ylabel('QC Pulse');
            xlim([0 300]);
            ylim([(min(QCs(m, 900:1200)-300)) (max(QCs(m, 900:1200)+300))]);
            hold on
        end
        subplot(3,2,2);
        plot(baselineRC, 'ko');
        ylabel('Baseline');
        xlim([0 size(QCs,1)])
        subplot(3,2,3);
        plot(Raccess,'ko');
        xlim([0 size(QCs,1)])
        ylabel('Raccess (MOhm)')
        subplot(3,2,4);
        plot(Rm,'bo');
        xlim([0 size(QCs,1)])
        ylabel('Rm (MOhm)')
        subplot(3,2,5);
        plot(tau,'ko');
        xlim([0 size(QCs,1)])
        ylabel('{\tau}_m (ms)')
        subplot(3,2,6);
        plot(Cm,'ko');
        xlim([0 size(QCs,1)])
        ylabel('Cm (pF)')
        xlabel('Sweep #')

    %Make it look pretty
        sgtitle(strcat(mouseID, ' cell ', celll, ', epoch ', num2str(i)));
        set(gcf, 'Color', 'w');

    %Adjust and save
        print(fullfile(savePath1, strcat(mouseID, '_cell', celll, '_epoch', epochh, '_QualityControl')), '-dpdf', '-fillpage', '-r1000');

    %Close everything
        close all
    
        
        clear names acsweeps nACQ savePath1 RCtrace steadyRC 
        clear peakRC peakloc baselineRC fittingTrace fittingTrace2 x f x2 fittedVal extrapolPkRC
        clear tauParams Raccess Rm tau tauCap Cm QCs 
    end

%% End display and turn back on warnings
disp('Aaaaand we''re out!');
disp('~~Analysis is complete!~~');            
warning('on', 'MATLAB:unknownObjectNowStruct');
warning('on', 'MATLAB:MKDIR:DirectoryExists');
end