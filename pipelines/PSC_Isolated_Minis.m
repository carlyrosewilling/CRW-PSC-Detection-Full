%% Pipeline for Processing of Isolated Mini Data%%

%Written by CRW, 27 May 2019
    
    %line 130 and 132 is really the important thing to change if your init
        %times seem wrong. I suggest plotting a refline at several
        %different points on a zscore of the filtered trace to make a more
        %educated decision for threshold. Decrease init_window if you have
        %spikes that are closer together.
    %change lines 26, 29, 118 and 124 to reflect local path.

function PSC_Isolated_Minis(PSCTableDate, datedfolder);
    
    %% Turn off dumb warning for loading "wave" struct and directory
warning('off', 'MATLAB:unknownObjectNowStruct');
warning('off', 'MATLAB:MKDIR:DirectoryExists');

%% Initialize %%
     [nrows ncolumns] = size(PSCTableDate);
%Makes input path given date information
    prepath = fullfile('//Volumes', 'Neurobio', 'MICROSCOPE', 'Kevin', '3-Experiments', '4-SliceEphys', '5-Isolated Mice', 'Preprocessed Data', datedfolder);

%Makes save path given date information
    savePath = fullfile('//Volumes', 'Neurobio', 'MICROSCOPE', 'Kevin', '3-Experiments', '4-SliceEphys', '5-Isolated Mice', '2-Output', strcat(datedfolder, '_output'));
    
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
        celll = num2str(PSCTableDate{i, 3});
        epochh = num2str(PSCTableDate{i, 17});
        
        %lists file names for every sweep in Epoch
        for nameindex = 1:nACQ
            names{nameindex} = strcat('AD0_', num2str(acqsweeps(nameindex)), '.mat');
        end

%% Set Parameters for Each Epoch and load files

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
           
            
            %Use event sign gathered from excel in line 52
            params.event_sign = event_sign;
            
            %Get a_min/max using baseline data and event sign. (these don't
                %really matter anymore since we aren't using sampler)
            if params.event_sign == -1 
                params.init_method.template_file = fullfile('//Volumes', 'Carly Rose', '2 - Code', '1 - MATLAB', 'CRW-PSC-Detection-master', 'template', 'epsc-template.mat');
            else
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
            params.init_method.ar_noise_params.phi = [1.000000000000000, 1.0, -0.25];
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
            name.QC = QC;
            
            raw_concatenated_traces = [raw_concatenated_traces name];
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
                    disp(['**On to the next! ' num2str(nACQ-file) ' traces left for epoch ' epochh '***']);
                    disp('------------------------------');
                end
            else
                disp(['Done with epoch #' epochh]);
            end
           
            %clear all variables from this trace to save memory
            clear(params.traces_file(1:end-4))
            clear name base_params params results traces trace times event_times event_amp filtered_trace raw_trace iQCSpikeBegin iQCSpikeEnd QC 
        
        end
        disp('Running Quality Control Check');
                QCs = zeros(length(raw_concatenated_traces), length(raw_concatenated_traces(1).QC));
                for QCCheck = 1:length(raw_concatenated_traces)
                    QCs(QCCheck,:) = raw_concatenated_traces(QCCheck).QC;
                end

                for l = 1:size(QCs,1)
                    RCtrace=QCs(l,1001:end-803); %isolate only the beginning and downward peak of the QC check
                    steadyRC = mean(RCtrace((end-40):end));   %find the mean steady state current of 4ms before upward peak
                    dV = 5; % voltage step in mV
                    [peakRC,peakloc] = min(RCtrace); %peak of downward spike 
                    baselineRC(l) = mean(QCs(l,1:999));  %find the mean baseline current

                    %Rs(l) = dV/abs(steadyRC(l)-maxpeak).*1000;

                    fittingTrace = [RCtrace((peakloc):(end-60))-steadyRC];  %fit the trace from the peak
                    x=linspace(0,length(fittingTrace)/10,length(fittingTrace)); %in ms

                    try
                        f=fit(x',fittingTrace','exp2','Robust','LAR','Upper',[-10 0 -5 0]); %fit two exponentials
                    catch
                        string1 = " " + raw_concatenated_traces(l).filename(1:end-4);
                        fprintf(strcat('Could not create fit for trace ', string1, '. QC pulse may be corrupted. Consider rejecting.'));
                        disp(' ');
                        tauParams(l,:) = [NaN, NaN];
                        Raccess(l) = NaN;
                        Rm(l) = NaN;
                        continue 
                    end 

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
                sgtitle(strcat(mouseID, ' cell ', celll, ', epoch ', epochh));
                set(gcf, 'Color', 'w');

            %Adjust and save
                print(fullfile(savePath1, strcat(mouseID, '_cell', celll, '_epoch', epochh, '_QualityControl')), '-dpdf', '-fillpage', '-r1000');

            %Close everything
                close all
    clear names acsweeps nACQ savePath1
    end

%% End display and turn back on warnings
disp('Aaaaand we''re out!');
disp('~~Analysis is complete!~~');            
warning('on', 'MATLAB:unknownObjectNowStruct');
warning('on', 'MATLAB:MKDIR:DirectoryExists');
