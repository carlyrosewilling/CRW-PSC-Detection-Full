%% EStim Analysis Pipeline%%

%Written by CRW, 17 June 2019
    
function PSC_Process_EStim(PSCTableDate, datedfolder);

%% Turn off dumb warning for loading "wave" struct and directory
warning('off', 'MATLAB:unknownObjectNowStruct');
warning('off', 'MATLAB:MKDIR:DirectoryExists');
warning('off', 'MATLAB:centerOfMass:neg');

%% Initialize %%

%User inputs date of acquisition
    prompt = {'Enter 1st stim pulse location:', 'Enter 2nd stim pulse location:'};
    dlgtitle = 'Inputs';
    dims = [1 75];
    definput = {'10005', '10505'};
    answer = inputdlg(prompt, dlgtitle, dims, definput);

    stimlocation = answer{1};
    stimlocation = str2num(stimlocation);
    PPRsecondstim = answer{2};
    PPRsecondstim = str2num(PPRsecondstim);
    [nrows ncolumns] = size(PSCTableDate);

%Makes input path given date information
    prepath = fullfile('//Volumes', 'Neurobio', 'MICROSCOPE', 'Kevin', '3-Experiments', '4-SliceEphys', '10-EStim', 'Preprocessed Data', datedfolder);

%Makes save path given date information
    savePath = fullfile('//Volumes', 'Neurobio', 'MICROSCOPE', 'Kevin', '3-Experiments', '4-SliceEphys', '10-EStim', '2-Output', strcat(datedfolder, '_output'));

%User enters which Epochs to run
    Epochs = input(['You have ', num2str(nrows) ' epochs. Input epochs to run in matrix form ']);
    disp(' ');
    
    outputpath = fullfile('//Volumes', 'Neurobio', 'MICROSCOPE', 'Kevin', '3-Experiments', '4-SliceEphys', '10-EStim', '2-Output');
    
%% Organize by Epoch 
    %input for acquisitions and creation of filename directories
    for i = Epochs;
        disp(['Beginning epoch #' num2str(PSCTableDate{i,18})]);
        disp('------------------------------');
        acqsweeps = PSCTableDate{i, 16}:PSCTableDate{i,17};
        nACQ = length(acqsweeps);
        experiment = PSCTableDate{i, 13};
        expdate = PSCTableDate{i,1};
        mouseID = PSCTableDate{i, 2};
        celll = num2str(PSCTableDate{i, 3});
        location = PSCTableDate{i, 4};
        epochh = num2str(PSCTableDate{i,18});
        magnitude = PSCTableDate{i,14};
        holding_current = PSCTableDate{i,15};
        names = cell(1, nACQ); %preallocate size of names
        savePath1 = fullfile(savePath, strcat('cell_', num2str(PSCTableDate{i, 3})), strcat('epoch_', num2str(PSCTableDate{i,18})));
        mkdir(savePath1);
        
        [~, ~, AllTracesTable] = xlsread(fullfile(outputpath, 'EStim_AllTraces2.xlsx'));
        [~, ~, AverageTraces] = xlsread(fullfile(outputpath, 'EStim_AverageTraces2.xlsx'));
        %prepath = fullfile(prepath, strcat('cell_', celll), strcat('epoch_', epochh));
        
        %lists file names for every sweep in Epoch
        for nameindex = 1:nACQ
            names{nameindex} = strcat('AD0_', num2str(acqsweeps(nameindex)), '.mat');
        end

%% Set Parameters for Each Epoch and load files
        figure
        raw_concatenated_traces = [];

%% CHR2 max                
            if isequal(experiment, 'max') == 1
                figure
                for file = 1:nACQ
                    %reset base_params and params for each wave file
                    params = [];
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

                    if holding_current == -70
                        params.event_sign = -1;
                         params.init_method.template_file = fullfile('//Volumes', 'Carly Rose', '2 - Code', '1 - MATLAB', 'CRW-PSC-Detection-master', 'template', 'epsc-template.mat');
                    else
                        params.event_sign = 1;
                         params.init_method.template_file = fullfile('//Volumes', 'Carly Rose', '2 - Code', '1 - MATLAB', 'CRW-PSC-Detection-master', 'template', 'ipsc-template.mat');
                    end

                    params.init_method.threshold = 1.8; %in std of filtered trace (uses the zscore, so any spike above 2.35 std gets marked as an
                                                           %event. This really depends on the noise ofthe data
                    params.init_method.min_interval = 100; %detects events only once in a # frame span.
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
                    params.magnitude = magnitude;

        %% Run Wiener Filter For Trace
                    disp(['Now running inference on ' params.traces_file(1:end-4) ' data...']);

                    %load template for init method
                    load_struct = load(params.init_method.template_file);
                    template = load_struct.template;

                    traces = name.data(stimlocation:stimlocation+15000); %preserve raw trace
                    trace = params.event_sign*traces;
                    trace = trace - min(trace);
                    nfft = length(trace) + length(template);
                    
                    %Run wiener filter to get event times and amplitudes
                    
                    [filtered_trace, event_times, event_amp] = wiener_filter_stim(trace*params.event_sign, params, template, nfft);
                    
                    if isempty(event_times) == 1
                        disp(strcat(params.traces_file(1:end-4), ' was automatically rejected.'));
                        clear name params
                        continue;
                    elseif event_times > 100
                        disp(strcat(params.traces_file(1:end-4), ' was automatically rejected.'));
                        clear name params
                        continue;
                    end
                    
                    plot(name.data(stimlocation-4:stimlocation+1000));
                    title(['Cell ' celll ' Epoch ' epochh ' After Stimulus']);
                    xlim([0 1000]);
                    xlabel('Frames after stimulation onset');
                    hold on
                    
                    name.latency = event_times;
                    name.event_times1 = event_times + stimlocation;
                    name.event_amp1 = event_amp;
                    
                    %These are not needed to find, as this is not PPR
                    name.event_times2 = NaN;
                    name.event_amp2 = NaN;
                    name.PPR = NaN;
                    
                    %Find AUC
                    name.AUC = sum(trace(30:230));
                    
                    %Find Local Extrema after stim
                    if params.event_sign == -1
                        [name.localextrema, name.localextremaloc] = min(name.data(stimlocation:stimlocation+100));
                    else
                        [name.localextrema, name.localextremaloc] = max(name.data(stimlocation:stimlocation+100));
                    end
                    
                    name.extremamagnitude = (name.localextrema-(mean(name.data(stimlocation+230:stimlocation+15230))))*params.event_sign;
                    
                    %Find centers of mass
                    COM = centerOfMass(name.data(stimlocation+30:stimlocation+230));
                    COMFull = centerOfMass(name.data(stimlocation+30:stimlocation+1030));
                    name.center_of_mass = COM(2);
                    name.center_of_mass_full = COMFull(2);
                    
                    baseline = mean(name.data(stimlocation+230:stimlocation+15230));
                    [peak, loc] = min(name.data(name.event_times1:name.event_times1+30));
                    fittingTracee = [name.data(name.event_times1+loc:name.event_times1+loc+200)-baseline];
                    xx = linspace(0, length(fittingTracee)/10, length(fittingTracee));
                    ff = fit(xx', fittingTracee', 'exp2', 'Robust', 'LAR', 'Upper', [-10 0 -5 0]);
                    name.tau1 = max(-ff.b, -ff.d);
                    name.tau2 = NaN;
                    
                    %store initialization results in structure and save
                    params.epoch = num2str(PSCTableDate{i,18});
                    params.cell = num2str(PSCTableDate{i, 3});
                    params.date = expdate;
                    params.location = location;
                    name.experiment = experiment;
                    params.mouseID = PSCTableDate{i, 2};
                    name.name = params.traces_file(1:end-4);
                    name.params = params;
                    name.filtered_trace = filtered_trace;
                    
                    %Assign in to variable and save 
                    assignin('caller', params.traces_file(1:end-4), name);
                    assignin('base', params.traces_file(1:end-4), name);
                    raw_concatenated_traces = [raw_concatenated_traces name];
                    save(params.full_save_string, params.traces_file(1:end-4));
                    
                    %% Write all traces
                    [rows cols] = size(AllTracesTable);
                    roww = rows+1;
                    AllTracesTable{roww, 1} = name.params.date;
                    AllTracesTable{roww, 2} = name.params.mouseID;
                    AllTracesTable{roww, 3} = name.params.location;
                    AllTracesTable{roww, 4} = name.params.cell;
                    AllTracesTable{roww, 5} = name.params.epoch;
                    AllTracesTable{roww, 6} = name.name;
                    AllTracesTable{roww, 7} = name.experiment;
                    AllTracesTable{roww, 8} = name.params.magnitude;
                    AllTracesTable{roww, 9} = name.localextrema;
                    AllTracesTable{roww, 10} = name.localextremaloc;
                    AllTracesTable{roww, 11} = name.extremamagnitude;
                    AllTracesTable{roww, 12} = name.latency;
                    AllTracesTable{roww, 13} = name.AUC;
                    AllTracesTable{roww, 14} = name.center_of_mass;
                    AllTracesTable{roww, 15} = name.center_of_mass_full;
                    AllTracesTable{roww, 16} = name.event_times1;
                    AllTracesTable{roww, 17} = name.event_amp1;
                    AllTracesTable{roww, 18} = name.tau1;
                    AllTracesTable{roww, 19} = name.event_times2;
                    AllTracesTable{roww, 20} = name.event_amp2;
                    AllTracesTable{roww, 21} = name.tau2;
                    AllTracesTable{roww, 22} = name.PPR;
                    
                    if file ~= nACQ
                        disp(' ');
                    else
                        disp(['Done with epoch #' num2str(PSCTableDate{i,18})]);
                    end

                    %clear all variables from this trace to save memory
                    clear(params.traces_file(1:end-4))
                    clear name base_params params results traces trace times event_times event_amp filtered_trace raw_trace iQCSpikeBegin iQCSpikeEnd QC
                    clear -regexp ^AD0_
                end 
                AllTracesTable = array2table(AllTracesTable);
                print(fullfile(savePath1, strcat(mouseID, '_cell', celll, '_epoch', epochh, '_Overlay')), '-dpdf', '-fillpage', '-r1000');
                close all 
                writetable(AllTracesTable, fullfile(outputpath, 'EStim_AllTraces2.xlsx'), 'WriteVariableNames', false);
                
                %% Write Average Traces
                
                if isempty(raw_concatenated_traces) == 1
                    disp(strcat('All of epoch ', epochh, ' was rejected.'));
                    continue;
                end
                
                %initialize variables
                localext = [];
                localextloc = [];
                localextmag = [];
                latencies = [];
                AUCs = [];
                times1 = [];
                amps1 = [];
                times2 = [];
                amps2 = [];
                taus1 = [];
                taus2 = [];
                PPRs =[];
                COMs = [];
                COMs2 = [];
                
                %concatenate
                for j = 1:length(raw_concatenated_traces)
                    localext = [localext raw_concatenated_traces(j).localextrema];
                    localextloc = [localextloc raw_concatenated_traces(j).localextremaloc];
                    localextmag = [localextmag raw_concatenated_traces(j).extremamagnitude];
                    latencies = [latencies raw_concatenated_traces(j).latency];
                    AUCs = [AUCs raw_concatenated_traces(j).AUC];
                    times1 = [times1 raw_concatenated_traces(j).event_times1];
                    amps1 = [amps1 raw_concatenated_traces(j).event_amp1];
                    times2 = [times2 raw_concatenated_traces(j).event_times2];
                    amps2 = [amps2 raw_concatenated_traces(j).event_amp2];
                    taus1 = [taus1 raw_concatenated_traces(j).tau1];
                    taus2 = [taus2 raw_concatenated_traces(j).tau2];
                    PPRs =[PPRs raw_concatenated_traces(j).PPR];
                    COMs = [COMs raw_concatenated_traces(j).center_of_mass];
                    COMs2 = [COMs2 raw_concatenated_traces(j).center_of_mass_full];
                end
                
                %find averages
                avglocalext = mean(localext);
                avglocalextloc = mean(localextloc);
                avglocalextmag = mean(localextmag);
                avglatencies = mean(latencies);
                avgAUCs = mean(AUCs);
                avgtimes1 = mean(times1);
                avgamps1 = mean(amps1);
                avgtimes2 = mean(times2);
                avgamps2 = mean(amps2);
                avgtaus1 = mean(taus1);
                avgtaus2 = mean(taus2);
                avgPPRs =mean(PPRs);
                avgCOMs = mean(COMs);
                avgCOMs2 = mean(COMs2);
                
                %write to excel
                [rows2 cols2] = size(AverageTraces);
                roww2 = rows2+1;
                AverageTraces{roww2+1, 1} = raw_concatenated_traces(1).params.date;
                AverageTraces{roww2+1, 2} = raw_concatenated_traces(1).params.mouseID;
                AverageTraces{roww2+1, 3} = raw_concatenated_traces(1).params.location;
                AverageTraces{roww2+1, 4} = raw_concatenated_traces(1).params.cell;
                AverageTraces{roww2+1, 5} = raw_concatenated_traces(1).params.epoch;
                AverageTraces{roww2+1, 6} = raw_concatenated_traces(1).experiment;
                AverageTraces{roww2+1, 7} = raw_concatenated_traces(1).params.magnitude;
                AverageTraces{roww2+1, 8} = avglocalext;
                AverageTraces{roww2+1, 9} = avglocalextloc;
                AverageTraces{roww2+1, 10} = avglocalextmag;
                AverageTraces{roww2+1, 11} = avglatencies;
                AverageTraces{roww2+1, 12} = avgAUCs;
                AverageTraces{roww2+1, 13} = avgCOMs;
                AverageTraces{roww2+1, 14} = avgCOMs2;
                AverageTraces{roww2+1, 15} = avgtimes1;
                AverageTraces{roww2+1, 16} = avgamps1;
                AverageTraces{roww2+1, 17} = avgtaus1;
                AverageTraces{roww2+1, 18} = avgtimes2;
                AverageTraces{roww2+1, 19} = avgamps2;
                AverageTraces{roww2+1, 20} = avgtaus2;
                AverageTraces{roww2+1, 21} = avgPPRs;
                    
                AverageTraces = array2table(AverageTraces);
                writetable(AverageTraces, fullfile(outputpath, 'EStim_AverageTraces2.xlsx'), 'WriteVariableNames', false);
                
                concatenated_traces = raw_concatenated_traces;
                cellepoch = strcat('concatenated_traces_cell_', celll, '_epoch_', epochh,'.mat');
                full_concatenate_name = fullfile(savePath1, cellepoch);
                save(full_concatenate_name,  'concatenated_traces');
                
                %% Run Quality Control Check
                
                disp('Running Quality Control Check');
                QCs = zeros(length(raw_concatenated_traces), length(raw_concatenated_traces(1).QC));
                for QCCheck = 1:length(raw_concatenated_traces)
                    QCs(QCCheck,:) = raw_concatenated_traces(QCCheck).QC;
                end

                for l = 1:size(QCs,1)
                    RCtrace=QCs(l,1001:end-803); %isolate only the beginning and downward peak of the QC check
                    steadyRC = mean(RCtrace((end-40):end));   %find the mean steady state current of 4ms before upward peak
                    dV = 5; % voltage step in mV
                    [peakRC,peakloc] = min(RCtrace(1:50)); %peak of downward spike 
                    baselineRC(l) = mean(QCs(l,1:999));  %find the mean baseline current

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
                clear AllTracesTable AverageTraces
                clear tau tauCap tauParams Cm Rm Raccess extrapolPkRC fittedVal f fittingTrace2 fittingTrace x x2 steadyRC baselineRC RCtrace locOnset
%% CHR2 PPR                 
            elseif isequal(experiment, 'PPR') == 1
                for file = 1:nACQ
                    %reset params
                    params = [];
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
                    name.peristim1 = name.data(stimlocation:PPRsecondstim-10);
                    name.peristim2 = name.data(PPRsecondstim:PPRsecondstim+490);
                    params.event_sign = -1;
                    subplot(2,1,1)
                    plot(name.peristim1);
                    title('After First Stimulus')
                    hold on;
                    subplot(2,1,2)
                    plot(name.peristim2);
                    title('After Second Stimulus');
                    xlim([0 490]);
                    hold on;

                    params.event_sign = -1;
                    params.init_method.template_file = fullfile('//Volumes', 'Carly Rose', '2 - Code', '1 - MATLAB', 'CRW-PSC-Detection-master', 'template', 'epsc-template.mat');
                    params.init_method.threshold = 1.8; %in std of filtered trace (uses the zscore, so any spike above 2.35 std gets marked as an
                                                       %event. This really depends on the noise ofthe data
                    params.init_method.min_interval = 100; %detects events only once in a 100 frame span.
                    params.dt = 1/10000; %time in seconds per sample
                    params.animal_information = PSCTableDate{i, 4};
                    params.init_method.ar_noise_params.sigma_sq = 3;
                    params.init_method.ar_noise_params.phi = [1.000000000000000, 1.0, -0.15];
                    params.magnitude = magnitude;

                    %update params to reflect current trace being analyzed 
                    params.traces_file = names{file};
                    params.savename = [params.traces_file(1:end-4) '-proc.mat'];
                    params.full_save_string = fullfile(params.savepath, params.savename);
                    params.number_of_trace = file;
                    params.nACQ = nACQ;
                    disp(['Now running inference on ' params.traces_file(1:end-4) ' data...']);

                    %load template for init method
                    load_struct = load(params.init_method.template_file);
                    template = load_struct.template;
                    baseline = mean(name.data(PPRsecondstim+500:PPRsecondstim+10500));
                    
                    trace1 = name.peristim1 * params.event_sign;
                    trace1 = trace1 - min(trace1);
                    nfft1 = length(trace1) + length(template);
                    
                    %run wiener filter on first peristimulus window
                    [filtered_trace1, event_times1] = wiener_filter_PPR(trace1*params.event_sign, params, template, nfft1);
                    event_times1(event_times1>300)=[];
                    name.event_times1 = event_times1;
                    name.filtered_trace1 = filtered_trace1;
                    name.tau1 = [];
                    
                    %these variables are not found because this is PPR
                    name.localextrema = NaN;
                    name.localextremaloc = NaN;
                    name.extremamagnitude = NaN;
                    name.latency = NaN;
                    name.center_of_mass = NaN;
                    name.center_of_mass_full = NaN;
                    name.latency = NaN;
                    name.AUC = NaN;
                    
                    trace2 = name.peristim2 * params.event_sign;
                    trace2 = trace2 - min(trace2);
                    nfft2 = length(trace2) + length(template);
                    
                    %run wiener filter on second peristimulus window
                    [filtered_trace2, event_times2] = wiener_filter_PPR(trace2*params.event_sign, params, template, nfft2);
                    event_times2(event_times2>300)=[];
                    name.event_times2 = event_times2;
                    name.filtered_trace2 = filtered_trace2;
                    name.tau2 = [];

                    if ~isempty(name.event_times2) == 1;
                        name.event_times2 = name.event_times2(1);
                        name.event_amp2 = (min(name.peristim2(1:name.event_times2+75))-baseline)*params.event_sign;
                        [peak, loc] = min(name.peristim2(1:name.event_times2+75));
                        if loc < name.event_times2 
                            name.event_times2 = loc;
                        end 
                        fittingTracee = [name.peristim2(loc:loc+250)-baseline];
                        xx = linspace(0, length(fittingTracee)/10, length(fittingTracee));
                        ff = fit(xx', fittingTracee', 'exp2', 'Robust', 'LAR', 'Upper', [-10 0 -5 0]);
                        name.tau2 = [name.tau2; ff.b ff.d];
                    else
                        name.tau2 = [name.tau2];
                        name.event_amp2 = [];
                    end 
                    name.tau2 = max(name.tau2*-1);
                    
                    clear xx ff peak loc fittingTracee 

                    if ~isempty(name.event_times1) == 1;
                        name.event_times1 = name.event_times1(1);
                        name.event_amp1 = (min(name.peristim1(1:name.event_times1+75))-baseline)*params.event_sign;
                        [peak, loc] = min(name.peristim1(1:name.event_times1+75));
                        if loc < name.event_times1
                            name.event_times1 = loc;
                        end
                        fittingTracee = [name.peristim1(loc:loc+250)-baseline];
                        xx = linspace(0, length(fittingTracee)/10, length(fittingTracee));
                        ff = fit(xx', fittingTracee', 'exp2', 'Robust', 'LAR', 'Upper', [-10 0 -5 0]);
                        name.tau1 = [name.tau1; ff.b ff.d];
                    else
                        name.tau1 = [name.tau1];
                        name.event_amp1 = [];
                    end
                    
                    name.tau1 = max(name.tau1*-1); 

                    clear xx ff peak loc fittingTracee 

                    if ~isempty(name.event_times1) == 1 && ~isempty(name.event_times2)==1
                        name.PPR = name.event_amp2/name.event_amp1;
                    else
                        name.PPR = [];
                    end
                    
                    if name.PPR == [];
                        continue;
                    end

                    %store initialization results in structure and save
                    params.epoch = num2str(PSCTableDate{i,18});
                    params.cell = num2str(PSCTableDate{i, 3});
                    params.date = expdate;
                    params.mouseID = PSCTableDate{i, 2};
                    name.name = params.traces_file(1:end-4);
                    params.location = location;
                    name.experiment = experiment;
                    name.params = params;

                    %assign in variables and save 
                    assignin('caller', params.traces_file(1:end-4), name);
                    assignin('base', params.traces_file(1:end-4), name);
                    raw_concatenated_traces = [raw_concatenated_traces name];
                    save(params.full_save_string, params.traces_file(1:end-4));
                    
                    %% Write all traces
                    [rows cols] = size(AllTracesTable);
                    roww = rows+1;
                    AllTracesTable{roww, 1} = name.params.date;
                    AllTracesTable{roww, 2} = name.params.mouseID;
                    AllTracesTable{roww, 3} = name.params.location;
                    AllTracesTable{roww, 4} = name.params.cell;
                    AllTracesTable{roww, 5} = name.params.epoch;
                    AllTracesTable{roww, 6} = name.name;
                    AllTracesTable{roww, 7} = name.experiment;
                    AllTracesTable{roww, 8} = name.params.magnitude;
                    AllTracesTable{roww, 9} = name.localextrema;
                    AllTracesTable{roww, 10} = name.localextremaloc;
                    AllTracesTable{roww, 11} = name.extremamagnitude;
                    AllTracesTable{roww, 12} = name.latency;
                    AllTracesTable{roww, 13} = name.AUC;
                    AllTracesTable{roww, 14} = name.center_of_mass;
                    AllTracesTable{roww, 15} = name.center_of_mass_full;
                    AllTracesTable{roww, 16} = name.event_times1;
                    AllTracesTable{roww, 17} = name.event_amp1;
                    AllTracesTable{roww, 18} = name.tau1;
                    AllTracesTable{roww, 19} = name.event_times2;
                    AllTracesTable{roww, 20} = name.event_amp2;
                    AllTracesTable{roww, 21} = name.tau2;
                    AllTracesTable{roww, 22} = name.PPR;

                    if file ~= nACQ
                        disp(' ');
                    else
                        disp(['Done with epoch #' num2str(PSCTableDate{i,18})]);
                    end

                    %clear all variables from this trace to save memory
                    clear(params.traces_file(1:end-4))
                    epoch = num2str(params.epoch);
                    clear name params results trace times event_times1 event_amp1 event_times2 event_amp2 filtered_trace1 filtered_trace2 QC 
                    clear -regexp ^AD0_
                    %save overlaid traces
                    print(fullfile(savePath1, strcat(mouseID, '_cell', celll, '_epoch', epochh, '_Overlay_Traces')), '-dpdf', '-fillpage', '-r1000');
                   
                end
                close all
                
                %% Write Average Traces
                if isempty(raw_concatenated_traces) == 1
                    disp(strcat('All of epoch ', epochh, ' was rejected.'));
                    continue;
                end
                
                localext = [];
                localextloc = [];
                localextmag = [];
                latencies = [];
                AUCs = [];
                times1 = [];
                amps1 = [];
                times2 = [];
                amps2 = [];
                taus1 = [];
                taus2 = [];
                PPRs =[];
                COMs = [];
                COMs2 = [];
                
                for j = 1:length(raw_concatenated_traces)
                    localext = [localext raw_concatenated_traces(j).localextrema];
                    localextloc = [localextloc raw_concatenated_traces(j).localextremaloc];
                    localextmag = [localextmag raw_concatenated_traces(j).extremamagnitude];
                    latencies = [latencies raw_concatenated_traces(j).latency];
                    AUCs = [AUCs raw_concatenated_traces(j).AUC];
                    times1 = [times1 raw_concatenated_traces(j).event_times1];
                    amps1 = [amps1 raw_concatenated_traces(j).event_amp1];
                    times2 = [times2 raw_concatenated_traces(j).event_times2];
                    amps2 = [amps2 raw_concatenated_traces(j).event_amp2];
                    taus1 = [taus1 raw_concatenated_traces(j).tau1];
                    taus2 = [taus2 raw_concatenated_traces(j).tau2];
                    PPRs =[PPRs raw_concatenated_traces(j).PPR];
                    COMs = [COMs raw_concatenated_traces(j).center_of_mass];
                    COMs2 = [COMs2 raw_concatenated_traces(j).center_of_mass_full];
                end

                avglocalext = mean(localext);
                avglocalextloc = mean(localextloc);
                avglocalextmag = mean(localextmag);
                avglatencies = mean(latencies);
                avgAUCs = mean(AUCs);
                avgtimes1 = mean(times1);
                avgamps1 = mean(amps1);
                avgtimes2 = mean(times2);
                avgamps2 = mean(amps2);
                avgtaus1 = mean(taus1);
                avgtaus2 = mean(taus2);
                avgPPRs =mean(PPRs);
                avgCOMs = mean(COMs);
                avgCOMs2 = mean(COMs2);

                [rows2 cols2] = size(AverageTraces);
                roww2 = rows2+1;
                AverageTraces{roww2+1, 1} = raw_concatenated_traces(1).params.date;
                AverageTraces{roww2+1, 2} = raw_concatenated_traces(1).params.mouseID;
                AverageTraces{roww2+1, 3} = raw_concatenated_traces(1).params.location;
                AverageTraces{roww2+1, 4} = raw_concatenated_traces(1).params.cell;
                AverageTraces{roww2+1, 5} = raw_concatenated_traces(1).params.epoch;
                AverageTraces{roww2+1, 6} = raw_concatenated_traces(1).experiment;
                AverageTraces{roww2+1, 7} = raw_concatenated_traces(1).params.magnitude;
                AverageTraces{roww2+1, 8} = avglocalext;
                AverageTraces{roww2+1, 9} = avglocalextloc;
                AverageTraces{roww2+1, 10} = avglocalextmag;
                AverageTraces{roww2+1, 11} = avglatencies;
                AverageTraces{roww2+1, 12} = avgAUCs;
                AverageTraces{roww2+1, 13} = avgCOMs;
                AverageTraces{roww2+1, 14} = avgCOMs2;
                AverageTraces{roww2+1, 15} = avgtimes1;
                AverageTraces{roww2+1, 16} = avgamps1;
                AverageTraces{roww2+1, 17} = avgtaus1;
                AverageTraces{roww2+1, 18} = avgtimes2;
                AverageTraces{roww2+1, 19} = avgamps2;
                AverageTraces{roww2+1, 20} = avgtaus2;
                AverageTraces{roww2+1, 21} = avgPPRs;
                
                AverageTraces = array2table(AverageTraces);
                writetable(AverageTraces, fullfile(outputpath, 'EStim_AverageTraces2.xlsx'), 'WriteVariableNames', false);

                AllTracesTable = array2table(AllTracesTable);
                writetable(AllTracesTable, fullfile(outputpath, 'EStim_AllTraces2.xlsx'), 'WriteVariableNames', false);
                
                concatenated_traces = raw_concatenated_traces;
                cellepoch = strcat('concatenated_traces_cell_', celll, '_epoch_', epochh,'.mat');
                full_concatenate_name = fullfile(savePath1, cellepoch);
                save(full_concatenate_name,  'concatenated_traces');
                
                %% Run Quality Control Check

                disp('Running Quality Control Check');
                QCs = zeros(length(raw_concatenated_traces), length(raw_concatenated_traces(1).QC));
                for QCCheck = 1:length(raw_concatenated_traces)
                    QCs(QCCheck,:) = raw_concatenated_traces(QCCheck).QC;
                end

                for l = 1:size(QCs,1)
                    RCtrace=QCs(l,1001:end-1003); %isolate only the beginning and downward peak of the QC check
                    steadyRC = mean(RCtrace((end-40):end));   %find the mean steady state current of 4ms before upward peak
                    dV = 5; % voltage step in mV
                    [peakRC,peakloc] = min(RCtrace(1:50)); %peak of downward spike 
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
                sgtitle(strcat(mouseID, ' cell ', celll, ', epoch ', epochh));
                set(gcf, 'Color', 'w');

            %Adjust and save
                print(fullfile(savePath1, strcat(mouseID, '_cell', celll, '_epoch', epochh, '_QualityControl')), '-dpdf', '-fillpage', '-r1000');

            %Close everything
                close all
                
            %Clear unnecessary variables
                clear tau tauCap tauParams Cm Rm Raccess extrapolPkRC peakRC peakloc names acsweeps nACQ fittedVal f fittingTrace2 fittingTrace x x2 steadyRC baselineRC RCtrace locOnset
                clear AllTracesTable AverageTraces

%% CHR2 I/0                
            elseif isequal(experiment, 'I/O') == 1
                figure
                for file = 1:nACQ
                    %reset params
                    params = [];
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
                    
                    if holding_current == -70
                        params.event_sign = -1;
                    else
                        params.event_sign = 1;
                    end

                    %Use event sign gathered from excel in line 52
                    params.init_method.template_file = fullfile('//Volumes', 'Carly Rose', '2 - Code', '1 - MATLAB', 'CRW-PSC-Detection-master', 'template', 'epsc-template.mat');

                    params.init_method.threshold = 1.8; %in std of filtered trace (uses the zscore, so any spike above 2.35 std gets marked as an
                                                           %event. This really depends on the noise ofthe data
                    params.init_method.min_interval = 100; %detects events only once in a 170 frame span.
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
                    params.magnitude = magnitude;

        %% Run Wiener Filter For Trace
                    disp(['Now running inference on ' params.traces_file(1:end-4) ' data...']);

                    %load template for init method
                    load_struct = load(params.init_method.template_file);
                    template = load_struct.template;

                    traces = name.data(stimlocation:stimlocation+15000); %preserve raw trace
                    trace = params.event_sign*traces;
                    trace = trace - min(trace);
                    nfft = length(trace) + length(template);
                    
                    %Run wiener filter to get event times and amplitudes
                    [filtered_trace, event_times, event_amp] = wiener_filter_stim(trace*params.event_sign, params, template, nfft);
                    
                    if isempty(event_times) == 1
                        disp(strcat(params.traces_file(1:end-4), ' was automatically rejected.'));
                        clear name params
                        continue;
                    elseif event_times > 100
                        disp(strcat(params.traces_file(1:end-4), ' was automatically rejected.'));
                        clear name params
                        continue;
                    end
                    
                    plot(name.data(stimlocation-4:stimlocation+1000));
                    title(['Cell ' celll ' Epoch ' epochh ' After Stimulus']);
                    xlim([0 1000]);
                    xlabel('Frames after stimulation onset');
                    hold on
                    
                    name.latency = event_times;
                    name.event_times1 = event_times + stimlocation;
                    name.event_amp1 = event_amp;
                    
                    %These are not needed to find, as this is not PPR
                    name.event_times2 = NaN;
                    name.event_amp2 = NaN;
                    name.PPR = NaN;
                    
                    %Find AUC
                    
                    name.AUC = sum(trace(30:230));
                    
                    %Find Local Extrema after stim
                    if params.event_sign == -1
                        [name.localextrema, name.localextremaloc] = min(name.data(stimlocation:stimlocation+100));
                    else
                        [name.localextrema, name.localextremaloc] = max(name.data(stimlocation:stimlocation+100));
                    end
                    
                    name.extremamagnitude = (name.localextrema-(mean(name.data(stimlocation+230:stimlocation+15230))))*params.event_sign;
                    
                    %Find centers of mass
                    COM = centerOfMass(name.data(stimlocation+30:stimlocation+230));
                    COMFull = centerOfMass(name.data(stimlocation+30:stimlocation+1030));
                    name.center_of_mass = COM(2);
                    name.center_of_mass_full = COMFull(2);
                    
                    baseline = mean(name.data(stimlocation+230:stimlocation+15230));
                    [peak, loc] = min(name.data(name.event_times1:name.event_times1+30));
                    fittingTracee = [name.data(name.event_times1+loc:name.event_times1+loc+200)-baseline];
                    xx = linspace(0, length(fittingTracee)/10, length(fittingTracee));
                    ff = fit(xx', fittingTracee', 'exp2', 'Robust', 'LAR', 'Upper', [-10 0 -5 0]);
                    name.tau1 = max(-ff.b, -ff.d);
                    name.tau2 = NaN;
                    
                    %store initialization results in structure and save
                    params.epoch = num2str(PSCTableDate{i,18});
                    params.cell = num2str(PSCTableDate{i, 3});
                    params.date = expdate;
                    params.mouseID = PSCTableDate{i, 2};
                    name.name = params.traces_file(1:end-4);
                    params.location = location;
                    name.experiment = experiment;
                    name.params = params;

                    name.filtered_trace = filtered_trace;
                    
                    %assign in variables and save
                    assignin('caller', params.traces_file(1:end-4), name);
                    assignin('base', params.traces_file(1:end-4), name);
                    raw_concatenated_traces = [raw_concatenated_traces name];
                    save(params.full_save_string, params.traces_file(1:end-4));
                    
                    %% Write All Traces
                    [rows cols] = size(AllTracesTable);
                    roww = rows+1;
                    AllTracesTable{roww, 1} = name.params.date;
                    AllTracesTable{roww, 2} = name.params.mouseID;
                    AllTracesTable{roww, 3} = name.params.location;
                    AllTracesTable{roww, 4} = name.params.cell;
                    AllTracesTable{roww, 5} = name.params.epoch;
                    AllTracesTable{roww, 6} = name.name;
                    AllTracesTable{roww, 7} = name.experiment;
                    AllTracesTable{roww, 8} = name.params.magnitude;
                    AllTracesTable{roww, 9} = name.localextrema;
                    AllTracesTable{roww, 10} = name.localextremaloc;
                    AllTracesTable{roww, 11} = name.extremamagnitude;
                    AllTracesTable{roww, 12} = name.latency;
                    AllTracesTable{roww, 13} = name.AUC;
                    AllTracesTable{roww, 14} = name.center_of_mass;
                    AllTracesTable{roww, 15} = name.center_of_mass_full;
                    AllTracesTable{roww, 16} = name.event_times1;
                    AllTracesTable{roww, 17} = name.event_amp1;
                    AllTracesTable{roww, 18} = name.tau1;
                    AllTracesTable{roww, 19} = name.event_times2;
                    AllTracesTable{roww, 20} = name.event_amp2;
                    AllTracesTable{roww, 21} = name.tau2;
                    AllTracesTable{roww, 22} = name.PPR;

                    if file ~= nACQ
                        if file == nACQ-1
                            disp(' ');
                        else
                            disp(' ');
                        end
                    else
                        disp(['Done with epoch #' num2str(PSCTableDate{i,18})]);
                    end

                    %clear all variables from this trace to save memory
                    clear(params.traces_file(1:end-4))
                    clear name base_params params results traces trace times event_times event_amp filtered_trace raw_trace iQCSpikeBegin iQCSpikeEnd QC
                    clear -regexp ^AD0_
                end 
                
                print(fullfile(savePath1, strcat(mouseID, '_cell', celll, '_epoch', epochh, '_Overlay')), '-dpdf', '-fillpage', '-r1000');
                close all
                
                %% Write Average Traces
                
                if isempty(raw_concatenated_traces) == 1
                    disp(strcat('All of epoch ', epochh, ' was rejected.'));
                    continue;
                end
                
                localext = [];
                localextloc = [];
                localextmag = [];
                latencies = [];
                AUCs = [];
                times1 = [];
                amps1 = [];
                times2 = [];
                amps2 = [];
                taus1 = [];
                taus2 = [];
                PPRs =[];
                COMs = [];
                COMs2 = [];
                
                for j = 1:length(raw_concatenated_traces)
                    localext = [localext raw_concatenated_traces(j).localextrema];
                    localextloc = [localextloc raw_concatenated_traces(j).localextremaloc];
                    localextmag = [localextmag raw_concatenated_traces(j).extremamagnitude];
                    latencies = [latencies raw_concatenated_traces(j).latency];
                    AUCs = [AUCs raw_concatenated_traces(j).AUC];
                    times1 = [times1 raw_concatenated_traces(j).event_times1];
                    amps1 = [amps1 raw_concatenated_traces(j).event_amp1];
                    times2 = [times2 raw_concatenated_traces(j).event_times2];
                    amps2 = [amps2 raw_concatenated_traces(j).event_amp2];
                    taus1 = [taus1 raw_concatenated_traces(j).tau1];
                    taus2 = [taus2 raw_concatenated_traces(j).tau2];
                    PPRs =[PPRs raw_concatenated_traces(j).PPR];
                    COMs = [COMs raw_concatenated_traces(j).center_of_mass];
                    COMs2 = [COMs2 raw_concatenated_traces(j).center_of_mass_full];
                end

                avglocalext = mean(localext);
                avglocalextloc = mean(localextloc);
                avglocalextmag = mean(localextmag);
                avglatencies = mean(latencies);
                avgAUCs = mean(AUCs);
                avgtimes1 = mean(times1);
                avgamps1 = mean(amps1);
                avgtimes2 = mean(times2);
                avgamps2 = mean(amps2);
                avgtaus1 = mean(taus1);
                avgtaus2 = mean(taus2);
                avgPPRs =mean(PPRs);
                avgCOMs = mean(COMs);
                avgCOMs2 = mean(COMs2);

                [rows2 cols2] = size(AverageTraces);
                roww2 = rows2+1;
                AverageTraces{roww2+1, 1} = raw_concatenated_traces(1).params.date;
                AverageTraces{roww2+1, 2} = raw_concatenated_traces(1).params.mouseID;
                AverageTraces{roww2+1, 3} = raw_concatenated_traces(1).params.location;
                AverageTraces{roww2+1, 4} = raw_concatenated_traces(1).params.cell;
                AverageTraces{roww2+1, 5} = raw_concatenated_traces(1).params.epoch;
                AverageTraces{roww2+1, 6} = raw_concatenated_traces(1).experiment;
                AverageTraces{roww2+1, 7} = raw_concatenated_traces(1).params.magnitude;
                AverageTraces{roww2+1, 8} = avglocalext;
                AverageTraces{roww2+1, 9} = avglocalextloc;
                AverageTraces{roww2+1, 10} = avglocalextmag;
                AverageTraces{roww2+1, 11} = avglatencies;
                AverageTraces{roww2+1, 12} = avgAUCs;
                AverageTraces{roww2+1, 13} = avgCOMs;
                AverageTraces{roww2+1, 14} = avgCOMs2;
                AverageTraces{roww2+1, 15} = avgtimes1;
                AverageTraces{roww2+1, 16} = avgamps1;
                AverageTraces{roww2+1, 17} = avgtaus1;
                AverageTraces{roww2+1, 18} = avgtimes2;
                AverageTraces{roww2+1, 19} = avgamps2;
                AverageTraces{roww2+1, 20} = avgtaus2;
                AverageTraces{roww2+1, 21} = avgPPRs;
                AverageTraces = array2table(AverageTraces);     
                writetable(AverageTraces, fullfile(outputpath, 'EStim_AverageTraces2.xlsx'), 'WriteVariableNames', false);

                AllTracesTable = array2table(AllTracesTable);
                writetable(AllTracesTable, fullfile(outputpath, 'EStim_AllTraces2.xlsx'), 'WriteVariableNames', false);
                
                concatenated_traces = raw_concatenated_traces;
                cellepoch = strcat('concatenated_traces_cell_', celll, '_epoch_', epochh,'.mat');
                full_concatenate_name = fullfile(savePath1, cellepoch);
                save(full_concatenate_name, 'concatenated_traces');
                
                %% Run Quality Control   
                
                disp('Running Quality Control Check');
                QCs = zeros(length(raw_concatenated_traces), length(raw_concatenated_traces(1).QC));
                for QCCheck = 1:length(raw_concatenated_traces)
                    QCs(QCCheck,:) = raw_concatenated_traces(QCCheck).QC;
                end

                for l = 1:size(QCs,1)
                    RCtrace=QCs(l,1001:end-803); %isolate only the beginning and downward peak of the QC check
                    steadyRC = mean(RCtrace((end-40):end));   %find the mean steady state current of 4ms before upward peak
                    dV = 5; % voltage step in mV
                    [peakRC,peakloc] = min(RCtrace(1:50)); %peak of downward spike 
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
                clear tau tauCap tauParams Cm Rm Raccess extrapolPkRC peakRC peakloc names acsweeps nACQ fittedVal f fittingTrace2 fittingTrace x x2 steadyRC baselineRC RCtrace locOnset
                clear AllTracesTable AverageTraces
%% CHR2 E-E/I                
            elseif isequal(experiment, 'E-E/I') == 1
                figure
                for file = 1:nACQ
                    %reset params
                    params = [];
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
                    
                    params.event_sign = -1;
                    params.init_method.template_file = fullfile('//Volumes', 'Carly Rose', '2 - Code', '1 - MATLAB', 'CRW-PSC-Detection-master', 'template', 'epsc-template.mat');

                    params.init_method.threshold = 1.8; %in std of filtered trace (uses the zscore, so any spike above 2.35 std gets marked as an
                                                           %event. This really depends on the noise ofthe data
                    params.init_method.min_interval = 100; %detects events only once in a 170 frame span.
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
                    params.magnitude = magnitude;

        %% Run Wiener Filter For Trace
                    disp(['Now running inference on ' params.traces_file(1:end-4) ' data...']);

                    %load template for init method
                    load_struct = load(params.init_method.template_file);
                    template = load_struct.template;

                    traces = name.data(stimlocation:stimlocation+15000); %preserve raw trace
                    trace = params.event_sign*traces;
                    trace = trace - min(trace);
                    nfft = length(trace) + length(template);
                    
                    %Run wiener filter to get event times and amplitudes
                    [filtered_trace, event_times, event_amp] = wiener_filter_stim(trace*params.event_sign, params, template, nfft);
                    
                    if isempty(event_times) == 1
                        disp(strcat(params.traces_file(1:end-4), ' was automatically rejected.'));
                        clear name params
                        continue;
                    elseif event_times > 100
                        disp(strcat(params.traces_file(1:end-4), ' was automatically rejected.'));
                        clear name params
                        continue;
                    end
                    
                    plot(name.data(stimlocation-4:stimlocation+1000));
                    title(['Cell ' celll ' Epoch ' epochh ' After Stimulus']);
                    xlim([0 1000]);
                    xlabel('Frames after stimulation onset');
                    hold on
                    
                    name.latency = event_times;
                    name.event_times1 = event_times + stimlocation;
                    name.event_amp1 = event_amp;
                    
                    %These are not needed to find, as this is not PPR
                    name.event_times2 = NaN;
                    name.event_amp2 = NaN;
                    name.PPR = NaN;
                    
                    %Find AUC
                    name.AUC = sum(trace(30:230));
                    
                    %Find Local Extrema after stim
                    if params.event_sign == -1
                        [name.localextrema, name.localextremaloc] = min(name.data(stimlocation:stimlocation+100));
                    else
                        [name.localextrema, name.localextremaloc] = max(name.data(stimlocation:stimlocation+100));
                    end
                    
                    name.extremamagnitude = (name.localextrema-(mean(name.data(stimlocation+230:stimlocation+15230))))*params.event_sign;
                    
                    %Find centers of mass
                    COM = centerOfMass(name.data(stimlocation+30:stimlocation+230));
                    COMFull = centerOfMass(name.data(stimlocation+30:stimlocation+1030));
                    name.center_of_mass = COM(2);
                    name.center_of_mass_full = COMFull(2);
                    
                    baseline = mean(name.data(stimlocation+230:stimlocation+15230));
                    [peak, loc] = min(name.data(name.event_times1:name.event_times1+30));
                    fittingTracee = [name.data(name.event_times1+loc:name.event_times1+loc+200)-baseline];
                    xx = linspace(0, length(fittingTracee)/10, length(fittingTracee));
                    ff = fit(xx', fittingTracee', 'exp2', 'Robust', 'LAR', 'Upper', [-10 0 -5 0]);
                    name.tau1 = max(-ff.b, -ff.d);
                    name.tau2 = NaN;
                    
                    %store initialization results in structure and save
                    params.epoch = num2str(PSCTableDate{i,18});
                    params.cell = num2str(PSCTableDate{i, 3});
                    params.date = expdate;
                    params.mouseID = PSCTableDate{i, 2};
                    name.name = params.traces_file(1:end-4);
                    params.location = location;
                    name.experiment = experiment;
                    name.params = params;
                    name.filtered_trace = filtered_trace;
                    
                    %Assign in variables and save
                    assignin('caller', params.traces_file(1:end-4), name);
                    assignin('base', params.traces_file(1:end-4), name);
                    raw_concatenated_traces = [raw_concatenated_traces name];
                    save(params.full_save_string, params.traces_file(1:end-4));
                    
                    %% Write All Traces
                    [rows cols] = size(AllTracesTable);
                    roww = rows+1;
                    AllTracesTable{roww, 1} = name.params.date;
                    AllTracesTable{roww, 2} = name.params.mouseID;
                    AllTracesTable{roww, 3} = name.params.location;
                    AllTracesTable{roww, 4} = name.params.cell;
                    AllTracesTable{roww, 5} = name.params.epoch;
                    AllTracesTable{roww, 6} = name.name;
                    AllTracesTable{roww, 7} = name.experiment;
                    AllTracesTable{roww, 8} = name.params.magnitude;
                    AllTracesTable{roww, 9} = name.localextrema;
                    AllTracesTable{roww, 10} = name.localextremaloc;
                    AllTracesTable{roww, 11} = name.extremamagnitude;
                    AllTracesTable{roww, 12} = name.latency;
                    AllTracesTable{roww, 13} = name.AUC;
                    AllTracesTable{roww, 14} = name.center_of_mass;
                    AllTracesTable{roww, 15} = name.center_of_mass_full;
                    AllTracesTable{roww, 16} = name.event_times1;
                    AllTracesTable{roww, 17} = name.event_amp1;
                    AllTracesTable{roww, 18} = name.tau1;
                    AllTracesTable{roww, 19} = name.event_times2;
                    AllTracesTable{roww, 20} = name.event_amp2;
                    AllTracesTable{roww, 21} = name.tau2;
                    AllTracesTable{roww, 22} = name.PPR;

                    if file ~= nACQ
                        disp(' ');
                    else
                        disp(['Done with epoch #' num2str(PSCTableDate{i,18})]);
                    end

                    %clear all variables from this trace to save memory
                    clear(params.traces_file(1:end-4))
                    clear name base_params params results traces trace times event_times event_amp filtered_trace raw_trace iQCSpikeBegin iQCSpikeEnd QC
                    clear -regexp ^AD0_
                end 
                print(fullfile(savePath1, strcat(mouseID, '_cell', celll, '_epoch', epochh, '_Overlay')), '-dpdf', '-fillpage', '-r1000');
                close all
                
                %% Write Average Traces
                
                if isempty(raw_concatenated_traces) == 1
                    disp(strcat('All of epoch ', epochh, ' was rejected.'));
                    continue;
                end
                
                localext = [];
                localextloc = [];
                localextmag = [];
                latencies = [];
                AUCs = [];
                times1 = [];
                amps1 = [];
                times2 = [];
                amps2 = [];
                taus1 = [];
                taus2 = [];
                PPRs =[];
                COMs = [];
                COMs2 = [];
                
                for j = 1:length(raw_concatenated_traces)
                    localext = [localext raw_concatenated_traces(j).localextrema];
                    localextloc = [localextloc raw_concatenated_traces(j).localextremaloc];
                    localextmag = [localextmag raw_concatenated_traces(j).extremamagnitude];
                    latencies = [latencies raw_concatenated_traces(j).latency];
                    AUCs = [AUCs raw_concatenated_traces(j).AUC];
                    times1 = [times1 raw_concatenated_traces(j).event_times1];
                    amps1 = [amps1 raw_concatenated_traces(j).event_amp1];
                    times2 = [times2 raw_concatenated_traces(j).event_times2];
                    amps2 = [amps2 raw_concatenated_traces(j).event_amp2];
                    taus1 = [taus1 raw_concatenated_traces(j).tau1];
                    taus2 = [taus2 raw_concatenated_traces(j).tau2];
                    PPRs =[PPRs raw_concatenated_traces(j).PPR];
                    COMs = [COMs raw_concatenated_traces(j).center_of_mass];
                    COMs2 = [COMs2 raw_concatenated_traces(j).center_of_mass_full];
                end

                avglocalext = mean(localext);
                avglocalextloc = mean(localextloc);
                avglocalextmag = mean(localextmag);
                avglatencies = mean(latencies);
                avgAUCs = mean(AUCs);
                avgtimes1 = mean(times1);
                avgamps1 = mean(amps1);
                avgtimes2 = mean(times2);
                avgamps2 = mean(amps2);
                avgtaus1 = mean(taus1);
                avgtaus2 = mean(taus2);
                avgPPRs =mean(PPRs);
                avgCOMs = mean(COMs);
                avgCOMs2 = mean(COMs2);

                [rows2 cols2] = size(AverageTraces);
                roww2 = rows2+1;
                AverageTraces{roww2+1, 1} = raw_concatenated_traces(1).params.date;
                AverageTraces{roww2+1, 2} = raw_concatenated_traces(1).params.mouseID;
                AverageTraces{roww2+1, 3} = raw_concatenated_traces(1).params.location;
                AverageTraces{roww2+1, 4} = raw_concatenated_traces(1).params.cell;
                AverageTraces{roww2+1, 5} = raw_concatenated_traces(1).params.epoch;
                AverageTraces{roww2+1, 6} = raw_concatenated_traces(1).experiment;
                AverageTraces{roww2+1, 7} = raw_concatenated_traces(1).params.magnitude;
                AverageTraces{roww2+1, 8} = avglocalext;
                AverageTraces{roww2+1, 9} = avglocalextloc;
                AverageTraces{roww2+1, 10} = avglocalextmag;
                AverageTraces{roww2+1, 11} = avglatencies;
                AverageTraces{roww2+1, 12} = avgAUCs;
                AverageTraces{roww2+1, 13} = avgCOMs;
                AverageTraces{roww2+1, 14} = avgCOMs2;
                AverageTraces{roww2+1, 15} = avgtimes1;
                AverageTraces{roww2+1, 16} = avgamps1;
                AverageTraces{roww2+1, 17} = avgtaus1;
                AverageTraces{roww2+1, 18} = avgtimes2;
                AverageTraces{roww2+1, 19} = avgamps2;
                AverageTraces{roww2+1, 20} = avgtaus2;
                AverageTraces{roww2+1, 21} = avgPPRs;

                AverageTraces = array2table(AverageTraces);
                writetable(AverageTraces, fullfile(outputpath, 'EStim_AverageTraces2.xlsx'), 'WriteVariableNames', false);

                AllTracesTable = array2table(AllTracesTable);
                writetable(AllTracesTable, fullfile(outputpath, 'EStim_AllTraces2.xlsx'), 'WriteVariableNames', false);
                
                concatenated_traces = raw_concatenated_traces;
                cellepoch = strcat('concatenated_traces_cell_', celll, '_epoch_', epochh,'.mat');
                full_concatenate_name = fullfile(savePath1, cellepoch);
                save(full_concatenate_name,  'concatenated_traces');
                
                %% Run Quality Control
                
                disp('Running Quality Control Check');
                QCs = zeros(length(raw_concatenated_traces), length(raw_concatenated_traces(1).QC));
                for QCCheck = 1:length(raw_concatenated_traces)
                    QCs(QCCheck,:) = raw_concatenated_traces(QCCheck).QC;
                end

                for l = 1:size(QCs,1)
                    RCtrace=QCs(l,1001:end-803); %isolate only the beginning and downward peak of the QC check
                    steadyRC = mean(RCtrace((end-40):end));   %find the mean steady state current of 4ms before upward peak
                    dV = 5; % voltage step in mV
                    [peakRC,peakloc] = min(RCtrace(1:50)); %peak of downward spike 
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
                clear tau tauCap tauParams Cm Rm Raccess extrapolPkRC fittedVal f fittingTrace2 fittingTrace x x2 steadyRC baselineRC RCtrace locOnset
                clear AverageTraces AllTracesTable
%% CHR2 I-E/I                
            elseif isequal(experiment, 'I-E/I') == 1
                figure
                for file = 1:nACQ
                    %reset params
                    params = [];
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
                    
                    params.event_sign = 1;
                    params.init_method.template_file = fullfile('//Volumes', 'Carly Rose', '2 - Code', '1 - MATLAB', 'CRW-PSC-Detection-master', 'template', 'ipsc-template.mat');

                    params.init_method.threshold = 1.8; %in std of filtered trace (uses the zscore, so any spike above 2.35 std gets marked as an
                                                           %event. This really depends on the noise ofthe data
                    params.init_method.min_interval = 100; %detects events only once in a 170 frame span.
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
                    params.magnitude = magnitude;

        %% Run Wiener Filter For Trace
                    disp(['Now running inference on ' params.traces_file(1:end-4) ' data...']);

                    %load template for init method
                    load_struct = load(params.init_method.template_file);
                    template = load_struct.template;

                    traces = name.data(stimlocation:stimlocation+15000); %preserve raw trace
                    trace = params.event_sign*traces;
                    trace = trace - min(trace);
                    nfft = length(trace) + length(template);
                    
                    %Run wiener filter to get event times and amplitudes
                    [filtered_trace, event_times, event_amp] = wiener_filter_stim(trace*params.event_sign, params, template, nfft);
                    
                    if isempty(event_times) == 1
                        disp(strcat(params.traces_file(1:end-4), ' was automatically rejected.'));
                        clear name params
                        continue;
                    elseif event_times > 100
                        disp(strcat(params.traces_file(1:end-4), ' was automatically rejected.'));
                        clear name params
                        continue;
                    end
                    
                    plot(name.data(stimlocation-4:stimlocation+1000));
                    title(['Cell ' celll ' Epoch ' epochh ' After Stimulus']);
                    xlim([0 1000]);
                    xlabel('Frames after stimulation onset');
                    hold on
                    
                    name.latency = event_times;
                    name.event_times1 = event_times + stimlocation;
                    name.event_amp1 = event_amp;
                    
                    %These are not needed to find, as this is not PPR
                    name.event_times2 = NaN;
                    name.event_amp2 = NaN;
                    name.PPR = NaN;
                    
                    %Find AUC
                    
                    name.AUC = sum(trace(30:230));
                    
                    %Find Local Extrema after stim
                    if params.event_sign == -1
                        [name.localextrema, name.localextremaloc] = min(name.data(stimlocation:stimlocation+100));
                    else
                        [name.localextrema, name.localextremaloc] = max(name.data(stimlocation:stimlocation+100));
                    end
                    
                    name.extremamagnitude = (name.localextrema-(mean(name.data(stimlocation+230:stimlocation+15230))))*params.event_sign;
                    
                    %Find centers of mass
                    COM = centerOfMass(name.data(stimlocation+30:stimlocation+230));
                    COMFull = centerOfMass(name.data(stimlocation+30:stimlocation+1030));
                    name.center_of_mass = COM(2);
                    name.center_of_mass_full = COMFull(2);
                    
                    baseline = mean(name.data(stimlocation+230:stimlocation+15230));
                    [peak, loc] = min(name.data(name.event_times1:name.event_times1+30));
                    fittingTracee = [name.data(name.event_times1+loc:name.event_times1+loc+200)-baseline];
                    xx = linspace(0, length(fittingTracee)/10, length(fittingTracee));
                    ff = fit(xx', fittingTracee', 'exp2', 'Robust', 'LAR', 'Upper', [-10 0 -5 0]);
                    name.tau1 = max(-ff.b, -ff.d);
                    name.tau2 = NaN;
                    
                    %store initialization results in structure and save
                    params.epoch = num2str(PSCTableDate{i,18});
                    params.cell = num2str(PSCTableDate{i, 3});
                    params.date =expdate;
                    params.mouseID = PSCTableDate{i, 2};
                    name.name = params.traces_file(1:end-4);
                    params.location = location;
                    name.experiment = experiment;
                    name.params = params;
                    name.filtered_trace = filtered_trace;
                    
                    %Assign in and save variables
                    assignin('caller', params.traces_file(1:end-4), name);
                    assignin('base', params.traces_file(1:end-4), name);
                    raw_concatenated_traces = [raw_concatenated_traces name];
                    save(params.full_save_string, params.traces_file(1:end-4));
                    
                    %% Write all traces
                    [rows cols] = size(AllTracesTable);
                    roww = rows+1;
                    AllTracesTable{roww, 1} = name.params.date;
                    AllTracesTable{roww, 2} = name.params.mouseID;
                    AllTracesTable{roww, 3} = name.params.location;
                    AllTracesTable{roww, 4} = name.params.cell;
                    AllTracesTable{roww, 5} = name.params.epoch;
                    AllTracesTable{roww, 6} = name.name;
                    AllTracesTable{roww, 7} = name.experiment;
                    AllTracesTable{roww, 8} = name.params.magnitude;
                    AllTracesTable{roww, 9} = name.localextrema;
                    AllTracesTable{roww, 10} = name.localextremaloc;
                    AllTracesTable{roww, 11} = name.extremamagnitude;
                    AllTracesTable{roww, 12} = name.latency;
                    AllTracesTable{roww, 13} = name.AUC;
                    AllTracesTable{roww, 14} = name.center_of_mass;
                    AllTracesTable{roww, 15} = name.center_of_mass_full;
                    AllTracesTable{roww, 16} = name.event_times1;
                    AllTracesTable{roww, 17} = name.event_amp1;
                    AllTracesTable{roww, 18} = name.tau1;
                    AllTracesTable{roww, 19} = name.event_times2;
                    AllTracesTable{roww, 20} = name.event_amp2;
                    AllTracesTable{roww, 21} = name.tau2;
                    AllTracesTable{roww, 22} = name.PPR;

                    if file ~= nACQ
                        disp(' ');
                    else
                        disp(['Done with epoch #' num2str(PSCTableDate{i,18})]);
                    end

                    %clear all variables from this trace to save memory
                    clear(params.traces_file(1:end-4))
                    
                    clear name base_params params results traces trace times event_times event_amp filtered_trace raw_trace iQCSpikeBegin iQCSpikeEnd QC
                    clear -regexp ^AD0_
                end 
                print(fullfile(savePath1, strcat(mouseID, '_cell', celll, '_epoch', epochh, '_Overlay')), '-dpdf', '-fillpage', '-r1000');
                close all
                
                %% Write Average Traces
                
                if isempty(raw_concatenated_traces) == 1
                    disp(strcat('All of epoch ', epochh, ' was rejected.'));
                    continue;
                end
                
                localext = [];
                localextloc = [];
                localextmag = [];
                latencies = [];
                AUCs = [];
                times1 = [];
                amps1 = [];
                times2 = [];
                amps2 = [];
                taus1 = [];
                taus2 = [];
                PPRs =[];
                COMs = [];
                COMs2 = [];
                
                for j = 1:length(raw_concatenated_traces)
                    localext = [localext raw_concatenated_traces(j).localextrema];
                    localextloc = [localextloc raw_concatenated_traces(j).localextremaloc];
                    localextmag = [localextmag raw_concatenated_traces(j).extremamagnitude];
                    latencies = [latencies raw_concatenated_traces(j).latency];
                    AUCs = [AUCs raw_concatenated_traces(j).AUC];
                    times1 = [times1 raw_concatenated_traces(j).event_times1];
                    amps1 = [amps1 raw_concatenated_traces(j).event_amp1];
                    times2 = [times2 raw_concatenated_traces(j).event_times2];
                    amps2 = [amps2 raw_concatenated_traces(j).event_amp2];
                    taus1 = [taus1 raw_concatenated_traces(j).tau1];
                    taus2 = [taus2 raw_concatenated_traces(j).tau2];
                    PPRs =[PPRs raw_concatenated_traces(j).PPR];
                    COMs = [COMs raw_concatenated_traces(j).center_of_mass];
                    COMs2 = [COMs2 raw_concatenated_traces(j).center_of_mass_full];
                end

                avglocalext = mean(localext);
                avglocalextloc = mean(localextloc);
                avglocalextmag = mean(localextmag);
                avglatencies = mean(latencies);
                avgAUCs = mean(AUCs);
                avgtimes1 = mean(times1);
                avgamps1 = mean(amps1);
                avgtimes2 = mean(times2);
                avgamps2 = mean(amps2);
                avgtaus1 = mean(taus1);
                avgtaus2 = mean(taus2);
                avgPPRs =mean(PPRs);
                avgCOMs = mean(COMs);
                avgCOMs2 = mean(COMs2);

                [rows2 cols2] = size(AverageTraces);
                roww2 = rows2+1;
                AverageTraces{roww2+1, 1} = raw_concatenated_traces(1).params.date;
                AverageTraces{roww2+1, 2} = raw_concatenated_traces(1).params.mouseID;
                AverageTraces{roww2+1, 3} = raw_concatenated_traces(1).params.location;
                AverageTraces{roww2+1, 4} = raw_concatenated_traces(1).params.cell;
                AverageTraces{roww2+1, 5} = raw_concatenated_traces(1).params.epoch;
                AverageTraces{roww2+1, 6} = raw_concatenated_traces(1).experiment;
                AverageTraces{roww2+1, 7} = raw_concatenated_traces(1).params.magnitude;
                AverageTraces{roww2+1, 8} = avglocalext;
                AverageTraces{roww2+1, 9} = avglocalextloc;
                AverageTraces{roww2+1, 10} = avglocalextmag;
                AverageTraces{roww2+1, 11} = avglatencies;
                AverageTraces{roww2+1, 12} = avgAUCs;
                AverageTraces{roww2+1, 13} = avgCOMs;
                AverageTraces{roww2+1, 14} = avgCOMs2;
                AverageTraces{roww2+1, 15} = avgtimes1;
                AverageTraces{roww2+1, 16} = avgamps1;
                AverageTraces{roww2+1, 17} = avgtaus1;
                AverageTraces{roww2+1, 18} = avgtimes2;
                AverageTraces{roww2+1, 19} = avgamps2;
                AverageTraces{roww2+1, 20} = avgtaus2;
                AverageTraces{roww2+1, 21} = avgPPRs;
                
                AverageTraces = array2table(AverageTraces);
                writetable(AverageTraces, fullfile(outputpath, 'EStim_AverageTraces2.xlsx'), 'WriteVariableNames', false);

                AllTracesTable = array2table(AllTracesTable);
                writetable(AllTracesTable, fullfile(outputpath, 'EStim_AllTraces2.xlsx'), 'WriteVariableNames', false);
                
                concatenated_traces = raw_concatenated_traces;
                cellepoch = strcat('concatenated_traces_cell_', celll, '_epoch_', epochh,'.mat');
                full_concatenate_name = fullfile(savePath1, cellepoch);
                save(full_concatenate_name,  'concatenated_traces');

                %% Run Quality Control 
                
                disp('Running Quality Control Check');
                QCs = zeros(length(raw_concatenated_traces), length(raw_concatenated_traces(1).QC));
                
                for QCCheck = 1:length(raw_concatenated_traces)
                    QCs(QCCheck,:) = raw_concatenated_traces(QCCheck).QC;
                end

                for l = 1:size(QCs,1)
                    RCtrace=QCs(l,1001:end-803); %isolate only the beginning and downward peak of the QC check
                    steadyRC = mean(RCtrace((end-40):end));   %find the mean steady state current of 4ms before upward peak
                    dV = 5; % voltage step in mV
                    [peakRC,peakloc] = min(RCtrace(1:50)); %peak of downward spike 
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
                clear tau tauCap tauParams Cm Rm Raccess extrapolPkRC fittedVal f fittingTrace2 fittingTrace x x2 steadyRC baselineRC RCtrace locOnset
                clear AllTracesTable AverageTraces
             
            end
            
            
warning('on', 'MATLAB:unknownObjectNowStruct');
warning('on', 'MATLAB:MKDIR:DirectoryExists');
warning('on', 'MATLAB:centerOfMass:neg');
    end
disp('Processing Complete');
end
