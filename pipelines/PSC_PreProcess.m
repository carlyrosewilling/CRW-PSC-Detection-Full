%% PreProcess Pipeline %%

%Written by CRW, 14 May 2019
    %last updated: 14 May 2019
    
    %lines 36 through 48 need to be updated to define the paths
        %based on if you're logged on to the server and what not
  	%you can also add new "Experiment" types to this code block, to access
        %the correct excel files and save appropriately.
    %change lines 68 and 71 to define proper paths.

    
%% Turn off dumb warning for loading "wave" struct and directory
warning('off', 'MATLAB:unknownObjectNowStruct');
warning('off', 'MATLAB:MKDIR:DirectoryExists');

%% Initialize %%

%User inputs date of acquisition
    prompt = {'Enter date of recording (i.e. 01/06/2019):', 'Enter experiment type (i.e. PPR, WT minis, Isolated minis, CHR2):', 'Enter Recorder:', 'Enter QC Pulse Beginning:', 'Enter QC Pulse End:'};
    dlgtitle = 'Inputs';
    dims = [1 75];
    definput = {'01/06/2019', 'PPR, WT minis, Isolated minis, CHR2, EStim', 'WW or KM', 'Index of first peak', 'Index of second peak'};
    answer = inputdlg(prompt, dlgtitle, dims, definput);

    date = answer{1};
    Experiment = answer{2};
    recorder = answer{3};
    QCBegin = str2num(answer{4});
    QCEnd = str2num(answer{5});
    
%User enters date of aquisition for parsing Excel Sheet
    tabledate = str2num(strcat(date(7:end), date(1:2), date(4:5)));
    
%Loads overview Excel sheet 
    %Uses user selected experiment type to read appropriate excel file
    if isequal(Experiment, 'PPR') == 1
        filename = 'Overview_wavebook_pairedpulse.xlsx';
        pathname = '//Volumes/Neurobio/MICROSCOPE/Kevin/3-Experiments/4-SliceEphys/7-Paired Pulse/1-Raw Data/';
        expfold = '7-Paired Pulse';
        ep = 13;
        ce = 3;
        acq1 = 11;
        acq2 = 12;
        exp = 5;
    elseif isequal(Experiment, 'WT minis') == 1
        filename = 'Overview_wavebook WT.xlsx';
        pathname = '//Volumes/Neurobio/MICROSCOPE/Kevin/3-Experiments/4-SliceEphys/6-WT Normal Development/1-Raw Data/';
        expfold = '6-WT Normal Development';
        ep = 17;
        ce = 3;
        acq1 = 15;
        acq2 = 16;
        exp = 4;
    elseif isequal(Experiment, 'Isolated minis') == 1
        filename = 'Overview_wavebook_isolatedmice copy.xlsx';
        pathname =  '//Volumes/Neurobio/MICROSCOPE/Kevin/3-Experiments/4-SliceEphys/5-Isolated Mice/1-Raw Data/';
        expfold = '5-Isolated Mice';
        ep = 17;
        ce = 3;
        acq1 = 15;
        acq2 = 16;
        exp = 4;
    elseif isequal(Experiment, 'CHR2') == 1
        filename = 'Overview_ChR2_copy.xlsx';
        pathname = '//Volumes/Neurobio/MICROSCOPE/Kevin/3-Experiments/4-SliceEphys/8-CHR2/1-Raw Data/';
        expfold = '8-CHR2';
        ep = 18;
        ce = 3;
        acq1 = 16;
        acq2 = 17;
        exp = 13; 
    elseif isequal(Experiment, 'EStim') == 1
        filename = 'Overview_wavebook_EStim3.xlsx';
        pathname = '//Volumes/Neurobio/MICROSCOPE/Kevin/3-Experiments/4-SliceEphys/10-EStim/1-Raw Data/';
        expfold = '10-EStim';
        ep = 18;
        ce = 3;
        acq1 = 16;
        acq2 = 17;
        exp = 13; 
    elseif isequal(Experiment, 'Plexicon Minis') == 1
        filename = 'Overview_wavebook_Plexicon.xlsx';
        pathname = '//Volumes/Neurobio/MICROSCOPE/Kevin/3-Experiments/4-SliceEphys/9-Plexicon/1-Raw Data/';
        expfold = '9-Plexicon';
        ep = 17;
        ce = 3;
        acq1 = 15;
        acq2 = 16;
        exp = 4; 
    else
        disp('Invalid Experiment Choice');    
    end
    
    PSC_LoadExcel(pathname, filename);
    
    recorderdate = strcat(recorder, date(1:2), date(4:5), date(9:10));
    
    PSCTableDate = {};
    [row column] = size(PSCTableRaw);
    for j = 2:row
        if PSCTableRaw{j,1} == tabledate 
            if PSCTableRaw{j,2} == recorderdate
                row = horzcat(PSCTableRaw(j,:));
                PSCTableDate = vertcat(PSCTableDate, row);
            else
                PSCTableDate = PSCTableDate;
            end
        else
            PSCTableDate = PSCTableDate;
        end
    end
    [nrows ncolumns] = size(PSCTableDate);
    
%Uses mouse ID to determine folder of data location and output    
    datedfolder = PSCTableDate{1,2};

%Makes input path given date information
    prepath = fullfile('//Volumes', 'Neurobio', 'MICROSCOPE', 'Kevin', '3-Experiments', '4-SliceEphys', expfold, '1-Raw Data', datedfolder);

%Makes save path given date information
    savePath = fullfile('//Volumes', 'Neurobio', 'MICROSCOPE', 'Kevin', '3-Experiments', '4-SliceEphys', expfold, 'Preprocessed Data', datedfolder);
    mkdir(savePath);

%User enters which Epochs to run
    Epochs = input(['You have ', num2str(nrows) ' epochs. Input epochs to run in matrix form ']);
    disp(' ');
%% Organize by Epoch 
    %input for acquisitions and creation of filename directories
    for i = Epochs;
        disp(['Beginning epoch #' num2str(PSCTableDate{i,ep})]); %18 for CHR2
        disp('------------------------------');
        acqsweeps = PSCTableDate{i, acq1}:PSCTableDate{i,acq2}; %16:17 for CHR2
        nACQ = length(acqsweeps);
        experiment = PSCTableDate{i, exp}; %13 for CHR2
        mouseID = PSCTableDate{i, 2};
        epoch = num2str(PSCTableDate{i,ep}); %18
        celll = num2str(PSCTableDate{i, ce});
        names = cell(1, nACQ); %preallocate size of names
        savePath1 = fullfile(savePath, strcat('cell_', num2str(PSCTableDate{i, ce})), strcat('epoch_', num2str(PSCTableDate{i,ep}))); %18 
        mkdir(savePath1);
        
        %lists file names for every sweep in Epoch
        for nameindex = 1:nACQ
            names{nameindex} = strcat('AD0_', num2str(acqsweeps(nameindex)), '.mat');
        end

%% Set Parameters for Each Epoch and load files
        raw_concatenated_traces = [];
        for file = 1:nACQ
            %initialize for get_params function
            traces_filename = strcat(prepath, '/', names{file});
            savepath = savePath1;
            traces_file = names{file};
            
            try 
                load(traces_filename);
            catch
                fprintf(strcat(names{file}, ' failed to load... may not exist... skipping'));
                disp(' ');
                continue;
            end 
            
            name = eval(names{file}(1:end-4));
            
            fields = {'xscale', 'yscale', 'zscale', 'plot', 'UserData', 'note', 'timeStamp', 'holdUpdates', 'needsReplot'};
            name = rmfield(name, fields);
            name.filename = traces_file;
            name.QC = name.data(QCBegin-1000:QCEnd+800);
            
            % Determines where raw data is based on location of QC Pulse
            if QCBegin < (1/3)*length(name.data);
                name.rawdata = name.data(QCEnd+997:end);
            else
                name.rawdata = name.data(1:QCBegin-1003);
            end

        
%% Save QC and data to new cell/epoch folder
            
            assignin('base', traces_file(1:end-4), name);
            
            raw_concatenated_traces = [raw_concatenated_traces name];

            save(fullfile(savePath1, traces_file), traces_file(1:end-4));
            
%% Clear Variables and Move on to Next Trace

            %displays to user based on current stage of analysis
            if file ~= nACQ
                if file == nACQ-1
                    disp('On to the next! Only one more to go!')
                    disp('------------------------------')
                else
                    disp(['**On to the next! ' num2str(nACQ-file) ' traces left for epoch ' num2str(PSCTableDate{i,ep}) '***']);
                    disp('------------------------------');
                end
            else
                disp(['Done with epoch #' num2str(PSCTableDate{i,ep})]); %18
            end
           
            %clear all variables from this trace to save memory
            clear(traces_file(1:end-4))
            clear name 
        
        end
        
%% Quality Control %%         
        
        disp('Running Quality Control Check');
        QCs = zeros(length(raw_concatenated_traces), length(raw_concatenated_traces(1).QC));
        for QCCheck = 1:length(raw_concatenated_traces)
            QCs(QCCheck,:) = raw_concatenated_traces(QCCheck).QC;
        end
        
        for l = 1:size(QCs,1)
            RCtrace=QCs(l,1001:end-803); %isolate only the beginning and downward peak of the QC check
            steadyRC = mean(RCtrace((end-40):end));   %find the mean steady state current of 4ms before upward peak
            dV = 5; % voltage step in mV
            [peakRC,peakloc] = min(RCtrace(1:900)); %peak of downward spike 
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
        sgtitle(strcat(mouseID, ' cell ', celll, ', epoch ', num2str(i)));
        set(gcf, 'Color', 'w');

    %Adjust and save
        print(fullfile(savePath1, strcat(mouseID, '_cell', celll, '_epoch', epoch, '_QualityControl')), '-dpdf', '-fillpage', '-r1000');

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
