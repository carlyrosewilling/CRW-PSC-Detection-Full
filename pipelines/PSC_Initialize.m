%% Shell for Processing Data %%

%Written by CRW, 27 May 2019
    %last updated: 6 March 2019
    
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

%User inputs date of acquisition
    prompt = {'Enter date of recording (i.e. 01/06/2019):', 'Enter Recorder:', 'Experiment:'};
    dlgtitle = 'Inputs';
    dims = [1 75];
    definput = {'01/06/2019', 'WW or KM', 'WT Minis, Isolated Minis, PPR, CHR2'};
    answer = inputdlg(prompt, dlgtitle, dims, definput);

    date = answer{1};
    recorder = answer{2};
    experiment = answer{3};

%User enters date of aquisition for parsing Excel Sheet
    tabledate = str2num(strcat(date(7:end), date(1:2), date(4:5)));
    
    %User enters date of aquisition for parsing Excel Sheet
    tabledate = str2num(strcat(date(7:end), date(1:2), date(4:5)));
    
%Loads overview Excel sheet 
    %Uses user selected experiment type to read appropriate excel file
    if isequal(experiment, 'PPR') == 1
        filename = 'Overview_wavebook_pairedpulse.xlsx';
        pathname = '//Volumes/Neurobio/MICROSCOPE/Kevin/3-Experiments/4-SliceEphys/7-Paired Pulse/1-Raw Data/';
        expfold = '7-Paired Pulse';
    elseif isequal(experiment, 'WT Minis') == 1
        filename = 'Overview_wavebook WT.xlsx';
        pathname = '//Volumes/Neurobio/MICROSCOPE/Kevin/3-Experiments/4-SliceEphys/6-WT Normal Development/1-Raw Data/';
        expfold = '6-WT Normal Development';
    elseif isequal(experiment, 'Isolated Minis') == 1
        filename = 'Overview_wavebook_isolatedmice copy.xlsx';
        pathname =  '//Volumes/Neurobio/MICROSCOPE/Kevin/3-Experiments/4-SliceEphys/5-Isolated Mice/1-Raw Data/';
        expfold = '5-Isolated Mice';
    elseif isequal(experiment, 'CHR2') == 1
        filename = 'Overview_ChR2.xlsx';
        pathname = '//Volumes/Neurobio/MICROSCOPE/Kevin/3-Experiments/4-SliceEphys/8-CHR2/1-Raw Data/';
        expfold = '8-CHR2';
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
    
    if isequal(experiment, 'PPR') == 1
        PSC_Process_PPR(PSCTableDate, datedfolder)
    elseif isequal(experiment, 'WT Minis') == 1
        PSC_WT_Minis(PSCTableDate, datedfolder)
    elseif isequal(experiment, 'Isolated Minis') == 1
        PSC_Isolated_Minis(PSCTableDate, datedfolder)
    elseif isequal(experiment, 'CHR2') == 1
        PSC_Process_CHR2(PSCTableDate, datedfolder)
    end 
    
warning('on', 'MATLAB:unknownObjectNowStruct');
warning('on', 'MATLAB:MKDIR:DirectoryExists');

    
