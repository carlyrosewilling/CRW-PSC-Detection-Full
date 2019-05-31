%% Pipeline for Bessel Filtering Data%%

%Written by CRW, 15 Oct 2018
    %last updated: 16 Oct 2018

warning('off', 'MATLAB:unknownObjectNowStruct');
warning('off', 'MATLAB:MKDIR:DirectoryExists');

%% Initialize %%

%User selects exact file where data is stored
    %prepath = fullfile('/Users', 'carlyrosewilling', 'Desktop', 'Data', 'Minis', 'KM092618/');
    [~,prepath] = uiputfile('*.*','Select data folder', 'datapath.mat');

    if isnumeric(prepath)
        disp('User must select an input');
        return
    end
	
%User selects where to save output files
    [~, savePath]=uiputfile('output.mat', 'Select path for output');
    %savePath = fullfile('/Users', 'carlyrosewilling', 'Desktop', 'Analysis', 'Minis', 'KM092618_output/');

%User enters date of aquisition for parsing Excel Sheet
    date = input('Input date of collection (i.e. 20180906) ');
    %date = 20180926;
    
%Loads overview Excel sheet 
    %User selects Excel file
    PSC_LoadExcel
    
    PSCTableDate = {};
    
    for j = 2:length(PSCTableRaw)
        if PSCTableRaw{j,1} == date
            row = horzcat(PSCTableRaw(j,:));
            PSCTableDate = vertcat(PSCTableDate, row);
        else
            PSCTableDate = PSCTableDate;
        end
    end
    [nrows ncolumns] = size(PSCTableDate);
    
%User enters which Epochs to run
    Epochs = input(['You have ', num2str(nrows) ' epochs. Input epochs to run in matrix form ']);
    %Epochs = [3];
    disp(' ');
    
    [b a] = besself(4, (3000/(1/(2*pi))), 'low'); %4-pole analog bessel lowpass with 3kHz cutoff (cutoff is in rad/s)
    [numd dend] = impinvar(b, a, 10000); %converts analog filter to digital filter, sampled at 10 kHz
    
%% Organize by Epoch 
    %input for acquisitions and creation of filename directories
    for i = Epochs;
        disp(['Filtering epoch #' num2str(PSCTableDate{i,17})]);
        disp('------------------------------');
        acqsweeps = PSCTableDate{i, 15}:PSCTableDate{i,16};
        nACQ = length(acqsweeps);
        names = cell(1, nACQ); %preallocate size of names
        %lists file names for every sweep in Epoch
        for nameindex = 1:nACQ
            names{nameindex} = strcat('AD0_', num2str(acqsweeps(nameindex)), '.mat');
        end

%% Set Parameters for Each Epoch and load files

        for file = 1:nACQ
            params.traces_filename = strcat(prepath, names{file});
            
            if exist(params.traces_filename) == 0
                continue;
            else 
                load(params.traces_filename);
            end 
              
            name = eval(names{file}(1:end-4));
            
            trace = name.data;
            params.traces_file = names{file};
            
            iQCSpikeBegin = find(name.data == min(name.data));
            iQCSpikeEnd = find(name.data == max(name.data));
            QC = name.data(iQCSpikeBegin-2:iQCSpikeEnd+250);
            name.rawQC = QC;
            
            name.data = filtfilt(numd,dend,trace); %applies filter
            savepath = fullfile(savePath, params.traces_file);
            
            assignin('base', params.traces_file(1:end-4), name);
            
            save(savepath, params.traces_file(1:end-4)); %saves filtered trace in output file with same name
            
            
            clear(params.traces_file(1:end-4))
            clear name params trace data
        end
        clear names acsweeps nACQ
    end
    disp('~~Done with all filtering!~~');
warning('on', 'MATLAB:unknownObjectNowStruct');
warning('on', 'MATLAB:MKDIR:DirectoryExists');