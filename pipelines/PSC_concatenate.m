%% Concatenate and Organize Accepted Traces %%
%Creates mat file containing only accepted traces and their params,
    %raw_traces, and results. To index, load concatenated_traces.mat and
    %index as concatenated_traces(1).params for the first traces params
    %etc. Change lines 83 and 84 to match local path.
    
%Written by CRW, 5 Oct 2018
    %Last Updated 10 Oct 2018

    
%User inputs date of acquisition, and cell and epoch to process
    date = input('Input date of recording (i.e. 01/06/2019): ', 's');
    recorder = input('KM or WW?', 's');
    cell = input('Input cell: ', 's');
    epoch = input('Input epoch: ', 's');
    
%Makes input path given date information
    datedfolder = strcat(recorder, date(1:2), date(4:5), date(9:end), '_output');
    cellfolder = strcat('cell_', cell);
    epochfolder = strcat('epoch_', epoch);
    prepath = fullfile('//Volumes', 'Neurobio', 'MICROSCOPE', 'Kevin', '3-Experiments', '4-SliceEphys', '9-Plexicon', '2-Output', datedfolder, cellfolder, epochfolder);
    filename = strcat('Accepted_Traces_cell', cell, '_epoch', epoch, '.xlsx');

    [~, ~, aPSCTableRaw] = xlsread(fullfile(prepath, filename));

%Concatenates traces
concatenated_traces = [];
for index = 2:length(aPSCTableRaw)
    name = aPSCTableRaw{index,1};
    load(fullfile(prepath, name));
    name2 = eval(name(1:end-9));
    if isfield(name2, 'QC') == 1
        name2.Resistance = PSC_QCCheck(name2.QC);
    else
        name2.Resistance = PSC_QCCheck(name2.raw_QC);
    end
    concatenated_traces = [concatenated_traces name2];
    clear(name(1:end-9));
end

%event_amp exclusion
for i = 1:length(concatenated_traces);
    concatenated_traces(i).event_times(concatenated_traces(i).event_amp < 7) =[];
    concatenated_traces(i).event_amp(concatenated_traces(i).event_amp < 7) =[];
end

x = concatenated_traces(1).params;

%Calculates number of events, average amplitude, freq of events, and
    %average ISI for each sweep.
for i = 1:length(concatenated_traces);
    concatenated_traces(i).num_events = length(concatenated_traces(i).event_times);
    concatenated_traces(i).avg_amplitude = mean(concatenated_traces(i).event_amp);
    concatenated_traces(i).frequency_of_events = (concatenated_traces(i).num_events/length(concatenated_traces(i).rawdata))*(1/concatenated_traces(i).params.dt);
    concatenated_traces(i).avg_ISI = mean(concatenated_traces(i).ISIs);
    concatenated_traces(i).average_resistance = mean(concatenated_traces(i).Resistance);
end 

%Saves concatenated trace
save(fullfile(prepath, strcat('Concatenated_Traces_cell', x.cell, '_epoch',x.epoch, '.mat')), 'concatenated_traces');


%Initialize to find averages
amps = [];
isi = [];
numev = [];
freqs = [];
resistance = [];
for j = 1:length(concatenated_traces)
    amps = [amps concatenated_traces(j).avg_amplitude];
    isi = [isi concatenated_traces(j).avg_ISI];
    numev = [numev concatenated_traces(j).num_events];
    freqs = [freqs concatenated_traces(j).frequency_of_events];
    resistance = [resistance concatenated_traces(j).average_resistance];
end

%Average amplitude for this cell =
avg_amplitude_for_cell = mean(amps);

%Average ISI for this cell =
avg_ISI_for_cell = mean(isi);

%Average number of events per sweep for this cell = 
avg_num_events_for_cell = mean(numev);

%total number of events for entire cell =
total_num_events_for_cell = sum(numev);

%Average frequency of events for this cell =
avg_freq_for_cell = mean(freqs);

%Average resistance of cell
avg_resistance = mean(resistance);
std_resistance = std(resistance);

%[filename2 pathname2 ~] = uigetfile({'*.xlsx', '*xls*'}, 'Select Excel file with cell summaries');
filename2 = 'Cell Summary-Plexicon.xlsx';
pathname2 = '/Volumes/Neurobio/MICROSCOPE/Kevin/3-Experiments/4-SliceEphys/9-Plexicon/2-Output/';
[~, ~, cellTableRaw] = xlsread(fullfile(pathname2, filename2));
[height width] = size(cellTableRaw);

cellTableNew = cellTableRaw;
cellTableNew{1+height, 1} = concatenated_traces(1).params.date;
cellTableNew{1+height, 2} = concatenated_traces(1).params.mouseID;
cellTableNew{1+height, 3} = concatenated_traces(1).params.location; 
cellTableNew{1+height, 4} = str2num(concatenated_traces(1).params.cell);
cellTableNew{1+height, 5} = str2num(concatenated_traces(1).params.epoch);
cellTableNew{1+height, 6} = length(concatenated_traces);
cellTableNew{1+height, 7} = avg_amplitude_for_cell;
cellTableNew{1+height, 8} = avg_ISI_for_cell;
cellTableNew{1+height, 9} = avg_freq_for_cell;
cellTableNew{1+height, 10} = avg_num_events_for_cell;
cellTableNew{1+height, 11} = total_num_events_for_cell;
cellTableNew{1+height, 12} = concatenated_traces(1).params.init_method.threshold;
cellTableNew{1+height, 13} = avg_resistance;
cellTableNew{1+height, 14} = std_resistance;
cellTableNew{1+height, 15} = x.event_sign;
cellTableFinal = array2table(cellTableNew);
writetable(cellTableFinal, fullfile(pathname2, filename2), 'WriteVariableNames', false);
clear all;


    
    