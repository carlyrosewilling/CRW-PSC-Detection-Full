%% Concatenate and Calculate PPR Info %%
%Creates mat file containing only accepted traces and their results. Plots
    %average traces for AMPA and NMDA currents, and writes important data
    %to excel sheet.
    
%Written by CRW, 8 March 2019
    %Last Updated 8 March 2019

    
%User inputs date of acquisition, and cell and epoch to process
    date = input('Input date of recording (i.e. 01/06/2019): ', 's');
    recorder = input('KM or WW?', 's');
    cell = input('Input cell: ', 's');
    epoch = input('Input epoch: ', 's');
    
%Makes input path given date information
    datedfolder = strcat(recorder, date(1:2), date(4:5), date(9:end), '_output');
    cellfolder = strcat('cell_', cell);
    epochfolder = strcat('epoch_', epoch);
    prepath = fullfile('//Volumes', 'Neurobio', 'MICROSCOPE', 'Kevin', '3-Experiments', '4-SliceEphys', '7-Paired Pulse', '2-Output', datedfolder, cellfolder, epochfolder);
    filename = strcat('Accepted_Traces_cell', cell, '_epoch', epoch, '.xlsx');

    [~, ~, aPSCTableRaw] = xlsread(fullfile(prepath, filename));
    
    %Concatenates traces
    concatenated_traces = [];
    for index = 2:length(aPSCTableRaw)-1 
        name = aPSCTableRaw{index,1};
        load(fullfile(prepath, name));
        name2 = eval(name(1:end-9));
        concatenated_traces = [concatenated_traces name2];
        clear(name(1:end-9));
    end


%Saves concatenated trace
save(fullfile(prepath, strcat('Proc_Concatenated_Traces_cell', cell, '_epoch',epoch, '.mat')), 'concatenated_traces');

%Initialize to find averages
peristims1 = [];
peristims2 = [];
eventtimes1 = [];
eventtimes2 = [];
eventamps1 = [];
eventamps2 = [];
PPR = [];
for j = 1:length(concatenated_traces)
    if concatenated_traces(1).experiment == -20
        peristims2 = [peristims2 ; concatenated_traces(j).peristim2];
        eventtimes2 = [eventtimes2 concatenated_traces(j).event_times2];
        eventamps2 = [eventamps2 concatenated_traces(j).event_amp2];
        PPR = [PPR concatenated_traces(j).PPR];
        
    end
    peristims1 = [peristims1 ; concatenated_traces(j).peristim1]; 
    eventtimes1 = [eventtimes1 concatenated_traces(j).event_times1];
    eventamps1 = [eventamps1 concatenated_traces(j).event_amp1];
end

%Average peristim 1 
avg_peristim1 = mean(peristims1);


%Average event time 1
avg_event_time1 = mean(eventtimes1);

%Average event amp 1
avg_event_amp1 = mean(eventamps1);


%[filename2 pathname2 ~] = uigetfile({'*.xlsx', '*xls*'}, 'Select Excel file with cell summaries');
filename2 = 'Cell Summary-Paired Pulse.xlsx';
pathname2 = '/Volumes/Neurobio/MICROSCOPE/Kevin/3-Experiments/4-SliceEphys/7-Paired Pulse/2-Output/';
[~, ~, cellTableRaw] = xlsread(fullfile(pathname2, filename2));
[height width] = size(cellTableRaw);

cellTableNew = cellTableRaw;
cellTableNew{1+height, 1} = concatenated_traces(1).params.date;
cellTableNew{1+height, 2} = concatenated_traces(1).params.mouseID;
cellTableNew{1+height, 3} = concatenated_traces(1).params.animal_information;
cellTableNew{1+height, 4} = concatenated_traces(1).experiment;
cellTableNew{1+height, 5} = cell;
cellTableNew{1+height, 6} = epoch;
cellTableNew{1+height, 7} = length(concatenated_traces);
cellTableNew{1+height, 8} = concatenated_traces(1).params.init_method.threshold;
cellTableNew{1+height, 9} = avg_event_time1;
cellTableNew{1+height, 10} = avg_event_amp1;
figure

if concatenated_traces(1).experiment == -20
    %Average event time 2
    avg_event_time2 = mean(eventtimes2);
    
    %Average event amp 2
    avg_event_amp2 = mean(eventamps2);

    %Average PPR 
    avg_PPR = mean(PPR);
    
    cellTableNew{1+height, 11} =  avg_event_time2;
    cellTableNew{1+height, 12} =  avg_event_amp2;
    cellTableNew{1+height, 13} = avg_PPR;
    
    %Average peristim 2
    avg_peristim2 = mean(peristims2);
    subplot(2,1,1)
    plot(avg_peristim1)
    xlim([0 500])
    title('Average Trace After First Stimulus')
    subplot(2,1,2)
    plot(avg_peristim2)
    xlim([0 500])
    title('Average Trace After Second Stimulus')
else
    cellTableNew{1+height, 11} =  'NA';
    cellTableNew{1+height, 12} =  'NA';
    cellTableNew{1+height, 13} = 'NA';
    plot(avg_peristim1)
    xlim([0 700])
    title('Average Trace');
end

cellTableFinal = array2table(cellTableNew);
writetable(cellTableFinal, fullfile(pathname2, filename2), 'WriteVariableNames', false);
print(fullfile(prepath, strcat(concatenated_traces(1).params.mouseID, '_cell', cell, '_epoch', epoch, '_Average_Traces')), '-dpdf', '-fillpage', '-r1000');
clear all;
close all;

