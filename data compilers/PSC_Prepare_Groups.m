%Prepares a structure with all traces for individual groups

%list group titles (also the file structure and naming structure)
    ages = {'P15', 'P20', 'P30', 'P45', 'P65', 'P81', 'P120'};
    
%containing folder of concatenated name trace
    prepath = fullfile('//Volumes', 'Neurobio', 'MICROSCOPE', 'Kevin', '3-Experiments', '4-SliceEphys', '6-WT Normal Development', '4-Good Traces');

%loads trace files from each group, and adds them to one structure
    for i = 1:length(ages)
        load(fullfile(prepath, ages{i}, strcat('full_', ages{i}, '.mat')));
        name2 = eval(strcat('full_', ages{i}));
        alltraces.(ages{i}) = name2;
    end
    
%save structure to output folder
    save(fullfile(prepath, 'alltraces.mat'), 'alltraces');
    
    