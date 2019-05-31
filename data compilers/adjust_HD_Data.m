
concatenated_traces = rmfield(concatenated_traces, 'fullbesseled');
concatenated_traces = rmfield(concatenated_traces, 'besseled');
concatenated_traces = rmfield(concatenated_traces, 'besseled_QC');
for i = 1:length(concatenated_traces)
    concatenated_traces(i).QC = concatenated_traces(i).raw_QC;
end
concatenated_traces = rmfield(concatenated_traces, 'raw_QC');
concatenated_traces = orderfields(concatenated_traces, full_group4);
full_group4 = [full_group4 concatenated_traces];