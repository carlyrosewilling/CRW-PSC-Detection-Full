function name = ReProcess(name, noiseparam, threshold)
	load_struct = load(name.params.init_method.template_file);
	template = load_struct.template;
            
	traces = name.raw_trace; %preserve raw trace
	trace = name.params.event_sign*traces;
	trace = trace - min(trace);

	%run wiener filter to find init event times
	nfft = length(trace) + length(template);
            
	name.params.init_method.threshold = threshold;
    name.params.init_method.ar_noise_params.phi(3) = noiseparam;
        
	[filtered_trace, event_times, event_amp] = wiener_filter(trace*name.params.event_sign, name.params, template, nfft);

	%store initialization results in structure and save
	name.event_times = event_times;
	name.event_amp = event_amp;
	if isempty(name.event_times) == 1
        name.ISIs = [];
    else
        name.ISIs = [event_times(1) diff(event_times)];
    end 
            
	name.SpikeTrain = zeros([1 length(name.raw_trace)]);
            
	for times = event_times;
        name.SpikeTrain(times) = 1;
    end

	name.filtered_trace = filtered_trace;
            
	%assignin('base', name.params.traces_file(1:end-4), name);
            
	%save(name.params.full_save_string, name.params.traces_file(1:end-4), '-append');
	%clear name
end
