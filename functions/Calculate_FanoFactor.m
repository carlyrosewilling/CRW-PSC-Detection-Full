all_event_times_P120 = [];
for i = 1:length(alltraces.P120)
    all_event_times_P120 = [all_event_times_P120 alltraces.P120(i).event_times];
end
    
all_event_times_P15 = [];
for j = 1:length(alltraces.P15)
    all_event_times_P15 = [all_event_times_P15 alltraces.P15(j).event_times];
end

P120_SpikeTrain = zeros(1, 100000);
for k = 1:length(all_event_times_P120)
    P120_SpikeTrain(all_event_times_P120(k)) = P120_SpikeTrain(all_event_times_P120(k))+1;
end

P15_SpikeTrain = zeros(1,100000);
for l = 1:length(all_event_times_P15)
    P15_SpikeTrain(all_event_times_P15(l)) = P15_SpikeTrain(all_event_times_P15(l))+1;
end

all_event_times_P20 = [];
for m = 1:length(alltraces.P20)
    all_event_times_P20 = [all_event_times_P20 alltraces.P20(m).event_times];
end

all_event_times_P30 = [];
for n = 1:length(alltraces.P30)
    all_event_times_P30 = [all_event_times_P30 alltraces.P30(n).event_times];
end

all_event_times_P45 = [];
for o = 1:length(alltraces.P45)
    all_event_times_P45 = [all_event_times_P45 alltraces.P45(o).event_times];
end

all_event_times_P65 = [];
for p = 1:length(alltraces.P65)
    all_event_times_P65 = [all_event_times_P65 alltraces.P65(p).event_times];
end

all_event_times_P81 = [];
for q = 1:length(alltraces.P81)
    all_event_times_P81 = [all_event_times_P81 alltraces.P81(q).event_times];
end


increments_P15 = hist(all_event_times_P15, bins);
FF_P15 = var(increments_P15)/mean(increments_P15);

increments_P20 = hist(all_event_times_P20, bins);
FF_P20 = var(increments_P20)/mean(increments_P20);

increments_P30 = hist(all_event_times_P30, bins);
FF_P30 = var(increments_P30)/mean(increments_P30);

increments_P45 = hist(all_event_times_P45, bins);
FF_P45 = var(increments_P45)/mean(increments_P45);

increments_P65 = hist(all_event_times_P65, bins);
FF_P65 = var(increments_P65)/mean(increments_P65);

increments_P81 = hist(all_event_times_P81, bins);
FF_P81 = var(increments_P81)/mean(increments_P81);

increments_P120 = hist(all_event_times_P120, bins);
FF_P120 = var(increments_P120)/mean(increments_P120);

    