MLDiff = 1/mean(P120ISIs)-1/mean(P15ISIs);
ISIs = [P120ISIs' P15ISIs'];
Nall = length(ISIs);
N15 = length(P15ISIs);
N120 = length(P120ISIs);

for i = 1:1000
    samp15 = ISIs(randsample(Nall, N15, 1));
    samp120 = ISIs(randsample(Nall, N120, 1));
    sampdiff(i) = 1/mean(samp120)-1/mean(samp15);
end
figure
hist(sampdiff, 30);
line([MLDiff MLDiff], [0 100])