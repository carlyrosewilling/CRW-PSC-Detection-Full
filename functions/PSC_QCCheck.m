function resistance = PSC_QCCheck(QC)

    start = find(QC == min(QC));
    startbase = start + 100;
    ends = find(QC == max(QC));
    endbase = ends - 10;
    baseline = mean(QC(startbase:endbase));
    idelta = baseline - max(QC);
    resistance = -5/idelta;
end
