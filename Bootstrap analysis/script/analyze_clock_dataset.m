function result = analyze_clock_dataset(filePath, sheetName, genotype, cfg)
%ANALYZE_CLOCK_DATASET Run model fit and both bootstrap analyses.

arguments
    filePath (1,1) string
    sheetName (1,1) string
    genotype (1,1) string
    cfg struct
end

fprintf("\n=== %s | %s cells | %s ===\n", genotype, upper(sheetName), cfg.modelName);
inputSheet = cfg.sheetMap.(lower(genotype)).(lower(sheetName));
data = read_clock_long_sheet(filePath, inputSheet, cfg);
fprintf("  n=%d, time points=%d, family=%s, denominator=%g\n", ...
    data.n, data.nTimePoints, data.family, data.denominator);

observedFit = fit_clock_model(data, cfg);
fprintf("  Observed deviance improvement = %.6g\n", observedFit.observedStatistic);

nullBootstrap = bootstrap_clock_rhythmicity(data, cfg, observedFit);
fprintf("  Rhythmicity bootstrap p = %.6g\n", nullBootstrap.pValue);

amplitudeBootstrap = bootstrap_clock_amplitude(data, cfg, observedFit);
fprintf("  Peak-to-trough amplitude = %.6g; one-sided U95 = %.6g\n", ...
    observedFit.amplitude, amplitudeBootstrap.amplitudeUpper95);

result = struct();
result.genotype = genotype;
result.signal = lower(sheetName);
result.data = data;
result.fit = observedFit;
result.nullBootstrap = nullBootstrap;
result.amplitudeBootstrap = amplitudeBootstrap;
end
