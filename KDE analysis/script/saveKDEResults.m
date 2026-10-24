function result = saveKDEResults(result, cfg)
% SAVEKDERESULTS Keep the effective parameters and numerical results together.
result.configuration = cfg;
result.matlabVersion = version;
save(fullfile(cfg.outputDir, cfg.outputStem + "_results.mat"), 'result', '-v7.3');
fid = fopen(fullfile(cfg.outputDir, cfg.outputStem + "_parameters.json"), 'w', 'n', 'UTF-8');
assert(fid >= 0, 'kde:OutputFile', 'Cannot create parameter file.');
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
parameters = cfg;
% JSON has no infinity literal. Keep the unbounded cutoff explicit in JSON;
% the MAT result above retains the original numeric Inf configuration.
if isfield(parameters,'maxTimeHours') && isinf(parameters.maxTimeHours)
    parameters.maxTimeHours = "Inf";
end
fprintf(fid, '%s\n', jsonencode(parameters, 'PrettyPrint', true));
end
