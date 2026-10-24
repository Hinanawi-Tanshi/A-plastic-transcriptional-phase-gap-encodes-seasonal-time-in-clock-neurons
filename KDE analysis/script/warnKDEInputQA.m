function warnKDEInputQA(qa,cfg)
% WARNKDEINPUTQA One summary warning; detailed source rows are in the QA CSV.
backwards = sum(qa.BackwardSteps > 0);
duplicates = sum(qa.DuplicateTimes > 0);
if backwards > 0 || duplicates > 0
    warning('kde:NativeTimestamps', ...
        ['Input QA: %d traces have backwards timestamps; %d have duplicates. ' ...
        'Source values and analysis rules are preserved. Details: %s'], ...
        backwards,duplicates,fullfile(cfg.outputDir,cfg.outputStem + "_input_QA.csv"));
end
end
