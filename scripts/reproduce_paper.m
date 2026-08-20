function report = reproduce_paper(outputDir)
%REPRODUCE_PAPER Reproduce and validate published Figures 5-13.
% Figure 4 is intentionally excluded from the reproduction workflow.

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repoRoot, 'src'));

if nargin < 1 || isempty(outputDir)
    outputDir = fullfile(repoRoot, 'output', 'paper_reproduction');
end
figureDir = fullfile(outputDir, 'figures');
validationDir = fullfile(outputDir, 'validation');
if ~isfolder(validationDir)
    mkdir(validationDir);
end

cases = load_paper_cases(repoRoot);
[series, metrics, diagnostics] = compute_paper_results(cases);
validation = validate_paper_results(metrics);
figureFiles = plot_paper_figures(series, figureDir);

writetable(validation, fullfile(validationDir, 'paper_baseline_validation.csv'));
save(fullfile(validationDir, 'paper_reproduction_metrics.mat'), ...
    'metrics', 'diagnostics', 'validation');

report.validation = validation;
report.metrics = metrics;
report.diagnostics = diagnostics;
report.figureFiles = figureFiles;
report.outputDir = outputDir;
report.passed = all(validation.passed);

fprintf('Paper baseline validation: %d/%d checks passed.\n', ...
    sum(validation.passed), height(validation));
fprintf('Figure 4 excluded by design. Generated Figures 5-13 in:\n  %s\n', figureDir);
fprintf('Eq. (30) corrected mean limits: discharge %.4f kW, charge %.4f kW.\n', ...
    diagnostics.correctedMeanDischargeLimit, diagnostics.correctedMeanChargeLimit);
fprintf('Mass SOC unweighted-vs-weighted max difference: %.6f.\n', ...
    diagnostics.massSocUnweightedVsWeightedMaxAbs);

if nargout == 0 && ~report.passed
    failed = validation(~validation.passed, :);
    disp(failed);
    error('Published baseline validation failed.');
end
end
