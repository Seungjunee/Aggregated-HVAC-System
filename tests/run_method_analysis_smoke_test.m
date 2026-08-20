function report = run_method_analysis_smoke_test()
%RUN_METHOD_ANALYSIS_SMOKE_TEST Exercise comparison using baseline as candidate.

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repoRoot, 'scripts'));
outputDir = fullfile(repoRoot, 'output', 'test_method_analysis');
report = analyze_method_case('D', ...
    fullfile(repoRoot, 'data', 'published', 'setpoint_control.mat'), ...
    fullfile(repoRoot, 'data', 'published', 'mass_flow_control.mat'), ...
    outputDir);

assert(report.figureExcluded == 4, 'Figure 4 exclusion was not recorded.');
assert(numel(report.figureFiles) == 9, ...
    'Method analysis did not generate Figures 5-13 exactly once.');
assert(report.maxAbsoluteMetricDelta < 1e-12, ...
    'Baseline self-comparison produced a nonzero metric difference.');
assert(~any(contains(report.figureFiles, 'figure_04')), ...
    'Figure 4 was generated unexpectedly.');
fprintf('PASS: method analysis self-comparison; Figure 4 excluded.\n');
end
