function report = analyze_method_case(caseName, setpointPath, massPath, outputDir)
%ANALYZE_METHOD_CASE Compare a corrected case with the published baseline.
% Figures 5-13 are generated; Figure 4 is intentionally excluded.

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repoRoot, 'src'));
config = method_correction_config(caseName);
caseLabel = lower(config.name);
methodDir = fullfile(repoRoot, 'output', 'method_corrections');

if nargin < 2 || isempty(setpointPath)
    setpointPath = fullfile(methodDir, ...
        sprintf('setpoint_control_%s.mat', caseLabel));
end
if nargin < 3 || isempty(massPath)
    massPath = fullfile(methodDir, ...
        sprintf('mass_flow_control_%s.mat', caseLabel));
end
if nargin < 4 || isempty(outputDir)
    outputDir = fullfile(methodDir, ['analysis_' caseLabel]);
end

baselineCases = load_paper_cases(repoRoot);
candidateCases = load_paper_cases(repoRoot, setpointPath, massPath);
[~, baselineMetrics] = compute_paper_results(baselineCases);
[series, candidateMetrics, candidateDiagnostics] = ...
    compute_paper_results(candidateCases);

comparison = compare_scalar_metrics(baselineMetrics, candidateMetrics);
figureDir = fullfile(outputDir, 'figures');
if ~isfolder(outputDir)
    mkdir(outputDir);
end
figureFiles = plot_paper_figures(series, figureDir);
writetable(comparison, fullfile(outputDir, 'metric_comparison.csv'));
save(fullfile(outputDir, 'method_case_analysis.mat'), 'config', ...
    'baselineMetrics', 'candidateMetrics', 'candidateDiagnostics', ...
    'comparison', 'setpointPath', 'massPath');

report.config = config;
report.comparison = comparison;
report.figureFiles = figureFiles;
report.outputDir = outputDir;
report.figureExcluded = 4;
report.maxAbsoluteMetricDelta = max(abs(comparison.delta));

fprintf('Method case %s: compared %d scalar metrics.\n', ...
    config.name, height(comparison));
fprintf('Figure 4 excluded. Generated Figures 5-13 in:\n  %s\n', figureDir);
end

function comparison = compare_scalar_metrics(baselineMetrics, candidateMetrics)
names = fieldnames(baselineMetrics);
keep = false(size(names));
for k = 1:numel(names)
    baselineValue = baselineMetrics.(names{k});
    candidateValue = candidateMetrics.(names{k});
    keep(k) = isnumeric(baselineValue) && isscalar(baselineValue) && ...
        isnumeric(candidateValue) && isscalar(candidateValue);
end
names = names(keep);
baseline = zeros(numel(names), 1);
candidate = zeros(numel(names), 1);
for k = 1:numel(names)
    baseline(k) = baselineMetrics.(names{k});
    candidate(k) = candidateMetrics.(names{k});
end
delta = candidate - baseline;
relativeDelta = delta ./ max(abs(baseline), eps);
comparison = table(names, baseline, candidate, delta, relativeDelta, ...
    'VariableNames', {'metric', 'baseline', 'candidate', 'delta', ...
    'relativeDelta'});
end
