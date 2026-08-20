function validation = run_paper_validation()
%RUN_PAPER_VALIDATION Run the numerical regression test without plotting.

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repoRoot, 'src'));

cases = load_paper_cases(repoRoot);
[~, metrics] = compute_paper_results(cases);
validation = validate_paper_results(metrics);

if ~all(validation.passed)
    disp(validation(~validation.passed, :));
end
assert(all(validation.passed), 'Published paper baseline regression failed.');
fprintf('PASS: %d/%d published baseline checks. Figure 4 excluded.\n', ...
    sum(validation.passed), height(validation));
end
