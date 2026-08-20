function results = run_method_correction_tests()
%RUN_METHOD_CORRECTION_TESTS Verify the three paper-method corrections.

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repoRoot, 'src'));

cases = load_paper_cases(repoRoot);
caseNames = {'PUBLISHED', 'A', 'B', 'C', 'D'};
expectedConfig = [false false 5; true false 5; false true 5; ...
    false false 1; true true 1];
for k = 1:numel(caseNames)
    config = method_correction_config(caseNames{k});
    actualConfig = [config.equation30Limits, ...
        config.weightedBuildingSoc, config.socLimit];
    assert(isequal(actualConfig, expectedConfig(k, :)), ...
        'Method correction case %s is configured incorrectly.', caseNames{k});
end

[legacyDischarge, legacyCharge] = vess_power_limits(cases.setpoint.hvac, ...
    cases.setpoint.a_agg, cases.setpoint.b_agg, false);
[dischargeLimit, chargeLimit] = vess_power_limits(cases.setpoint.hvac, ...
    cases.setpoint.a_agg, cases.setpoint.b_agg, true);

results.legacyMeanDischargeLimit = mean(legacyDischarge);
results.legacyMeanChargeLimit = mean(legacyCharge);
results.meanDischargeLimit = mean(dischargeLimit);
results.meanChargeLimit = mean(chargeLimit);
assert(abs(results.legacyMeanDischargeLimit - 73.390175652) < 1e-8, ...
    'Published implementation discharge limit changed unexpectedly.');
assert(abs(results.legacyMeanChargeLimit - 25.452018097) < 1e-8, ...
    'Published implementation charge limit changed unexpectedly.');
assert(abs(results.meanDischargeLimit - 73.170965568) < 1e-8, ...
    'Equation (30a) discharge limit does not match the corrected baseline.');
assert(abs(results.meanChargeLimit - 24.908490202) < 1e-8, ...
    'Equation (30b) charge limit does not match the corrected baseline.');

weights = zone_soc_weights([2; 1], true);
uniformWeights = zone_soc_weights([2; 1], false);
results.weightedSoc = aggregate_zone_soc([1; -1], [2; 1], true);
results.unweightedSoc = aggregate_zone_soc([1; -1], [2; 1], false);
assert(max(abs(weights - [1/3; 2/3])) < 1e-12, ...
    'Zone SOC weights do not implement delta_i / b_i.');
assert(max(abs(uniformWeights - [1/2; 1/2])) < 1e-12, ...
    'Published implementation does not reproduce the unweighted mean.');
assert(abs(results.weightedSoc + 1/3) < 1e-12, ...
    'Weighted building SOC is inconsistent with Eqs. (13) and (38).');
assert(abs(results.unweightedSoc) < 1e-12, ...
    'Unweighted building SOC changed unexpectedly.');

[~, publishedSoc] = schedule_demand_response(1, 1, 0, 10 * ones(1, 2), ...
    10 * ones(1, 2), 2, -ones(2, 1), zeros(1, 3), 5);
[~, scheduledSoc] = schedule_demand_response(1, 1, 0, 10 * ones(1, 2), ...
    10 * ones(1, 2), 2, -ones(2, 1), zeros(1, 3), 1);
results.publishedSocMaxAbs = max(abs(publishedSoc));
results.correctedSocMaxAbs = max(abs(scheduledSoc));
assert(results.publishedSocMaxAbs > 1 + 1e-8, ...
    'SOC-bound test case does not distinguish published and corrected limits.');
results.scheduledSocMin = min(scheduledSoc);
results.scheduledSocMax = max(scheduledSoc);
assert(results.scheduledSocMin >= -1 - 1e-8 && ...
    results.scheduledSocMax <= 1 + 1e-8, ...
    'Demand-response scheduling violates the paper SOC bounds.');

fprintf(['PASS: PUBLISHED and A-D configurations, Eq. (30) limits, ' ...
    'weighted building SOC, and scheduling bounds.\n']);
end
