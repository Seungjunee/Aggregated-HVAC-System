function [series, metrics, diagnostics] = compute_paper_results(cases)
%COMPUTE_PAPER_RESULTS Recompute the published Figure 5-13 data and metrics.

setpoint = cases.setpoint;
mass = cases.mass;
hvac = setpoint.hvac;
hvacMass = mass.hvacMass;

numBuildings = numel(hvac);
assert(numBuildings == 20 && numel(hvacMass) == 20, ...
    'The published case is expected to contain 20 buildings.');

day = 9;
drStart = 96 * day + 10 * 4;
drEnd = 96 * day + 18 * 4;
temperatureIndex = drStart:drEnd;
powerIndex = (drStart + 1):(drEnd + 1);

series.time = (10:0.25:18)';
series.lambdaTime = (10:0.25:18.25)';
series.powerSignal = setpoint.P_sched{1}(:);
series.lambda = setpoint.soc_sched{1}(:);
series.baseline = sum_baseline_power(hvac);
series.actualPowerSetpoint = sum_actual_power(hvac, powerIndex);
series.actualPowerMass = sum_actual_power(hvacMass, powerIndex);
series.referencePower = series.baseline + series.powerSignal;
series.dischargeLimit = setpoint.Pbaggmin(:);
series.chargeLimit = setpoint.Pbaggmax(:);

representativeBuildings = [1, 11];
for k = 1:numel(representativeBuildings)
    n = representativeBuildings(k);
    label = sprintf('b%d', n);
    series.vb.(label).baseline = hvac(n).Pbase(:);
    series.vb.(label).baselineActual = hvac(n).Pbase_true(temperatureIndex)';
    series.vb.(label).reference = mass.PrefMass(:, n);
    series.vb.(label).actualSetpoint = actual_power(hvac(n), powerIndex);
    series.vb.(label).actualMass = actual_power(hvacMass(n), powerIndex);
    series.temperature.setpoint.(label) = zone_temperature(hvac(n), temperatureIndex);
    series.temperature.mass.(label) = zone_temperature(hvacMass(n), temperatureIndex);
    series.soz.setpoint.(label) = zone_state(hvac(n), temperatureIndex);
    series.soz.mass.(label) = zone_state(hvacMass(n), temperatureIndex);
end

series.socByBuildingSetpoint = building_soc(hvac, temperatureIndex, true);
series.socByBuildingMass = building_soc(hvacMass, temperatureIndex, true);
series.socByBuildingMassUnweighted = building_soc(hvacMass, temperatureIndex, false);
series.socVessSetpoint = aggregate_soc(hvac, series.socByBuildingSetpoint);
series.socVessMass = aggregate_soc(hvacMass, series.socByBuildingMass);

metrics.baselineRmseB1 = rmse(series.vb.b1.baselineActual, series.vb.b1.baseline);
metrics.baselineRmseB11 = rmse(series.vb.b11.baselineActual, series.vb.b11.baseline);

referenceRange = range(series.referencePower);
metrics.vessRmseSetpoint = rmse(series.actualPowerSetpoint, series.referencePower);
metrics.vessRmseMass = rmse(series.actualPowerMass, series.referencePower);
metrics.vessNrmseSetpoint = metrics.vessRmseSetpoint / referenceRange;
metrics.vessNrmseMass = metrics.vessRmseMass / referenceRange;

for n = representativeBuildings
    label = sprintf('b%d', n);
    localRange = range(series.vb.(label).reference);
    metrics.(['vb' num2str(n) 'RmseSetpoint']) = rmse( ...
        series.vb.(label).actualSetpoint, series.vb.(label).reference);
    metrics.(['vb' num2str(n) 'RmseMass']) = rmse( ...
        series.vb.(label).actualMass, series.vb.(label).reference);
    metrics.(['vb' num2str(n) 'NrmseSetpoint']) = ...
        metrics.(['vb' num2str(n) 'RmseSetpoint']) / localRange;
    metrics.(['vb' num2str(n) 'NrmseMass']) = ...
        metrics.(['vb' num2str(n) 'RmseMass']) / localRange;
end

sozSetpoint = [series.soz.setpoint.b1, series.soz.setpoint.b11];
sozMass = [series.soz.mass.b1, series.soz.mass.b11];
metrics.sozVarianceVessSetpoint = mean(var(sozSetpoint, 0, 2));
metrics.sozVarianceThreeZoneSetpoint = mean(var(sozSetpoint(:, 1:3), 0, 2));
metrics.sozVarianceFiveZoneSetpoint = mean(var(sozSetpoint(:, 4:8), 0, 2));
metrics.sozVarianceVessMass = mean(var(sozMass, 0, 2));
metrics.sozVarianceThreeZoneMass = mean(var(sozMass(:, 1:3), 0, 2));
metrics.sozVarianceFiveZoneMass = mean(var(sozMass(:, 4:8), 0, 2));

metrics.socVarianceSetpoint = mean(var(series.socByBuildingSetpoint, 0, 1));
metrics.socVarianceMass = mean(var(series.socByBuildingMass, 0, 1));
metrics.socVessRmseSetpoint = rmse(series.socVessSetpoint, series.lambda(1:end-1));
metrics.socVessRmseMass = rmse(series.socVessMass, series.lambda(1:end-1));

buildingWeights = 1 ./ [hvac.bhat]';
weightedStdSetpoint = sqrt(sum( ...
    (series.socByBuildingSetpoint - series.socVessSetpoint').^2 .* buildingWeights, 1) ...
    / sum(buildingWeights));
weightedStdMass = sqrt(sum( ...
    (series.socByBuildingMass - series.socVessMass').^2 .* buildingWeights, 1) ...
    / sum(buildingWeights));
metrics.stdAlpha = sqrt(sum(([hvac.a]' - setpoint.a_agg).^2 ./ [hvac.bhat]') ...
    / sum(1 ./ [hvac.bhat]'));
metrics.socApproximationBoundSetpoint = mean(weightedStdSetpoint) * metrics.stdAlpha;
metrics.socApproximationBoundMass = mean(weightedStdMass) * metrics.stdAlpha;

metrics.vbParameters = vb_parameters(hvac);
metrics.vessParameters = [setpoint.a_agg, setpoint.b_agg, setpoint.C_agg, ...
    mean(setpoint.Pbaggmin), mean(setpoint.Pbaggmax)];

[correctedDischarge, correctedCharge] = vess_power_limits(hvac, ...
    setpoint.a_agg, setpoint.b_agg);
correctedDischarge = correctedDischarge(:);
correctedCharge = correctedCharge(:);
diagnostics.legacyMeanDischargeLimit = mean(series.dischargeLimit);
diagnostics.legacyMeanChargeLimit = mean(series.chargeLimit);
diagnostics.correctedMeanDischargeLimit = mean(correctedDischarge);
diagnostics.correctedMeanChargeLimit = mean(correctedCharge);
diagnostics.meanDischargeLimitReduction = mean(series.dischargeLimit - correctedDischarge);
diagnostics.meanChargeLimitReduction = mean(series.chargeLimit - correctedCharge);

socDifference = series.socByBuildingMassUnweighted - series.socByBuildingMass;
diagnostics.massSocUnweightedVsWeightedMae = mean(abs(socDifference), 'all');
diagnostics.massSocUnweightedVsWeightedRmse = sqrt(mean(socDifference.^2, 'all'));
diagnostics.massSocUnweightedVsWeightedMaxAbs = max(abs(socDifference), [], 'all');
diagnostics.lambdaMin = min(series.lambda);
diagnostics.lambdaMax = max(series.lambda);
end

function baseline = sum_baseline_power(objects)
baseline = zeros(numel(objects(1).Pbase), 1);
for n = 1:numel(objects)
    baseline = baseline + objects(n).Pbase(:);
end
end

function total = sum_actual_power(objects, index)
total = zeros(numel(index), 1);
for n = 1:numel(objects)
    total = total + actual_power(objects(n), index);
end
end

function power = actual_power(object, index)
mask = contains(object.logTable.Properties.VariableNames, 'Electric_Power');
power = sum(table2array(object.logTable(index, mask)), 2) / 1e3;
end

function temperature = zone_temperature(object, index)
mask = contains(object.logTable.Properties.VariableNames, 'Zone_Air_Temperature');
temperature = table2array(object.logTable(index, mask));
end

function soz = zone_state(object, index)
temperature = zone_temperature(object, index);
soz = (object.Tset(:)' - temperature) ./ object.delta(:)';
end

function soc = building_soc(objects, index, weighted)
soc = zeros(numel(objects), numel(index));
for n = 1:numel(objects)
    soz = zone_state(objects(n), index);
    if weighted
        weights = 1 ./ objects(n).Btild(:);
        weights = weights / sum(weights);
        soc(n, :) = (soz * weights)';
    else
        soc(n, :) = mean(soz, 2)';
    end
end
end

function socVess = aggregate_soc(objects, socByBuilding)
weights = 1 ./ [objects.bhat]';
weights = weights / sum(weights);
socVess = (weights' * socByBuilding)';
end

function value = rmse(actual, expected)
actual = actual(:);
expected = expected(:);
value = sqrt(mean((actual - expected).^2));
end

function parameters = vb_parameters(objects)
parameters = zeros(numel(objects), 5);
for n = 1:numel(objects)
    parameters(n, :) = [objects(n).a, objects(n).bhat, ...
        0.25 / objects(n).bhat, mean(objects(n).Pbmin), mean(objects(n).Pbmax)];
end
end
