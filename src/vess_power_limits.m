function [dischargeLimit, chargeLimit] = vess_power_limits(objects, alphaVess, betaVess, applyCorrection)
%VESS_POWER_LIMITS Compute the state-aware inner limits in Eq. (30).

assert(~isempty(objects), 'At least one virtual building is required.');
assert(betaVess > 0, 'The aggregate beta parameter must be positive.');
if nargin < 4 || isempty(applyCorrection)
    applyCorrection = true;
end

numBuildings = numel(objects);
numSteps = numel(objects(1).Pbmin);
dischargeCandidates = zeros(numBuildings, numSteps);
chargeCandidates = zeros(numBuildings, numSteps);

for n = 1:numBuildings
    assert(numel(objects(n).Pbmin) == numSteps && ...
        numel(objects(n).Pbmax) == numSteps, ...
        'All virtual-building power limits must have the same horizon.');
    correction = applyCorrection * abs(alphaVess - objects(n).a) / betaVess;
    scale = objects(n).bhat / betaVess;
    dischargeCandidates(n, :) = objects(n).Pbmin(:)' * scale - correction;
    chargeCandidates(n, :) = objects(n).Pbmax(:)' * scale - correction;
end

dischargeLimit = min(dischargeCandidates, [], 1);
chargeLimit = min(chargeCandidates, [], 1);
end
