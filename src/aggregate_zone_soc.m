function soc = aggregate_zone_soc(soz, bTilde, applyThermalWeights)
%AGGREGATE_ZONE_SOC Aggregate zone SOZ using selected experiment weights.

if nargin < 3 || isempty(applyThermalWeights)
    applyThermalWeights = true;
end
soz = soz(:);
weights = zone_soc_weights(bTilde, applyThermalWeights);
assert(numel(soz) == numel(weights), ...
    'SOZ values and normalized zone coefficients must have the same size.');
soc = weights' * soz;
end
