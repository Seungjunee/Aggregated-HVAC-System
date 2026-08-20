function weights = zone_soc_weights(bTilde, applyThermalWeights)
%ZONE_SOC_WEIGHTS Return normalized delta_i / b_i weights from Eq. (13).

bTilde = bTilde(:);
assert(~isempty(bTilde) && all(bTilde > 0), ...
    'All normalized zone input coefficients must be positive.');
if nargin < 2 || isempty(applyThermalWeights)
    applyThermalWeights = true;
end
if applyThermalWeights
    weights = 1 ./ bTilde;
else
    weights = ones(size(bTilde));
end
weights = weights / sum(weights);
end
