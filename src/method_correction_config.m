function config = method_correction_config(caseName)
%METHOD_CORRECTION_CONFIG Define published and A-D correction variants.

if nargin < 1 || isempty(caseName)
    caseName = 'D';
end
caseName = upper(char(caseName));

config.name = caseName;
config.equation30Limits = false;
config.weightedBuildingSoc = false;
config.socLimit = 5;

switch caseName
    case 'PUBLISHED'
        % Implementation used for the saved published cases.
    case 'A'
        config.equation30Limits = true;
    case 'B'
        config.weightedBuildingSoc = true;
    case 'C'
        config.socLimit = 1;
    case 'D'
        config.equation30Limits = true;
        config.weightedBuildingSoc = true;
        config.socLimit = 1;
    otherwise
        error('Unknown method case "%s". Use PUBLISHED, A, B, C, or D.', caseName);
end

config.unitSocBounds = config.socLimit == 1;
end
