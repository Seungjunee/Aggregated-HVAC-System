function cases = load_paper_cases(repoRoot, setpointPath, massPath)
%LOAD_PAPER_CASES Load the published setpoint and mass-flow case data.

if nargin < 1 || isempty(repoRoot)
    repoRoot = fileparts(fileparts(mfilename('fullpath')));
end

if nargin < 2 || isempty(setpointPath)
    setpointPath = fullfile(repoRoot, 'data', 'published', ...
        'setpoint_control.mat');
end
if nargin < 3 || isempty(massPath)
    massPath = fullfile(repoRoot, 'data', 'published', ...
        'mass_flow_control.mat');
end

assert(isfile(setpointPath), 'Missing published result file: %s', setpointPath);
assert(isfile(massPath), 'Missing published result file: %s', massPath);

warningState = warning;
warningCleanup = onCleanup(@() warning(warningState));
warning('off', 'all');

cases.setpoint = load(setpointPath, ...
    'a_agg', 'b_agg', 'C_agg', 'hvac', 'P_sched', 'Pnsched', ...
    'Pbaggmax', 'PbaggmaxRelax', 'Pbaggmin', 'PbaggminRelax', ...
    'soc_sched', 'Pref');
cases.mass = load(massPath, ...
    'a_agg', 'b_agg', 'C_agg', 'hvacMass', 'P_sched', ...
    'Pbaggmax', 'PbaggmaxRelax', 'Pbaggmin', 'PbaggminRelax', ...
    'PrefMass', 'soc_sched');
clear warningCleanup;

requiredSetpoint = {'hvac', 'P_sched', 'soc_sched', 'Pbaggmin', 'Pbaggmax'};
requiredMass = {'hvacMass', 'PrefMass'};
assert(all(isfield(cases.setpoint, requiredSetpoint)), ...
    'The setpoint-control data does not contain the expected variables.');
assert(all(isfield(cases.mass, requiredMass)), ...
    'The mass-flow data does not contain the expected variables.');

cases.repoRoot = repoRoot;
cases.setpointPath = setpointPath;
cases.massPath = massPath;
cases.paperCommit = 'c51d05e42cc7446b44e2d43924960ce1953e0be5';
cases.paperFigureExcluded = 4;
end
