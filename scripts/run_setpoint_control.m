%% Matlab <-> EnergyPlus co-simulation
clc
close all
clearvars -except methodCase modelRoot

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repoRoot, 'src'));
if ~exist('methodCase', 'var')
    methodCase = 'D';
end
methodCorrections = method_correction_config(methodCase);
if ~exist('modelRoot', 'var')
    modelRoot = repoRoot;
end
threeZoneDir = fullfile(modelRoot, 'Threezone_buildings');
fiveZoneDir = fullfile(modelRoot, 'Fivezone_buildings');
assert(isfolder(threeZoneDir), 'Missing model directory: %s', threeZoneDir);
assert(isfolder(fiveZoneDir), 'Missing model directory: %s', fiveZoneDir);

%   Construct HVAC systems
day = 9;
NumB1 = 10;
NumB2 = 10;
NumB = NumB1 + NumB2;
DRstart = 10*4;
DRend = 18*4;
DRsize = DRend - DRstart + 1;

for bldg = 1:NumB1
    hvac(bldg) = hvac_model_temperature(1/4,bldg,threeZoneDir,'Temperature'); %#ok<SAGROW>
end
for bldg = (NumB1+1):NumB
    hvac(bldg) = hvac_model_temperature(1/4,bldg-NumB1,fiveZoneDir,'Temperature');
end
b_agg = 1/sum(1./[hvac.bhat]);
a_agg = sum([hvac.a]./[hvac.bhat])*b_agg;
C_agg = 0.25/b_agg;
Std_alpha = sqrt(sum(([hvac.a]-a_agg).^2./[hvac.bhat])/sum(1./[hvac.bhat]));

numPowerSteps = numel(hvac(1).Pbmin);
Pbmaxhvac = zeros(NumB, numPowerSteps);
Pbminhvac = zeros(NumB, numPowerSteps);

for bldg = 1:NumB
    Pbminhvac(bldg, :) = hvac(bldg).Pbmin(:)';
    Pbmaxhvac(bldg, :) = hvac(bldg).Pbmax(:)';
end
[Pbaggmin, Pbaggmax] = vess_power_limits(hvac, a_agg, b_agg, ...
    methodCorrections.equation30Limits);
PbaggminRelax = sum(Pbminhvac,2);
PbaggmaxRelax = sum(Pbmaxhvac,2);
bhatTable = round([hvac.bhat],2);
CapacityTable = round(0.25./[hvac.bhat],2);
alphaTable = round([hvac.a],2);
PbminTable = round(mean(Pbminhvac,2),2);
PbmaxTable = round(mean(Pbmaxhvac,2),2);

tauagg = [150*ones(1,16) 100*ones(1,18)]*5;

DR_signal1 = -[-20 -35 -50 -65 -75 -85 -90 -95,-90,-85,-80, -65, -50, -35, -20, -5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 90, 80, 70, 55, 40, 25, 10]*0.18;
% [DR_signal1,soc_sched_act] = schedule_demand_response(a_agg,b_agg,0,Pbaggmin,Pbaggmax,DRsize,RTP(40:72),tauagg);
% DR_signal1 = [zeros(1,8), 150*ones(1,25)]*0.1;
DR_signal2 = [zeros(1,10), Pbaggmax(11)*2 Pbaggmax(12)*2 Pbaggmax(13) zeros(1,10), Pbaggmax(24)*2 Pbaggmax(25)*2 Pbaggmax(26)*2 zeros(1,7)];
DR_signal3 = [zeros(1,8), 150*ones(1,25)]*0.1;
DR_signal4 = [-80 -90 -120 -140 -200 -170 -100 -100,-100,-90,-80, -65, -50, -35, -20, -10, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 80, 80, 80, 90, 90, 90, 90]*2;
[P_sched{1},soc_sched{1}] = schedule_demand_response(a_agg,b_agg,0,-DR_signal1,DR_signal1,73-40,zeros(33,1),tauagg,methodCorrections.socLimit);
[P_sched{3},soc_sched{3}] = schedule_demand_response(a_agg,b_agg,0,-DR_signal3,DR_signal3,73-40,zeros(33,1),tauagg,methodCorrections.socLimit);

for sig = 1
    for bldg = 1:NumB
        hvac(bldg).ep.start
    end
    %% Simulate
    % Specify simulation duration
    endTime = 10*24*60*60; %[s]
    iLog = 1;
    DRlog = 1;
    t = 0;
    % Start the co-simulation process and communication.
    Pref = zeros(DRsize,NumB);
    Pnsched = zeros(DRsize,NumB);
    Coileff = 0.68;
    ChillerCOP = 3.9;
    Tc = 13;
    cair = 1.03;
    soc_table = zeros(1,NumB);
    soc_agg = 0;
    hist_soc_agg = zeros(1,32);
    for bldg = 1:NumB
%         Pref(:,bldg) = hvac(bldg).Pbase_true(96*day + [DRstart:DRend])';
        Pref(:,bldg) = hvac(bldg).Pbase';
        Pref(:,bldg) = Pref(:,bldg) + P_sched{sig}' * b_agg/hvac(bldg).bhat;
        Pnsched(:,bldg) = P_sched{sig}' * b_agg/hvac(bldg).bhat;
    end
    
    hist_mdot = [];
    socTable = zeros(1,NumB);
    kk = 1;
    
    while t < endTime-900
        for bldg = 1:NumB
            % Get outputs from EnergyPlus
            [y, t] = hvac(bldg).ep.read;
            % Log
            hvac(bldg).logTable(iLog, :) = num2cell([t y(:)']);
        end
        for bldg = 1:NumB
            tempEnd = table2array(hvac(bldg).logTable(end,contains(hvac(bldg).logTable.Properties.VariableNames,'Zone_Air_Temperature')));
            sozMeasured = (hvac(bldg).Tset' - tempEnd) ./ hvac(bldg).delta';
            socTable(bldg) = aggregate_zone_soc(sozMeasured, ...
                hvac(bldg).Btild, methodCorrections.weightedBuildingSoc);
        end
        for bldg = 1:NumB
            % Prepare inputs (for next step t+1)
            if (t < day*24*60*60 + 3600*9.75) || (t >= day*24*60*60 + 3600*18)
                u = hvac(bldg).Tset;
                hvac(bldg).ep.write(u,t);
            else
%                 u = hvac(bldg).Tset - hvac(bldg).delta*(a_agg * mean(socTable(bldg)) + b_agg * P_sched{sig}((t-day*24*60*60)/900-DRstart+2));
                u = hvac(bldg).Tset - hvac(bldg).delta*soc_sched{sig}((t-day*24*60*60)/900-DRstart+3);
                hvac(bldg).ep.write(u,t);
            end
        end
        kk = kk + 1;
        iLog = iLog + 1;
    end
    % Stop co-simulation process
    for bldg = 1:NumB
        hvac(bldg).ep.stop
    end
end
resultDir = fullfile(repoRoot, 'output', 'method_corrections');
if ~isfolder(resultDir)
    mkdir(resultDir);
end
save(fullfile(resultDir, sprintf('setpoint_control_%s.mat', ...
    lower(methodCorrections.name))), ...
    'a_agg', 'b_agg', 'C_agg', 'hvac', 'P_sched', 'Pnsched', ...
    'Pbaggmax', 'PbaggmaxRelax', 'Pbaggmin', 'PbaggminRelax', ...
    'soc_sched', 'Pref', 'methodCorrections')



