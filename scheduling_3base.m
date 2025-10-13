%% Simple Matlab <-> EnergyPlus co-simulation example
% Demonstrates the functionality of the mlep (MatLab-EnergyPlus) tool in 
% a small office building simulation scenario.
%
% Note that a start of the simulation period as well as a timestep and
% an input/output configuration is defined by the the EnergyPlus simulation
% configuration file (.IDF). Climatic conditions are obtained from a
% EnergyPlus Weather data file (.EPW). 
%
% See also: mlepMatlab_so_example.m, mlepSimulink_example.slx
clc
close all
clear all

set(0, 'DefaultAxesFontName', 'Times');
cd Threezone_buildings

day = 9;
NumB = 10;
NumZ = 3;
%% Create mlep instance and configure it
for buildidx=1:NumB
load(strcat('coefficients',int2str(buildidx),'.mat'))

ep{buildidx} = mlep;
ep{buildidx}.idfFile = strcat('BuildingTemperature',int2str(buildidx));
ep{buildidx}.epwFile = 'USA_FL_Miami.Intl.AP.722020_TMY3';
ep{buildidx}.outputDirName = strcat('3zone_base',int2str(buildidx));
ep{buildidx}.initialize;

%% Simulate

% Specify simulation duration
endTime = 10*24*60*60; %[s]

% Prepare data logging
nRows = ceil(endTime / ep{buildidx}.timestep); %Query timestep after mlep initialization
logTable = table('Size',[0, 1 + ep{buildidx}.nOut],...
    'VariableTypes',repmat({'double'},1,1 + ep{buildidx}.nOut),...
    'VariableNames',[{'Time'}; ep{buildidx}.outputSigName]);
iLog = 1;

% Start the co-simulation process and communication. 
ep{buildidx}.start

% The simulation loop
t = 0;
while t < endTime
    % Prepare inputs (possibly from last outputs)
    u = Tset';
    % Get outputs from EnergyPlus
    [y, t] = ep{buildidx}.read;
    % Send inputs to EnergyPlus
    ep{buildidx}.write(u,t); 
    % Log
    logTable(iLog, :) = num2cell([t y(:)']);
    iLog = iLog + 1;
end
% Stop co-simulation process
ep{buildidx}.stop;

Massflowbase = table2array(logTable(2:end,contains(logTable.Properties.VariableNames,'Mass_Flow_Rate')));
Pbase_true = sum(table2array(logTable(2:end,contains(logTable.Properties.VariableNames,'Electric_Power'))),2)/1e3;
qbase_true = table2array(logTable(2:end,contains(logTable.Properties.VariableNames,'Zone_Air_System_Sensible_Cooling_Rate')))'/1e3;

qbase = zeros(NumZ,DRsize-1);
Tohist = table2array(logTable(96*day+40:96*day+73,contains(logTable.Properties.VariableNames,'Drybulb_Temperature')));
for t=1:DRsize-1
    qbase(:,t) = -inv(diag(valB))*((eye(NumZ)-valA)*Tset - valC*Tohist(t+1)-valD(:,t));
end
qbmax = qmax - qbase';
qbmin = qbase' - qmin;

Qbmax = sum(qmax) - sum(qbase);
Qbmin = sum(qbase) - sum(qmin);

Qbase = sum(qbase);
Pbase = coeff_1 *sum(qbase) + coeff_2;
Pbmax = coeff_1 *Qbmax';
Pbmin = coeff_1 *Qbmin';
%% Plot results
figure(1)
Tistep = table2array(logTable(:,contains(logTable.Properties.VariableNames,'Zone_Air_Temperature')));
plot(seconds(table2array(logTable(:,1))),Tistep);
xtickformat('hh:mm:ss');
legend(logTable.Properties.VariableNames(contains(logTable.Properties.VariableNames,'Zone_Air_Temperature')),'Interpreter','none');
xlabel('Time [hh:mm:ss]');
ylabel('Temperature [C]');

figure(3)
plot(0.25:0.25:24,sum(qbase_true(:,96*day+1:96*(day+1))));
hold on
plot(10:0.25:18,sum(qbase));
hold off
xlim([10 18])
% ylim([0 12])
xlabel('Time');
ylabel('Cooling rate [kW]');
legend(["Eplus","Estimated"],'Interpreter','none');

figure(4)
plot(0.25:0.25:24,Pbase_true(96*day+[1:96]));
hold on
plot(10:0.25:18,Pbase);
hold off
xlim([10 18])
xlabel('Time');
ylabel('Cooling rate [kW]');
legend(["Eplus","Estimated"],'Interpreter','none');

delete(strcat('3zone_base',int2str(buildidx),'/*'));
save(strcat('Baseline_info',int2str(buildidx),'.mat'), 'Massflowbase','qbase_true','Pbase_true','Pbase','Pbmax','Pbmin','Qbase','qbase');
end
