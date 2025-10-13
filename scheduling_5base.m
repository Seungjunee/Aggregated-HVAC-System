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
cd Fivezone_buildings

day  = 9;
NumB = 10;
NumZ = 5;
%% Create mlep instance and configure it
for bldg=1:NumB
load(strcat('coefficients',int2str(bldg),'.mat'))

ep{bldg} = mlep;
ep{bldg}.idfFile = strcat('BuildingTemperature',int2str(bldg));
ep{bldg}.epwFile = 'USA_FL_Miami.Intl.AP.722020_TMY3';
ep{bldg}.outputDirName = strcat('5zone_base',int2str(bldg));
ep{bldg}.initialize;

%% Simulate

% Specify simulation duration
endTime = 10*24*60*60; %[s]

% Prepare data logging
nRows = ceil(endTime / ep{bldg}.timestep); %Query timestep after mlep initialization
logTable = table('Size',[0, 1 + ep{bldg}.nOut],...
    'VariableTypes',repmat({'double'},1,1 + ep{bldg}.nOut),...
    'VariableNames',[{'Time'}; ep{bldg}.outputSigName]);
iLog = 1;

% Start the co-simulation process and communication. 
ep{bldg}.start

% The simulation loop
t = 0;
while t < endTime
    % Prepare inputs (possibly from last outputs)
    u = Tset';
    % Get outputs from EnergyPlus
    [y, t] = ep{bldg}.read;
    % Send inputs to EnergyPlus
    ep{bldg}.write(u,t); 
    % Log
    logTable(iLog, :) = num2cell([t y(:)']);
    iLog = iLog + 1;
end
% Stop co-simulation process
ep{bldg}.stop;

Massflowbase = table2array(logTable(2:end,2:6));
Chillerbase = table2array(logTable(2:end,19))/1e3;
Fanbase = table2array(logTable(2:end,7))/1e3;
qbase_true = table2array(logTable(2:end,[8,14:17]))'/1e3;

Tohist = table2array(logTable(96*day+40:96*day+73,18));
qbase = zeros(NumZ,DRsize-1);
for t=1:(DRsize-1)
    qbase(:,t) = -inv(diag(valB))*((eye(NumZ)-valA)*Tset - valC*Tohist(t+1)-valD(:,t));
end
qbmax = qmax - qbase';
qbmin = qbase' - qmin;

Qbmax = sum(qmax) - sum(qbase);
Qbmin = sum(qbase) - sum(qmin);

Qbase = sum(qbase);
Pbase = coeff_1 *sum(qbase) + coeff_2; %
Pbmax = coeff_1 *Qbmax';
Pbmin = coeff_1 *Qbmin';

%% Plot results
figure(1)
plot(seconds(table2array(logTable(:,1))),table2array(logTable(:,9:13)));
xtickformat('hh:mm:ss');
legend(logTable.Properties.VariableNames(9:13),'Interpreter','none');
xlabel('Time [hh:mm:ss]');
ylabel('Temperature [C]');

load(strcat('coefficients',int2str(bldg),'.mat'))

figure(2)
Pbase_true = Chillerbase+Fanbase;
plot(0.25:0.25:24,Pbase_true(96*day+1:96*(day+1)),'b','LineWidth',1.5);
hold on
plot(10:0.25:18,Pbase,'r--','LineWidth',1.5);
hold off
xlim([10 18])
xlabel('Time [h]');
ylabel('Baseline Power [kW]');
legend(["Eplus","Estimated"],'Interpreter','latex','fontname','Times New Roman');
xlim([10 18])
ylim([0 8])
set(gca,'FontSize',18)
% 
figure(3)
plot(0.25:0.25:24,sum(qbase_true(:,96*day+1:96*(day+1))));
hold on
plot(10:0.25:18,sum(qbase));
hold off
xlim([10 18])
ylim([0 12])
xlabel('Time');
ylabel('Cooling rate [kW]');
legend(["Eplus","Estimated"],'Interpreter','none');
% 
% figure(4)
% newcolors = [0.83 0.14 0.14
%              1.00 0.54 0.00
%              0.47 0.25 0.80
%              0.25 0.80 0.54
%              0 0 0];
% plot(0.25:0.25:24,qbase_true,'--');
% colororder(newcolors)
% hold on
% plot(10:0.25:18,qbase);
% colororder(newcolors)
% hold off
% xlim([10 18])
% ylim([0 12])
% xlabel('Time');
% ylabel('Cooling rate [kW]');
% legend(["Zone1","Zone2","Zone3","Zone4","Zone5"],'Interpreter','none');

delete(strcat('5zone_base',int2str(bldg),'/*'));
save(strcat('Baseline_info',int2str(bldg),'.mat'), 'Massflowbase','qbase_true','Pbase_true','Pbase','Pbmax','Pbmin','Qbase','qbase');
end