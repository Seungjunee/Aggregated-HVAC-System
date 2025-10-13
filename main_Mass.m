%% Matlab <-> EnergyPlus co-simulation
clc
close all
clear all

%   Construct HVAC systems
day = 9;
NumB1 = 10;
NumB2 = 10;
NumB = NumB1 + NumB2;
DRstart = 10*4;
DRend = 18*4;
DRsize = DRend - DRstart + 1;

for bldg = 1:NumB1
    hvacMass(bldg) = hvac_model_temperature(1/4,bldg,"Threezone_buildings",'Mass');
end
for bldg = (NumB1+1):NumB
    hvacMass(bldg) = hvac_model_temperature(1/4,bldg-NumB1,"Fivezone_buildings",'Mass');
end
b_agg = 1/sum(1./[hvacMass.bhat]);
a_agg = sum([hvacMass.a]./[hvacMass.bhat])*b_agg;
C_agg = 0.25/b_agg;
Pbmaxhvac = [];
Pbminhvac = [];

for bldg = 1:NumB
    Pbminhvac = [Pbminhvac; (hvacMass(bldg).Pbmin)'];
    Pbmaxhvac = [Pbmaxhvac; (hvacMass(bldg).Pbmax)'];
    %     Pbaggmin = Pbaggmin + hvac(bldg).Pbmin;
    %     Pbaggmax = Pbaggmax + hvac(bldg).Pbmax;
end
Pbaggmin = min(Pbminhvac .* [hvacMass.bhat]' / b_agg,[],1);
Pbaggmax = min(Pbmaxhvac .* [hvacMass.bhat]' / b_agg,[],1);
PbaggminRelax = sum(Pbminhvac,2);
PbaggmaxRelax = sum(Pbmaxhvac,2);
bhatTable = round([hvacMass.bhat],2);
CapacityTable = round(0.25./[hvacMass.bhat],2);
alphaTable = round([hvacMass.a],2);
PbminTable = round(mean(Pbminhvac,2),2);
PbmaxTable = round(mean(Pbmaxhvac,2),2);

tauagg = [150*ones(1,16) 100*ones(1,16)];

DR_signal1 = [zeros(1,8), zeros(1,25)];
DR_signal1 = -[-20 -35 -50 -65 -80 -90 -100 -100,-100,-90,-80, -65, -50, -35, -20, -10, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 80, 70, 60, 50, 40, 25, 10]*0.15;
DR_signal1 = -[-20 -35 -50 -65 -75 -85 -90 -95,-90,-85,-80, -65, -50, -35, -20, -5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 90, 80, 70, 55, 40, 25, 10]*0.18;
% [DR_signal1,~] = demand_response(a_agg,b_agg,0,Pbaggmin,Pbaggmax,DRsize,RTP(DRstart:DRend),[tauagg, 150, 150]);
% DR_signal1 = [zeros(1,8), 150*ones(1,25)]*0.1;
DR_signal2 = [zeros(1,10), Pbaggmax(11)*2 Pbaggmax(12)*2 Pbaggmax(13) zeros(1,10), Pbaggmax(24)*2 Pbaggmax(25)*2 Pbaggmax(26)*2 zeros(1,7)];
DR_signal3 = [zeros(1,8), 150*ones(1,25)]*0.1;
DR_signal4 = [-80 -90 -120 -140 -200 -170 -100 -100,-100,-90,-80, -65, -50, -35, -20, -10, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 80, 80, 80, 90, 90, 90, 90]*2;
[P_sched{1},soc_sched{1}] = demand_response(a_agg,b_agg,0,-DR_signal1,DR_signal1,73-40,zeros(33,1),[tauagg, 0, 0]);
% [P_sched{1},soc_sched{1}] = demand_response(a_agg,b_agg,0,-DR_signal2,DR_signal2,73-40,zeros(33,1),[tauagg, 0, 0]);
[P_sched{3},soc_sched{3}] = demand_response(a_agg,b_agg,0,-DR_signal3,DR_signal3,73-40,zeros(33,1),[tauagg, 0, 0]);

for sig = 1
    for bldg = 1:NumB
        hvacMass(bldg).ep.start
    end
    %% Simulate
    % Specify simulation duration
    endTime = 10*24*60*60; %[s]
    iLog = 1;
    t = 0;
    record_time = [];
    % Start the co-simulation process and communication.
    PrefMass = zeros(DRsize,NumB);
    Coileff = 0.68; ChillerCOP = 3.9; Tc = 13; cair = 1.03;
    soc_table = zeros(1,NumB);
    for bldg = 1:NumB
%         Pref(:,bldg) = hvac(bldg).Pbase_true(96*1 + [DRstart:DRend])';
        PrefMass(:,bldg) = hvacMass(bldg).Pbase';
        
    end
    
    while t < endTime-900
        for bldg = 1:NumB
            % Get outputs from EnergyPlus
            [y, t] = hvacMass(bldg).ep.read;
            % Log
            hvacMass(bldg).logTable(iLog, :) = num2cell([t y(:)']);
        end
        for bldg = 1:NumB
            tempEnd = table2array(hvacMass(bldg).logTable(end,contains(hvacMass(bldg).logTable.Properties.VariableNames,'Zone_Air_Temperature')));
            socTable(bldg) = mean((hvacMass(bldg).Tset' - tempEnd)./(hvacMass(bldg).delta'));
        end
        for bldg = 1:NumB
            Tm = table2array(hvacMass(bldg).logTable(end,contains(hvacMass(bldg).logTable.Properties.VariableNames,'Mix')));
            Tc = table2array(hvacMass(bldg).logTable(end,contains(hvacMass(bldg).logTable.Properties.VariableNames,'Coil_1_Outlet_Node')));
%             Tc = 13;
            Coileff = table2array(hvacMass(bldg).logTable(end,contains(hvacMass(bldg).logTable.Properties.VariableNames,'Cooling_Coil_Sensible_Cooling_Rate')))...
                /table2array(hvacMass(bldg).logTable(end,contains(hvacMass(bldg).logTable.Properties.VariableNames,'Chiller_Evaporator_Cooling_Rate')));
            if bldg<11
                Tc = Tc-0.7;
            end
            ChillerCOP = table2array(hvacMass(bldg).logTable(end,contains(hvacMass(bldg).logTable.Properties.VariableNames,'Chiller_COP')));
            Ti = hvacMass(bldg).logTable(end,contains(hvacMass(bldg).logTable.Properties.VariableNames,'Zone_Air_Temperature'));
            Ti = table2array(Ti)';
            soz = (hvacMass(bldg).Tset - Ti)./hvacMass(bldg).delta;
            
            % Prepare inputs (for next step t+1)
            if (t < day*24*60*60 + 3600*9.75) || (t >= day*24*60*60 + 3600*18)
                mdot = hvacMass(bldg).Massflowbase(t/900+1,:);
            if bldg<11
                mdot = mdot;
            end
            else
                PrefMass((t-day*24*60*60)/900-DRstart+2,bldg) = PrefMass((t-day*24*60*60)/900-DRstart+2,bldg) + P_sched{sig}((t-day*24*60*60)/900-DRstart+2) * b_agg/hvacMass(bldg).bhat + (a_agg*sum(socTable./[hvacMass.bhat])/sum(1./[hvacMass.bhat]) - hvacMass(bldg).a*socTable(bldg))/hvacMass(bldg).bhat;
                tempb = cair*(Tm-Tc)/(Coileff*ChillerCOP) + hvacMass(bldg).Fancoeff_2;
                mtot = (-tempb + sqrt(tempb^2-4*hvacMass(bldg).Fancoeff_1*(-PrefMass((t-day*24*60*60)/900-DRstart+2,bldg)+hvacMass(bldg).Fancoeff_3))) / (2*hvacMass(bldg).Fancoeff_1);
                mtot(isnan(mtot))=0;
                mdotvar = optimvar('mdotvar',hvacMass(bldg).numZ,1,'LowerBound',0,'UpperBound',hvacMass(bldg).m_high'*1.2); %
                soz_star = optimvar('soz_star',hvacMass(bldg).numZ,1); soc_star = optimvar('soc_star',1);
                prob = optimproblem;
                prob.Objective = sum((soc_star - soz_star).^2);
                mtot = max(sum(hvacMass(bldg).m_low),min(mtot,sum(hvacMass(bldg).m_high*1.2)));
                prob.Constraints = sum(mdotvar) == mtot;
                for z = 1:hvacMass(bldg).numZ
                    prob.Constraints = [prob.Constraints; soz_star(z) == hvacMass(bldg).A(z,:)*soz + hvacMass(bldg).Btild(z)*(mdotvar(z)*(Ti(z)-Tc)-hvacMass(bldg).qbase(z,t/900 - 96*day - DRstart + 2))];
%                     prob.Constraints = [prob.Constraints; soz_star(z) == hvac(bldg).A(z,:)*soz + hvac(bldg).Btild(z)*(mdotvar(z)*(Ti(z)-Tc)-hvac(bldg).qbase_true(z,t/900+1))];
                end
                prob.Constraints = [prob.Constraints; soc_star == mean(soz_star)];
                tic
                [sol,fval,exitflag,output] = solve(prob);
                record_time = [record_time, toc];
                mdot = sol.mdotvar';
%                 mdot = hvac(bldg).Massflowbase(t/900+1,:);
            end
            % Send inputs to EnergyPlus
            hvacMass(bldg).ep.write(mdot,t);
        end
        iLog = iLog + 1;
    end
    % Stop co-simulation process
    for bldg = 1:NumB
        hvacMass(bldg).ep.stop
    end
end
save(strcat('CaseMass',int2str(sig),'.mat'),'a_agg','b_agg','C_agg','hvacMass','P_sched','Pbaggmax','PbaggmaxRelax','Pbaggmin','PbaggminRelax','PrefMass','soc_sched')
