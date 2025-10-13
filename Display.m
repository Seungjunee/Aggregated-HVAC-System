
warning('off', 'all');
load('CaseMass1.mat');
load('Case1.mat'); sig = 1;
day = 9;
NumB = 20;
DRstart = 96*day + 10*4;
DRend = 96*day + 18*4;
DRsize = DRend - DRstart + 1;

Tistep5 = [hvac(11).logTable.EP_SPACE1_1__Zone_Air_Temperature, hvac(11).logTable.EP_SPACE2_1__Zone_Air_Temperature, ...
    hvac(11).logTable.EP_SPACE3_1__Zone_Air_Temperature, hvac(11).logTable.EP_SPACE4_1__Zone_Air_Temperature, ...
    hvac(11).logTable.EP_SPACE5_1__Zone_Air_Temperature];
Tistep3 = table2array(hvac(1).logTable(:,contains(hvac(1).logTable.Properties.VariableNames,'Zone_Air_Temperature')));
Paggbase = 0;
for bldg = 1:NumB
    Paggbase = Paggbase + hvac(bldg).Pbase';
end
PaggbaseTrue = 0;
for bldg = 1:NumB
    PaggbaseTrue = PaggbaseTrue + hvac(bldg).Pbase_true(DRstart:DRend);
end

f34 = figure(34);
p1 = plot(10:0.25:18,hvac(1).Pbase_true(DRstart:DRend),'-',"LineWidth",3,'Color',[0.85,0.33,0.10],'MarkerSize',3);
hold on
p2 = plot(10:0.25:18,hvac(11).Pbase_true(DRstart:DRend),'-',"LineWidth",3,'Color',[0.00,0.45,0.75],'MarkerSize',3);
p5 = plot(10:0.25:18,hvac(1).Pbase,':',"LineWidth",3,'Color',[0.85,0.33,0.10],'MarkerSize',3);
p6 = plot(10:0.25:18,hvac(11).Pbase,':',"LineWidth",3,'Color',[0.00,0.45,0.75],'MarkerSize',3);
hold off
set(gca,'FontSize',12,'fontname','Times New Roman')
legend([p1,p2,p5,p6],"$P_{1}^\mathrm{base}$ (Eplus)","$P_{11}^\mathrm{base}$ (Eplus)","$P_{1}^\mathrm{base}$ (predict)","$P_{11}^\mathrm{base}$ (predict)",...
    'Interpreter','latex','fontname','Times New Roman','FontSize',14,'NumColumns',2,'Location','SouthEast')
xlabel('Time [h]','Interpreter','latex','fontname','Times New Roman')
ylabel(['Baseline  Power [kW]'])
xticks([0:1:24])
yticks([0:1:10])
set(gcf,'position',[0,0,600,230])

f3 = figure(35);
s1 = plot(10:0.25:18,Tistep3(DRstart:DRend,1),'-',"LineWidth",3,'Color',[1.00,0.41,0.16],'MarkerSize',3);
hold on
s2 = plot(10:0.25:18,Tistep3(DRstart:DRend,2),'-',"LineWidth",3,'Color',[0.93,0.69,0.13],'MarkerSize',3);
s3 = plot(10:0.25:18,Tistep3(DRstart:DRend,3),'-',"LineWidth",3,'Color',[0.80,0.50,0.10],'MarkerSize',3);
p1 = plot(10:0.25:18,Tistep5(DRstart:DRend,1),'-',"LineWidth",3,'Color',[0.00,0.45,0.75],'MarkerSize',3);
p2 = plot(10:0.25:18,Tistep5(DRstart:DRend,2),'-',"LineWidth",3,'Color',[0.00,0.70,0.70],'MarkerSize',3);
p3 = plot(10:0.25:18,Tistep5(DRstart:DRend,3),'-',"LineWidth",3,'Color',[0.00,0.60,1.00],'MarkerSize',3);
p4 = plot(10:0.25:18,Tistep5(DRstart:DRend,4),'-',"LineWidth",3,'Color',[0.00,0.80,1.00],'MarkerSize',3);
p5 = plot(10:0.25:18,Tistep5(DRstart:DRend,5),'-',"LineWidth",3,'Color',[0.00,0.00,1.00],'MarkerSize',3);
hold off
legend([s1,s2,s3,p1,p2,p3,p4,p5],"$T_{1,1}$","$T_{2,1}$","$T_{3,1}$","$T_{1,11}$","$T_{2,11}$","$T_{3,11}$","$T_{4,11}$","$T_{5,11}$",'Interpreter','latex','fontname','Times New Roman','FontSize',13,'Orientation','horizontal','NumColumns',3,'Location','SouthEast')
set(gca,'FontSize',12,'fontname','Times New Roman')
xlim([10 18])
ylim([20 30])
xticks([0:1:24])
xlabel('Time [h]','Interpreter','latex','fontname','Times New Roman')
ylabel(['Temperature [' char(176) 'C]'])
set(gcf,'position',[10,10,600,230])

f1 = figure(1);
l3 = plot(10:0.25:18,-Pbaggmin,'Color',[0.6,0.6,0.6],'LineWidth',0.5);
hold on
plot(10:0.25:18,Pbaggmax,'Color',[0.6,0.6,0.6],'LineWidth',0.5);
base1 = plot(10:0.25:18,Paggbase,'-','Color',[0,0,0],"LineWidth",3,'MarkerSize',3);
EplusPagg = 0;
EplusQagg = 0;
for bldg = 1:NumB
    Ptable = table2array(hvac(bldg).logTable((DRstart+1):(DRend+1),contains(hvac(bldg).logTable.Properties.VariableNames,'Electric_Power')));
    EplusPagg = EplusPagg + sum(Ptable,2)/1e3;
    Qtable = table2array(hvacMass(bldg).logTable((DRstart+1):(DRend+1),contains(hvacMass(bldg).logTable.Properties.VariableNames,'Electric_Power')));
    EplusQagg = EplusQagg + sum(Qtable,2)/1e3;
end
Eplus = plot(10:0.25:18,EplusPagg,'-','Color',[0.00,0.45,0.74],"LineWidth",3,'MarkerSize',3);
EplusMass = plot(10:0.25:18,EplusQagg,'-','Color',[0.93,0.69,0.13],"LineWidth",3,'MarkerSize',3);
p2 = plot(10:0.25:18,Paggbase+P_sched{sig}','--r',"LineWidth",2,'MarkerSize',3);
p1 = plot(10:0.25:18,P_sched{sig}','-','Color',[0.58,0.26,0.96],'LineWidth',3,'MarkerSize',3);
RMSEbaseline1 = sqrt(mean((hvac(1).Pbase_true(DRstart:DRend)'-hvac(1).Pbase).^2));
RMSEbaseline11 = sqrt(mean((hvac(11).Pbase_true(DRstart:DRend)'-hvac(11).Pbase).^2));
fprintf('Baseline RMSE for building 1: %f [kW]\n', RMSEbaseline1);
fprintf('Baseline RMSE for building 11: %f [kW]\n', RMSEbaseline11);
hold off
xlabel('Time [h]','Interpreter','latex','fontname','Times New Roman')
ylabel('Power [kW]','Interpreter','latex','fontname','Times New Roman')
set(gca,'FontSize',12,'fontname','Times New Roman')
legend([p1,base1,p2,Eplus,EplusMass,l3],["$P^{\mathrm{ch/dch}}_\mathrm{VESS}$","$P^{\mathrm{baseline}}_\mathrm{VESS}$","$P^{\mathrm{ref}}_\mathrm{VESS}$","$P^{\mathrm{setpoint}}_\mathrm{VESS}$","$P^{\mathrm{mass}}_\mathrm{VESS}$","$\mathrm{Limits}$"],'Interpreter','latex','fontname','Times New Roman','NumColumns',3,'FontSize',14,'Location','southoutside'); % ,
set(gcf,'position',[10,10,600,340])
msetemp = mean((Paggbase+P_sched{sig}'-EplusPagg).^2);
rmsetemp = sqrt(msetemp);                % Root Mean Squared Error
msemass = mean((Paggbase+P_sched{sig}'-EplusQagg).^2);
rmsemass = sqrt(msemass);                % Root Mean Squared Error
range_y = max(Paggbase+P_sched{sig}') - min(Paggbase+P_sched{sig}');
nrmsetemp = rmsetemp / range_y; % NRMSE
nrmsemass = rmsemass / range_y; % NRMSE
fprintf('Power tracking RMSE for all buildings (setpoint control): %f\n', rmsetemp);
fprintf('Power tracking RMSE for all buildings (mass flow rate control): %f\n', rmsemass);
fprintf('Power tracking nRMSE for all buildings (setpoint control): %f\n', nrmsetemp);
fprintf('Power tracking nRMSE for all buildings (mass flow rate control): %f\n', nrmsemass);
ylim([-100 200])

f111 = figure(111);
base1 = plot(10:0.25:18,Paggbase,'-o','Color',[0,0,0],"LineWidth",3,'MarkerSize',3);
hold on
plot(10:0.25:18,Paggbase'-Pbaggmin,'Color',[0.6,0.6,0.6],'LineWidth',0.5);
plot(10:0.25:18,Paggbase'+Pbaggmax,'Color',[0.6,0.6,0.6],'LineWidth',0.5);
EplusPagg = 0;
EplusQagg = 0;
for bldg = 1:NumB
    Ptable = table2array(hvacMass(bldg).logTable((DRstart+1):(DRend+1),contains(hvacMass(bldg).logTable.Properties.VariableNames,'Electric_Power')));
    EplusPagg = EplusPagg + sum(Ptable,2)/1e3;
    Qtable = table2array(hvacMass(bldg).logTable((DRstart+1):(DRend+1),contains(hvacMass(bldg).logTable.Properties.VariableNames,'Sensible_Cooling_Rate')));
    EplusQagg = EplusQagg + sum(Qtable,2)/1e3;
end
EplusCOP = EplusQagg./EplusPagg;
Eplus = plot(10:0.25:18,EplusPagg,'-o','Color',[0.00,0.45,0.74],"LineWidth",3,'MarkerSize',3);
p2 = plot(10:0.25:18,Paggbase+P_sched{sig}','--ro',"LineWidth",2,'MarkerSize',3);
p1 = plot(10:0.25:18,P_sched{sig}','-o','Color',[0.8,0.2,1],'LineWidth',3,'MarkerSize',3);
hold off
xlabel('Time [h]','Interpreter','latex','fontname','Times New Roman')
ylabel('Power [kW]','Interpreter','latex','fontname','Times New Roman')
set(gca,'FontSize',12,'fontname','Times New Roman')
legend([base1,p2,Eplus,p1],["Baseline","$P_{ref}$","Actual (Eplus)","DR signal"],'Interpreter','latex','fontname','Times New Roman','NumColumns',2,'FontSize',12); % ,
set(gcf,'position',[10,10,600,230])
msetemp = mean((Paggbase+P_sched{sig}'-EplusPagg).^2);
rmse = sqrt(msetemp);                % Root Mean Squared Error
range_y = max(Paggbase+P_sched{sig}') - min(Paggbase+P_sched{sig}'); % 실제 값의 범위


f2 = figure(2);
bldg = 11;
base1 = plot(10:0.25:18,hvac(bldg).Pbase,'-','Color',[0,0,0],"LineWidth",3,'MarkerSize',3);
hold on
EplusPagg = 0;
EplusQagg = 0;
Ptable = table2array(hvac(bldg).logTable((DRstart+1):(DRend+1),contains(hvac(bldg).logTable.Properties.VariableNames,'Electric_Power')));
EplusPagg = EplusPagg + sum(Ptable,2)/1e3;
Qtable = table2array(hvacMass(bldg).logTable((DRstart+1):(DRend+1),contains(hvacMass(bldg).logTable.Properties.VariableNames,'Electric_Power')));
EplusQagg = EplusQagg + sum(Qtable,2)/1e3;
Eplus = plot(10:0.25:18,EplusPagg,'-','Color',[0.00,0.45,0.74],"LineWidth",3,'MarkerSize',3);
p1 = plot(10:0.25:18,EplusQagg,'-','Color',[0.93,0.69,0.13],'LineWidth',3,'MarkerSize',3);
p2 = plot(10:0.25:18,PrefMass(:,bldg),'--r',"LineWidth",2,'MarkerSize',3);
hold off
xlabel('Time [h]','Interpreter','latex','fontname','Times New Roman')
ylabel('Power [kW]','Interpreter','latex','fontname','Times New Roman')
set(gca,'FontSize',12,'fontname','Times New Roman')
legend([base1,p2,Eplus,p1],["$P^{\mathrm{baseline}}_{11}$","$P^{\mathrm{ref}}_{11}$","$P^{\mathrm{setpoint}}_{11}$","$P^{\mathrm{mass}}_{11}$"],'Interpreter','latex','fontname','Times New Roman','NumColumns',2,'FontSize',12,'Location','NorthEast'); % ,
set(gcf,'position',[10,10,600,230])
xlim([10 18])
ylim([5 9])
msetemp = mean((EplusPagg-PrefMass(:,bldg)).^2,'all');

range_y5 = max(PrefMass(:,bldg)) - min(PrefMass(:,bldg)); % 실제 값의 범위
rmse5 = sqrt(msetemp);
nrmse5 = rmse5/range_y5;
rmseMass5 = sqrt(mean((EplusQagg-PrefMass(:,bldg)).^2,'all'));
nrmseMass5 = rmseMass5/range_y5;
fprintf('Power tracking RMSE for 5-zone buildings (setpoint control): %f\n', rmse5);
fprintf('Power tracking RMSE for 5-zone buildings (mass flow rate control): %f\n', rmseMass5);
fprintf('Power tracking nRMSE for 5-zone buildings (setpoint control): %f\n', nrmse5);
fprintf('Power tracking nRMSE for 5-zone buildings (mass flow rate control): %f\n', nrmseMass5);


f22 = figure(22);
bldg = 1;
base1 = plot(10:0.25:18,hvac(bldg).Pbase,'-','Color',[0,0,0],"LineWidth",3,'MarkerSize',3);
hold on
EplusPagg = 0;
EplusQagg = 0;
Ptable = table2array(hvac(bldg).logTable((DRstart+1):(DRend+1),contains(hvac(bldg).logTable.Properties.VariableNames,'Electric_Power')));
EplusPagg = EplusPagg + sum(Ptable,2)/1e3;
Qtable = table2array(hvacMass(bldg).logTable((DRstart+1):(DRend+1),contains(hvacMass(bldg).logTable.Properties.VariableNames,'Electric_Power')));
EplusQagg = EplusQagg + sum(Qtable,2)/1e3;
Eplus = plot(10:0.25:18,EplusPagg,'-','Color',[0.00,0.45,0.74],"LineWidth",3,'MarkerSize',3);
p1 = plot(10:0.25:18,EplusQagg,'-','Color',[0.93,0.69,0.13],'LineWidth',3,'MarkerSize',3);
p2 = plot(10:0.25:18,PrefMass(:,bldg),'--r',"LineWidth",2,'MarkerSize',3);
hold off
xlabel('Time [h]','Interpreter','latex','fontname','Times New Roman')
ylabel('Power [kW]','Interpreter','latex','fontname','Times New Roman')
set(gca,'FontSize',12,'fontname','Times New Roman')
legend([base1,p2,Eplus,p1],["$P^{\mathrm{baseline}}_{1}$","$P^{\mathrm{ref}}_{1}$","$P^{\mathrm{setpoint}}_{1}$","$P^{\mathrm{mass}}_{1}$"],'Interpreter','latex','fontname','Times New Roman','NumColumns',2,'FontSize',14,'Location','NorthEast'); % ,
set(gcf,'position',[10,10,600,230])
xlim([10 18])
ylim([5 9])
range_y3 = max(PrefMass(:,bldg)) - min(PrefMass(:,bldg)); % 실제 값의 범위
msetemp3 = mean((EplusPagg-PrefMass(:,bldg)).^2,'all');
rmse3 = sqrt(msetemp3);
nrmse3 = rmse3/range_y3;
rmseMass3 = sqrt(mean((EplusQagg-PrefMass(:,bldg)).^2,'all'));
nrmseMass3 = rmseMass3/range_y3;
fprintf('Power tracking RMSE for 3-zone buildings (setpoint control): %f\n', rmse3);
fprintf('Power tracking RMSE for 3-zone buildings (mass flow rate control): %f\n', rmseMass3);
fprintf('Power tracking nRMSE for 3-zone buildings (setpoint control): %f\n', nrmse3);
fprintf('Power tracking nRMSE for 3-zone buildings (mass flow rate control): %f\n', nrmseMass3);

f5 = figure(5);
Qtable = table2array(hvac(bldg).logTable(DRstart:DRend,contains(hvac(bldg).logTable.Properties.VariableNames,'Sensible_Cooling_Rate')));
socTable = zeros(NumB, DRsize);
socTableMass = zeros(NumB, DRsize);
for MBidx=1:NumB
    Type2Bldg = table2array(hvac(MBidx).logTable(DRstart:DRend,contains(hvac(MBidx).logTable.Properties.VariableNames,'Zone_Air_Temperature')));
    TypeMassBldg = table2array(hvacMass(MBidx).logTable(DRstart:DRend,contains(hvacMass(MBidx).logTable.Properties.VariableNames,'Zone_Air_Temperature')));
    socTable(MBidx,:) = (hvac(MBidx).Tset'-Type2Bldg)./hvac(MBidx).delta'*(ones(1,hvac(MBidx).numZ)*inv(diag(hvac(MBidx).Btild)))'/sum(inv(diag(hvac(MBidx).Btild)),'all');
    socTableMass(MBidx,:) = (hvacMass(MBidx).Tset'-TypeMassBldg)./hvacMass(MBidx).delta'*(ones(1,hvacMass(MBidx).numZ)*inv(diag(hvacMass(MBidx).Btild)))'/sum(inv(diag(hvacMass(MBidx).Btild)),'all');
end
l1 = plot(10:0.25:18,-ones(1,length(10:0.25:18)),'r--','LineWidth',1.5);
hold on
plot(10:0.25:18,ones(1,length(10:0.25:18)),'r--','LineWidth',2);
SOCVESStemp = 1./[hvac.bhat]*(socTable)/sum(1./[hvac.bhat]);
SOCVESSmass = 1./[hvac.bhat]*(socTableMass)/sum(1./[hvac.bhat]);
a1 = plot(10:0.25:18,1./[hvac.bhat]*(socTable)/sum(1./[hvac.bhat]),'LineWidth',2,'Color',[0.00,0.45,0.74]);
a2 = plot(10:0.25:18,1./[hvac.bhat]*(socTableMass)/sum(1./[hvac.bhat]),'LineWidth',2,'Color',[0.93,0.69,0.13]);
pp = plot(10:0.25:18+0.25,soc_sched{sig},'r--','LineWidth',2,'Color',[0,0,0]);
hold off
xlim([10, 18])
ylim([-2 2])
ylabel(['SOC_{VESS}'])
xlabel(['Time [h]'])
set(gca,'FontSize',12,'fontname','Times New Roman')
legend([pp,l1,a1,a2],"$\lambda$","$\mathrm{Limits}$","$soc^\mathrm{setpoint}_\mathrm{VESS}$","$soc^\mathrm{mass}_\mathrm{VESS}$",'Interpreter','latex','fontname','Times New Roman','NumColumns',1,'FontSize',14,'Location','NorthEast'); % ,
set(gcf,'position',[10,10,600,230])
xticks(1:1:24)
SOCtempRMSE = sqrt(mean((soc_sched{sig}(1:end-1)-1./[hvac.bhat]*(socTable)/sum(1./[hvac.bhat])).^2));
SOCmassRMSE = sqrt(mean((soc_sched{sig}(1:end-1)-1./[hvac.bhat]*(socTableMass)/sum(1./[hvac.bhat])).^2));
Std_SOCtemp = sqrt(sum((socTable-SOCVESStemp).^2./[hvac.bhat]',1)/sum(1./[hvac.bhat]));
Std_SOCmass = sqrt(sum((socTableMass-SOCVESSmass).^2./[hvac.bhat]',1)/sum(1./[hvac.bhat]));
meanStd_SOCtemp = mean(Std_SOCtemp);
meanStd_SOCmass = mean(Std_SOCmass);
meanStd_SOCtempStdalpha = meanStd_SOCtemp*0.0273;
meanStd_SOCmassStdalpha = meanStd_SOCmass*0.0273;
RMSEbaseline1 = sqrt(mean((hvac(1).Pbase_true(DRstart:DRend)'-hvac(1).Pbase).^2));
RMSEbaseline11 = sqrt(mean((hvac(11).Pbase_true(DRstart:DRend)'-hvac(11).Pbase).^2));

f6 = figure(6);
pp = plot(10:0.25:18+0.25,soc_sched{sig},'k--','LineWidth',2);
xlim([10, 18])
ylim([-1.5 1.5])
ylabel(['$\mathrm{SOC}_\mathrm{VESS}$'],'Interpreter','latex','fontname','Times New Roman','FontSize',12)
xlabel(['Time [h]'])
set(gca,'FontSize',12,'fontname','Times New Roman')
legend(pp,"$\lambda$",'Interpreter','latex','fontname','Times New Roman','NumColumns',1,'FontSize',14,'Location','NorthEast'); % ,
set(gcf,'position',[10,10,600,230])
xticks(1:1:24)
ylim([-2 2])

f55 = figure(55);
plot(10:0.25:18,-ones(1,length(10:0.25:18)),'r--','LineWidth',1.5)
hold on
plot(10:0.25:18,ones(1,length(10:0.25:18)),'r--','LineWidth',1.5)
for MBidx=1:NumB
    plot(10:0.25:18,socTable(MBidx,:),'LineWidth',2)
end
hold off
xlim([10, 18])
ylim([-2 2])
ylabel(['SOC'])
xlabel(['Time [h]'])
box on;
legend('Limits','FontSize',12)
set(gca,'FontSize',12,'fontname','Times New Roman')
set(gcf,'position',[10,10,600,230])
xticks(1:1:24)
SOCvariance = mean(var(socTable, 0, 1));

f9 = figure(9);
bldg = 1;
Tistep3 = table2array(hvac(bldg).logTable(:,contains(hvac(bldg).logTable.Properties.VariableNames,'Zone_Air_Temperature')));
plot(10:0.25:18,((hvac(bldg).Tset(1)-Tistep3(DRstart:DRend,1))./hvac(bldg).delta(1)),'-','LineWidth',2,'Color',[1.00,0.41,0.16],'MarkerSize',3)
hold on
plot(10:0.25:18,((hvac(bldg).Tset(2)-Tistep3(DRstart:DRend,2))./hvac(bldg).delta(2)),'-','LineWidth',2,'Color',[0.93,0.69,0.13],'MarkerSize',3)
plot(10:0.25:18,((hvac(bldg).Tset(3)-Tistep3(DRstart:DRend,3))./hvac(bldg).delta(3)),'-','LineWidth',2,'Color',[0.80,0.50,0.10],'MarkerSize',3)
bldg = 11;
Tistep5 = table2array(hvac(bldg).logTable(:,contains(hvac(bldg).logTable.Properties.VariableNames,'Zone_Air_Temperature')));
plot(10:0.25:18,((hvac(bldg).Tset(1)-Tistep5(DRstart:DRend,1))./hvac(bldg).delta(1)),'-','LineWidth',2,'Color',[0.00,0.45,0.75],'MarkerSize',3)
plot(10:0.25:18,((hvac(bldg).Tset(2)-Tistep5(DRstart:DRend,2))./hvac(bldg).delta(2)),'-','LineWidth',2,'Color',[0.00,0.70,0.70],'MarkerSize',3)
plot(10:0.25:18,((hvac(bldg).Tset(3)-Tistep5(DRstart:DRend,3))./hvac(bldg).delta(3)),'-','LineWidth',2,'Color',[0.00,0.60,1.00],'MarkerSize',3)
plot(10:0.25:18,((hvac(bldg).Tset(4)-Tistep5(DRstart:DRend,4))./hvac(bldg).delta(4)),'-','LineWidth',2,'Color',[0.00,0.80,1.00],'MarkerSize',3)
plot(10:0.25:18,((hvac(bldg).Tset(5)-Tistep5(DRstart:DRend,5))./hvac(bldg).delta(5)),'-','LineWidth',2,'Color',[0.00,0.00,1.00],'MarkerSize',3)
plot(10:0.25:18,-ones(1,length(10:0.25:18)),'r--','LineWidth',1.5)
plot(10:0.25:18,ones(1,length(10:0.25:18)),'r--','LineWidth',1.5)
set(gca,'FontSize',12,'fontname','Times New Roman')
legend("$soz_{1,1}$","$soz_{2,1}$","$soz_{3,1}$","$soz_{1,11}$","$soz_{2,11}$","$soz_{3,11}$","$soz_{4,11}$","$soz_{5,11}$","$\mathrm{Limits}$",...
    'Interpreter','latex','fontname','Times New Roman','FontSize',12,'Orientation','horizontal','NumColumns',3,'Location','SouthWest')
xlim([10, 18])
ylim([-2 2])
ylabel(['SOZ'],'Interpreter','latex','fontname','Times New Roman')
xlabel(['Time [h]'],'Interpreter','latex','fontname','Times New Roman')
xticks(0:1:24)
yticks([-2 -1 0 1 2])

ax = gca;
ax.GridAlpha = 1;
ax.GridColor = [0.8 0.8 0.8];
set(gcf,'position',[10,10,600,230])
SOZtable = [(hvac(1).Tset' - Tistep3(DRstart:DRend,:))./hvac(1).delta',(hvac(11).Tset' - Tistep5(DRstart:DRend,:))./hvac(11).delta'];
SOZvariance = mean(var(SOZtable, 0, 2));

if ~exist('result', 'dir')
    mkdir('result');
end

% 모든 figure 저장
figHandles = findall(0, 'Type', 'figure');  % 열려있는 모든 figure 찾기
for i = 1:length(figHandles)
    fig = figHandles(i);
    filename = fullfile('result', sprintf('figure_%d.png', i));  % 저장 경로 및 파일명
    saveas(fig, [filename, '.fig'])
end
warning('on', 'all');