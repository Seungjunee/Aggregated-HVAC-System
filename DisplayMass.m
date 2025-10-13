load('CaseMass1.mat'); sig = 1;
NumB = 20;
DRstart = 24*4*9 + 10*4;
DRend = 24*4*9 + 18*4;
DRsize = DRend - DRstart + 1;

Tistep5 = [hvacMass(11).logTable.EP_SPACE1_1__Zone_Air_Temperature, hvacMass(11).logTable.EP_SPACE2_1__Zone_Air_Temperature, ...
    hvacMass(11).logTable.EP_SPACE3_1__Zone_Air_Temperature, hvacMass(11).logTable.EP_SPACE4_1__Zone_Air_Temperature, ...
    hvacMass(11).logTable.EP_SPACE5_1__Zone_Air_Temperature];
Tistep3 = table2array(hvacMass(1).logTable(:,contains(hvacMass(1).logTable.Properties.VariableNames,'Zone_Air_Temperature')));
Paggbase = 0;
for bldg = 1:NumB
    Paggbase = Paggbase + hvacMass(bldg).Pbase';
end
PaggbaseTrue = 0;
for bldg = 1:NumB
    PaggbaseTrue = PaggbaseTrue + hvacMass(bldg).Pbase_true(DRstart:DRend);
end

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
% xlim([10 18])
ylim([20 30])
xticks([0:1:24])
xlabel('Time [h]','Interpreter','latex','fontname','Times New Roman')
ylabel(['Temperature [' char(176) 'C]'])
set(gcf,'position',[10,10,600,230])

f1 = figure(1);
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
set(gca,'FontSize',13,'fontname','Times New Roman')
legend([base1,p2,Eplus,p1],["Baseline","$P_{ref}$","Actual (Eplus)","DR signal"],'Interpreter','latex','fontname','Times New Roman','NumColumns',2,'FontSize',13); % ,
set(gcf,'position',[10,10,600,230])
mse = mean((Paggbase+P_sched{sig}'-EplusPagg).^2);
rmse = sqrt(mse);                % Root Mean Squared Error
range_y = max(Paggbase+P_sched{sig}') - min(Paggbase+P_sched{sig}'); % 실제 값의 범위
nrmse = rmse / range_y; % NRMSE 계산

f11 = figure(11);
plot(10:0.25:18,P_sched{sig}','-o','Color',[0.8,0.2,1],'LineWidth',3,'MarkerSize',3);
hold on
plot(10:0.25:18,-Pbaggmin,'Color',[0.6,0.6,0.6],'LineWidth',0.5);
plot(10:0.25:18,Pbaggmax,'Color',[0.6,0.6,0.6],'LineWidth',0.5);
% plot(DRstart*0.25:0.25:DRend*0.25,-PbaggminRelax,'Color',[1,0,0],'LineWidth',0.5);
% plot(DRstart*0.25:0.25:DRend*0.25,PbaggmaxRelax,'Color',[1,0,0],'LineWidth',0.5);
hold off
xlabel('Time [h]','Interpreter','latex','fontname','Times New Roman')
ylabel('Power [kW]','Interpreter','latex','fontname','Times New Roman')
set(gca,'FontSize',12,'fontname','Times New Roman')
set(gcf,'position',[10,10,600,230])
xlim([10 18])
% ylim([-230 1000])

f2 = figure(2);
bldg = 11;
% base1 = plot(10:0.25:18,hvacMass(bldg).Pbase_true(DRstart:DRend),'-o','Color',[0,0,0],"LineWidth",3,'MarkerSize',3);
base1 = plot(10:0.25:18,hvacMass(bldg).Pbase,'-o','Color',[0,0,0],"LineWidth",3,'MarkerSize',3);
hold on
EplusPagg = 0;
EplusQagg = 0;
Ptable = table2array(hvacMass(bldg).logTable((DRstart+1):(DRend+1),contains(hvacMass(bldg).logTable.Properties.VariableNames,'Electric_Power')));
EplusPagg = EplusPagg + sum(Ptable,2)/1e3;
Qtable = table2array(hvacMass(bldg).logTable((DRstart+1):(DRend+1),contains(hvacMass(bldg).logTable.Properties.VariableNames,'Sensible_Cooling_Rate')));
EplusQagg = EplusQagg + sum(Qtable,2)/1e3;
Eplus = plot(10:0.25:18,EplusPagg,'-o','Color',[0.00,0.45,0.74],"LineWidth",3,'MarkerSize',3);
p2 = plot(10:0.25:18,Pref(:,bldg),'--ro',"LineWidth",2,'MarkerSize',3);
% p1 = plot(10:0.25:18,P_sched,'-o','Color',[0.8,0.2,1],'LineWidth',3,'MarkerSize',3);
plot([10,10],[-2,10],'Color',[0.6,0.6,0.6],'LineWidth',0.5);
plot([18,18],[-2,10],'Color',[0.6,0.6,0.6],'LineWidth',0.5);
hold off
xlabel('Time [h]','Interpreter','latex','fontname','Times New Roman')
ylabel('Power [kW]','Interpreter','latex','fontname','Times New Roman')
set(gca,'FontSize',12,'fontname','Times New Roman')
legend([base1,p2,Eplus],["$P^{\mathrm{baseline}}_{11}$","$P^{\mathrm{ref}}_{11}$","$P^{\mathrm{actual}}_{11}$"],'Interpreter','latex','fontname','Times New Roman','NumColumns',2,'FontSize',13,'Location','SouthEast'); % ,
set(gcf,'position',[10,10,600,230])
xlim([10 18])
ylim([4 9])
mse = mean((EplusPagg-Pref(:,bldg)).^2,'all');
rmse = sqrt(mse);

f22 = figure(22);
bldg = 1;
% base1 = plot(10:0.25:18,hvacMass(bldg).Pbase_true(DRstart:DRend),'-o','Color',[0,0,0],"LineWidth",3,'MarkerSize',3);
base1 = plot(10:0.25:18,hvacMass(bldg).Pbase,'-o','Color',[0,0,0],"LineWidth",3,'MarkerSize',3);
hold on
EplusPagg = 0;
EplusQagg = 0;
Ptable = table2array(hvacMass(bldg).logTable((DRstart+1):(DRend+1),contains(hvacMass(bldg).logTable.Properties.VariableNames,'Electric_Power')));
EplusPagg = EplusPagg + sum(Ptable,2)/1e3;
Qtable = table2array(hvacMass(bldg).logTable((DRstart+1):(DRend+1),contains(hvacMass(bldg).logTable.Properties.VariableNames,'Sensible_Cooling_Rate')));
EplusQagg = EplusQagg + sum(Qtable,2)/1e3;
Eplus = plot(10:0.25:18,EplusPagg,'-o','Color',[0.00,0.45,0.74],"LineWidth",3,'MarkerSize',3);
p2 = plot(10:0.25:18,PrefMass(:,bldg),'--ro',"LineWidth",2,'MarkerSize',3);
% p1 = plot(10:0.25:18,P_sched,'-o','Color',[0.8,0.2,1],'LineWidth',3,'MarkerSize',3);
plot([10,10],[-2,10],'Color',[0.6,0.6,0.6],'LineWidth',0.5);
plot([18,18],[-2,10],'Color',[0.6,0.6,0.6],'LineWidth',0.5);
hold off
xlabel('Time [h]','Interpreter','latex','fontname','Times New Roman')
ylabel('Power [kW]','Interpreter','latex','fontname','Times New Roman')
set(gca,'FontSize',13,'fontname','Times New Roman')
legend([base1,p2,Eplus],["$P^{\mathrm{baseline}}_{1}$","$P^{\mathrm{ref}}_{1}$","$P^{\mathrm{actual}}_{1}$"],'Interpreter','latex','fontname','Times New Roman','NumColumns',2,'FontSize',13,'Location','SouthEast'); % ,
set(gcf,'position',[10,10,600,230])
xlim([10 18])
ylim([5 7])
mse = mean((EplusPagg-PrefMass(:,bldg)).^2,'all');
rmse = sqrt(mse);

% f6 = figure(6);
% bldg = 1;
% base1 = plot(0.25:0.25:24,hvacMass(bldg).Pbase_true,'-o','Color',[0,0,0],"LineWidth",3,'MarkerSize',3);
% hold on
% base2 = plot(10:0.25:18,hvacMass(bldg).Pbase,'-o','Color',[0.5,0.5,0.5],"LineWidth",3,'MarkerSize',3);
% EplusPagg = 0;
% EplusQagg = 0;
% 
% Ptable = table2array(hvacMass(bldg).logTable(:,contains(hvacMass(bldg).logTable.Properties.VariableNames,'Electric_Power')));
% EplusPagg = EplusPagg + sum(Ptable,2)/1e3;
% Qtable = table2array(hvacMass(bldg).logTable(:,contains(hvacMass(bldg).logTable.Properties.VariableNames,'Sensible_Cooling_Rate')));
% EplusQagg = EplusQagg + sum(Qtable,2)/1e3;
% 
% Eplus = plot(0:0.25:23.75,EplusPagg,'-o','Color',[0.00,0.45,0.74],"LineWidth",3,'MarkerSize',3);
% p2 = plot(DRstart*0.25:0.25:DRend*0.25,PrefMass(:,bldg),'--ro',"LineWidth",2,'MarkerSize',3);
% % p1 = plot(10:0.25:18,P_sched,'-o','Color',[0.8,0.2,1],'LineWidth',3,'MarkerSize',3);
% plot([10,10],[-2,10],'Color',[0.6,0.6,0.6],'LineWidth',0.5);
% plot([18,18],[-2,10],'Color',[0.6,0.6,0.6],'LineWidth',0.5);
% hold off
% xlabel('Time [h]','Interpreter','latex','fontname','Times New Roman')
% ylabel('Power [kW]','Interpreter','latex','fontname','Times New Roman')
% set(gca,'FontSize',13,'fontname','Times New Roman')
% legend([base1,p2,Eplus,base2],["Baseline","$P_{ref}$","Actual (Eplus)"],'Interpreter','latex','fontname','Times New Roman','NumColumns',2,'FontSize',12,'Location','SouthEast'); % ,
% set(gcf,'position',[10,10,600,230])
% xlim([10 18])
% % ylim([4 8])


f5 = figure(5);
Qtable = table2array(hvacMass(bldg).logTable(DRstart:DRend,contains(hvacMass(bldg).logTable.Properties.VariableNames,'Sensible_Cooling_Rate')));
socTable = zeros(NumB, DRsize);
for MBidx=1:NumB
    Type2Bldg = table2array(hvacMass(MBidx).logTable(DRstart:DRend,contains(hvacMass(MBidx).logTable.Properties.VariableNames,'Zone_Air_Temperature')));
    socTable(MBidx,:) = (hvacMass(MBidx).Tset'-Type2Bldg)./hvacMass(MBidx).delta'*(ones(1,hvacMass(MBidx).numZ)*inv(diag(hvacMass(MBidx).Btild)))'/sum(inv(diag(hvacMass(MBidx).Btild)),'all');
end
plot(10:0.25:18,mean(socTable),'LineWidth',1.5)
hold on
pp = plot(10:0.25:18+0.25,soc_sched{sig},'r--','LineWidth',2);
hold off
xlim([10, 18])
ylim([-3 3])
ylabel(['SOC'])
xlabel(['Time [h]'])
set(gca,'FontSize',13,'fontname','Times New Roman')
legend(pp,"$\lambda$",'Interpreter','latex','fontname','Times New Roman','NumColumns',1,'FontSize',13,'Location','SouthEast'); % ,
set(gcf,'position',[10,10,600,230])
xticks(1:1:24)

f55 = figure(55);
SOCtable = [];
plot(10:0.25:18,-ones(1,length(10:0.25:18)),'r--','LineWidth',1.5)
hold on
plot(10:0.25:18,ones(1,length(10:0.25:18)),'r--','LineWidth',1.5)
for MBidx=1:NumB
    plot(10:0.25:18,socTable(MBidx,:),'LineWidth',2)
end
hold off
xlim([10, 18])
ylim([-3 3])
ylabel(['SOC'])
xlabel(['Time [h]'])
box on;
set(gca,'FontSize',12,'fontname','Times New Roman')
set(gcf,'position',[10,10,600,230])
% xticks(1:1:24)
ylim([-2 2])
% legend(string(1:NumB),'NumColumns',4)
legend('Limits','FontSize',12)
SOCvariance = mean(var(socTable, 0, 1));

f8 = figure(8);
plot(10:0.25:18,((hvacMass(11).Tset'-Tistep5(DRstart:DRend,:))./hvacMass(11).delta'),'-o','LineWidth',3,'MarkerSize',3)
% legend("SOC (T_{real})","SOC (T_{est})","SOC (27)", "SOC (35)")
legend("$soz_{1,11}$","$soz_{2,11}$","$soz_{3,11}$","$soz_{4,11}$","$soz_{5,11}$",'Interpreter','latex','fontname','Times New Roman','FontSize',14,'Orientation','horizontal','NumColumns',2,'Location','SouthWest')
% xticks([8:4:20])
xlim([10, 18])
ylim([-2 2])
ylabel(['SOZ'],'Interpreter','latex','fontname','Times New Roman')
xlabel(['Time [h]'],'Interpreter','latex','fontname','Times New Roman')
xticks(0:1:24)
yticks([-2 -1 0 1 2])
set(gca,'FontSize',13,'fontname','Times New Roman')
ax = gca;
ax.GridAlpha = 1;
ax.GridColor = [0.8 0.8 0.8];
set(gcf,'position',[10,10,600,230])

f9 = figure(9);
bldg = 1;
Tistep3 = table2array(hvacMass(bldg).logTable(:,contains(hvacMass(bldg).logTable.Properties.VariableNames,'Zone_Air_Temperature')));
plot(10:0.25:18,((hvacMass(bldg).Tset(1)-Tistep3(DRstart:DRend,1))./hvacMass(bldg).delta(1)),'-','LineWidth',2,'Color',[1.00,0.41,0.16],'MarkerSize',3)
hold on
plot(10:0.25:18,((hvacMass(bldg).Tset(2)-Tistep3(DRstart:DRend,2))./hvacMass(bldg).delta(2)),'-','LineWidth',2,'Color',[0.93,0.69,0.13],'MarkerSize',3)
plot(10:0.25:18,((hvacMass(bldg).Tset(3)-Tistep3(DRstart:DRend,3))./hvacMass(bldg).delta(3)),'-','LineWidth',2,'Color',[0.80,0.50,0.10],'MarkerSize',3)
bldg = 11;
Tistep5 = table2array(hvacMass(bldg).logTable(:,contains(hvacMass(bldg).logTable.Properties.VariableNames,'Zone_Air_Temperature')));
plot(10:0.25:18,((hvacMass(bldg).Tset(1)-Tistep5(DRstart:DRend,1))./hvacMass(bldg).delta(1)),'-','LineWidth',2,'Color',[0.00,0.45,0.75],'MarkerSize',3)
plot(10:0.25:18,((hvacMass(bldg).Tset(2)-Tistep5(DRstart:DRend,2))./hvacMass(bldg).delta(2)),'-','LineWidth',2,'Color',[0.00,0.70,0.70],'MarkerSize',3)
plot(10:0.25:18,((hvacMass(bldg).Tset(3)-Tistep5(DRstart:DRend,3))./hvacMass(bldg).delta(3)),'-','LineWidth',2,'Color',[0.00,0.60,1.00],'MarkerSize',3)
plot(10:0.25:18,((hvacMass(bldg).Tset(4)-Tistep5(DRstart:DRend,4))./hvacMass(bldg).delta(4)),'-','LineWidth',2,'Color',[0.00,0.80,1.00],'MarkerSize',3)
plot(10:0.25:18,((hvacMass(bldg).Tset(5)-Tistep5(DRstart:DRend,5))./hvacMass(bldg).delta(5)),'-','LineWidth',2,'Color',[0.00,0.00,1.00],'MarkerSize',3)
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
SOZtable = [(hvac(1).Tset' - Tistep3(DRstart:DRend,:))./hvacMass(1).delta',(hvac(11).Tset' - Tistep5(DRstart:DRend,:))./hvacMass(11).delta'];
SOZvariance = mean(var(SOZtable, 0, 2));