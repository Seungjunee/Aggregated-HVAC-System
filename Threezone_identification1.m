%% Simple Matlab <-> EnergyPlus co-simulation example
% Demonstrates the functionality of the mlep (MatLab-EnergyPlus) tool in
% a small office building simulation scenario.
%
% Note that a start of the simulation period as well as a timestep and
% an input/output configuration is defined by the the EnergyPlus simulation
% configuration file (.IDF). Climatic conditions are obtained from a
% EnergyPlus Weather data file (.EPW).

clc
clear

rng(2)
set(0, 'DefaultAxesFontName', 'Times');

cd Threezone_buildings
NumB = 10;
NumZ = 3;
DRstart = 10*4;
DRend = 18*4 + 1;
DRsize = DRend - DRstart + 1;
building_results = struct('validTin', [], 'Tistep', [], 'mse', []); 
%% Create mlep instance and configure it
for bldg = 1:NumB
    Tset = 25 * ones(NumZ,1) + rand(NumZ,1);
    delta = 1 + 0.5*rand(NumZ,1);
%     Tset = 25* ones(NumZ,1);
%     delta = 1* ones(NumZ,1);
    ep{bldg} = mlep;
    ep{bldg}.idfFile = strcat('BuildingID',int2str(bldg));
    ep{bldg}.epwFile = 'USA_FL_Miami.Intl.AP.722020_TMY3';
    ep{bldg}.outputDirName = strcat('3zoneID',int2str(bldg));
    ep{bldg}.initialize;
    
    %% Simulate
    % Specify simulation duration
    endTime = 19*24*60*60; %[s]
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
    k1 = reshape(repmat(24+3*rand(1,endTime/900),900,1),[],1);
    k2 = reshape(repmat(24+3*rand(1,endTime/900),900,1),[],1);
    k3 = reshape(repmat(24+3*rand(1,endTime/900),900,1),[],1);
    while t < endTime
        % Prepare inputs (possibly from last outputs)
        u = [k1(t+1) k2(t+1) k3(t+1)];
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
    
    % logTable_15min = logTable(1:15:end,:);
    % for i = 2:38
    %     temp = reshape(logTable{2:end,i},15,[]);
    %     logTable_15min{1:end-1,i} = mean(temp)';
    % end
    %% optimization for identification
    weekday_time = [];
    for week = 0:2
        weekday_time = [weekday_time;96*7*week + 96*0 + [DRstart:DRend]];
        weekday_time = [weekday_time;96*7*week + 96*1 + [DRstart:DRend]];
        weekday_time = [weekday_time;96*7*week + 96*2 + [DRstart:DRend]];
        weekday_time = [weekday_time;96*7*week + 96*3 + [DRstart:DRend]];
        weekday_time = [weekday_time;96*7*week + 96*4 + [DRstart:DRend]];
    end
%     weekday_time = 1:96*5;
%     for week = 1:2
%         weekday_time = [weekday_time;96*7*week + [4*6:96*5]];
%     end
    linear_time = reshape(weekday_time',1,[]);
    Tistep = table2array(logTable(linear_time,contains(logTable.Properties.VariableNames,'Zone_Air_Temperature')));
    Coolingrate = table2array(logTable(linear_time,contains(logTable.Properties.VariableNames,'Zone_Air_System_Sensible_Cooling_Rate')))./1e3; % [kW]
    Tohist = table2array(logTable(linear_time,contains(logTable.Properties.VariableNames,'Drybulb_Temperature')));
    
    A = sdpvar(NumZ,NumZ,'full');
    B = sdpvar(NumZ,1);
    C = sdpvar(NumZ,1);
    D = sdpvar(NumZ,DRsize);
    Tinvar = sdpvar(length(linear_time),NumZ);
    Constraints = [];
    for days = 1:size(weekday_time,1)
        Constraints = [Constraints, Tinvar(((days-1)*DRsize+1),:) == Tistep(((days-1)*DRsize+1),:)];
        for t = ((days-1)*DRsize+1):days*DRsize
            if rem(t,DRsize) ~= 0
                Constraints = [Constraints, Tinvar(t+1,:)' == A * Tistep(t,:)' - B.*Coolingrate(t+1,:)' + C * Tohist(t+1) + D(:,rem(t,DRsize))];
            end
        end
    end
%     Conn = ones(NumZ,NumZ);
%     for i=1:NumZ
%         for j=1:NumZ
%             if Conn(i,j) ~= 1
%                 Constraints = [Constraints, A(i,j) == 0];
%             end
%         end
%     end
    options = sdpsettings('solver','gurobi');
    sol = optimize(Constraints,sum((Tinvar - Tistep).^2,'all'),options);
    
    valTinvar = value(Tinvar);
    valA = value(A);
    valB = value(B);
    Btild = value(B)./delta;
    valC = value(C);
    valD = value(D);
    valD = valD(:,1:end-1);
    
    nn = 0;
    validTin = zeros(DRsize,NumZ);
    validTin(1,:) = Tistep((nn*DRsize+1),:);
    for t = 1:DRsize-1
        validTin(t+1,:) = valA * validTin(t,:)' - valB.*Coolingrate(t+(nn*DRsize+1),:)' + valC * Tohist(t+(nn*DRsize+1)) + valD(:,t);
    end
    %% MSTB modeling
    a_MB = sum(inv(diag(Btild))*valA,'all')/sum(inv(diag(Btild)),'all');
    b_MB = 1/sum(inv(diag(Btild)),'all');
    
    fig1 = figure(1);
    for i=1:NumZ
        subplot(3,3,i)
        plot(DRstart*0.25:0.25:DRend*0.25,validTin(:,i))
        hold on
        plot(DRstart*0.25:0.25:DRend*0.25,Tistep((nn*DRsize+1):(nn+1)*DRsize,i))
        hold off
        xticks([0:4:24])
        xlim([0, 24])
        ylabel(['Indoor temperature [' char(176) 'C]'])
        xlabel(['Time [h]'])
        set(gca,'FontSize',15)
    end
    fig2 = figure(2);
    plot(0:0.1:100, 0:0.1:100,'k','LineWidth',2)
    hold on
    plot(Tistep((nn*DRsize+1):(nn+1)*DRsize,:),validTin,'.','MarkerSize',12);
    hold off
    xlim([24 27])
    ylim([24 27])
    ylabel(['Estimated Temp. [' char(176) 'C]'])
    xlabel(['Actual Temp. [' char(176) 'C]'])
    set(gca,'FontSize',22)
    set(gcf,'position',[10,10,600,600])
    
    fig3 = figure(3);
    plot(10:0.25:18.25,Tistep((nn*DRsize+1):(nn+1)*DRsize,1),'MarkerSize',12);
    hold on
    plot(10:0.25:18.25,validTin(:,1),'MarkerSize',12);
    hold off
    
    mse = mean((Tistep((nn*DRsize+1):(nn+1)*DRsize,:)-validTin).^2,'all');
    rmse = sqrt(mse);                % Root Mean Squared Error
    operate_time = 2+10*4:17*4;
    for week = [1:4, 7:11, 14:18]
        operate_time = [operate_time;96*week + [2+10*4:17*4]];
    end
    operate_time = reshape(operate_time',1,[]);
    
    Chillerpower = sum(table2array(logTable(operate_time,contains(logTable.Properties.VariableNames,'Chiller_Electric_Power'))),2)/1e3; % [kW]
    Fanpower = sum(table2array(logTable(operate_time,contains(logTable.Properties.VariableNames,'Fan_Electric_Power'))),2)/1e3; % [kW]
    SumCoolingrate = sum(logTable{operate_time,contains(logTable.Properties.VariableNames,'Zone_Air_System_Sensible_Cooling_Rate')},2)/1e3; % [kW]
    Totalpower = Chillerpower+Fanpower;
    
    figure(4)
    mdl = fitlm(SumCoolingrate,Chillerpower+Fanpower); % ,'intercept', false
    coeff_1 = mdl.Coefficients{2,1};
    coeff_2 = mdl.Coefficients{1,1};
    plot(mdl)
    hold on
    plot(linspace(0,14,100),linspace(0,14,100)*coeff_1+coeff_2,'k')
    hold off
    xlabel('Q(k) [kW]')
    ylabel('P(k) [kW]')
    set(gca,'FontSize',15)
    
    mtot = sum(logTable{operate_time,contains(logTable.Properties.VariableNames,'System_Node_Mass_Flow_Rate')},2);
    
    figure(5)
    mdl2 = polyfit(mtot,Fanpower,2); % ,'intercept', false
    Fancoeff_1 = mdl2(1);
    Fancoeff_2 = mdl2(2);
    Fancoeff_3 = mdl2(3);
    scatter(mtot,Fanpower)
    hold on
    plot(linspace(0,1.3,100),mdl2(1)*linspace(0,1.3,100).^2 + mdl2(2)*linspace(0,1.3,100) + mdl2(3))
    hold off
    xlabel('Q(k) [kW]')
    ylabel('P(k) [kW]')
    set(gca,'FontSize',15)
        
    folderName = strcat('3zoneID',int2str(bldg));         % 하위 폴더 이름
    fileName   = strcat('BuildingID',int2str(bldg),'.eio');     % 대상 텍스트 파일 이름
    fullPath   = fullfile(folderName, fileName);
    fid = fopen(fullPath, 'r');
    if fid == -1
        error('파일을 열 수 없습니다: %s', fullPath);
    end
    m_high = [];
    tline = fgetl(fid);
    while ischar(tline)
        if contains(tline, 'User-Specified Maximum Air Flow Rate [m3/s]')
            splitLine = strsplit(tline, ',');
            lastValueStr = splitLine{end};
            lastValueNum = str2double(lastValueStr);
            if ~isnan(lastValueNum)
                m_high(end+1) = lastValueNum; %#ok<SAGROW>
            else
                fprintf('문자열을 숫자로 변환할 수 없습니다: %s\n', lastValueStr);
            end
        end
        tline = fgetl(fid);
    end
    m_low = m_high*0;
    qmax = m_high * (25-13)*1.2;
    qmin = m_low * (25-13)*1.2;
    
    building_results(bldg).validTin = validTin;
    building_results(bldg).Tistep = Tistep((nn*DRsize+1):(nn+1)*DRsize,:);
    building_results(bldg).mse = mse;
    
    delete(strcat('3zoneID',int2str(bldg),'/*'));
%     save(strcat('coefficients',int2str(bldg),'.mat'),'Tset','delta','m_high','m_low','valA','valB','valC','valD','coeff_1','coeff_2','Fancoeff_1','Fancoeff_2','Fancoeff_3',...
%         'a_MB','b_MB','validTin','Tistep','DRsize','qmin','qmax')
end
% nn=0;
% NumB = 10;
% NumZ = 3;
% DRstart = 10*4;
% DRend = 18*4 + 1;
% DRsize = DRend - DRstart + 1;
fig22 = figure(22);
% subplot(1,2,1);
hold on
% validTin = zeros(DRsize,NumZ);
% validTin(1,:) = Tistep((nn*DRsize+1),:);
all_mse = [];
for bldg = 1:NumB
    plot(building_results(bldg).Tistep, building_results(bldg).validTin, '.k', 'MarkerSize', 16);
    all_mse = [all_mse; (building_results(bldg).Tistep - building_results(bldg).validTin).^2];
end

% 전체 RMSE 계산
mse = mean(all_mse, 'all');
rmse = sqrt(mse);
fprintf('Total RMSE for all buildings: %f\n', rmse);

plot(20:0.1:30, 20:0.1:30,'r','LineWidth',2)
% fivezone
% subplot(1,2,2);
% cd ../Fivezone_buildings
% plot(Tistep((nn*DRsize+1):(nn+1)*DRsize,:),validTin,'.k','MarkerSize',12);
% for bldg = 2:10
%     load(strcat('coefficients',int2str(bldg),'.mat'));
%     plot(Tistep((nn*DRsize+1):(nn+1)*DRsize,:),validTin,'.k','MarkerSize',12);
% end
% plot(20:0.1:26, 20:0.1:26,'r','LineWidth',2)
hold off
xlim([24 27])
ylim([24 27])
yticks([24:1:27])
xlabel(['Actual Temperature [' char(176) 'C]'])
ylabel(['Predicted Temperature [' char(176) 'C]'])

set(gca,'FontSize',24)
set(gcf,'position',[10,10,600,600])
box on;
cd ../