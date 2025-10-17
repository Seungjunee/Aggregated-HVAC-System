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

rng(1)
set(0, 'DefaultAxesFontName', 'Times');
cd Fivezone_buildings

NumB = 10;
NumZ = 5;
DRstart = 10*4;
DRend = 18*4+1;
DRsize = DRend - DRstart + 1;
building_results = struct('validTin', [], 'Tistep', [], 'mse', []); 
%% Create mlep instance and configure it
for bldg=1:NumB
    Tset = [23 + rand(5,1)];
    delta = 0.5 + 1*rand(5,1);
    ep{bldg} = mlep;
    ep{bldg}.idfFile = strcat('BuildingID',int2str(bldg));
    ep{bldg}.epwFile = 'USA_FL_Miami.Intl.AP.722020_TMY3';
    ep{bldg}.outputDirName = strcat('5zoneID',int2str(bldg));
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
    k1 = reshape(repmat(23+3*rand(1,endTime/900),900,1),[],1);
    k2 = reshape(repmat(23+3*rand(1,endTime/900),900,1),[],1);
    k3 = reshape(repmat(23+3*rand(1,endTime/900),900,1),[],1);
    k4 = reshape(repmat(23+3*rand(1,endTime/900),900,1),[],1);
    k5 = reshape(repmat(23+3*rand(1,endTime/900),900,1),[],1);
    while t < endTime
        % Prepare inputs (possibly from last outputs)
        u = [k1(t+1) k2(t+1) k3(t+1) k4(t+1) k5(t+1)];
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
    
    %% optimization for identification
    weekday_time = [];
    for week = 0:2
        weekday_time = [weekday_time;96*7*week + 96*0 + [DRstart:DRend]];
        weekday_time = [weekday_time;96*7*week + 96*1 + [DRstart:DRend]];
        weekday_time = [weekday_time;96*7*week + 96*2 + [DRstart:DRend]];
        weekday_time = [weekday_time;96*7*week + 96*3 + [DRstart:DRend]];
        weekday_time = [weekday_time;96*7*week + 96*4 + [DRstart:DRend]];
    end
    linear_time = reshape(weekday_time',1,[]);
    
    Tistep = table2array(logTable(linear_time,[9:13]));
    Coolingrate = table2array(logTable(linear_time,[8,14:17]))./1e3; % [kW]
    Tohist = table2array(logTable(linear_time,18));
    
    A = sdpvar(NumZ,NumZ);
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
    Constraints = [Constraints, A(1,3)==0, A(3,1)==0, A(2,4)==0, A(4,2)==0];
    options = sdpsettings('solver','gurobi');
    sol = optimize(Constraints,sum((Tinvar - Tistep).^2,'all'),options);
    
    valTinvar = value(Tinvar);
    valA = value(A);
    valB = value(B);
    Btild = value(B)./delta;
    valC = value(C);
    valD = value(D);
    valD = valD(:,1:end-1);
    
    figure(20)
    plot(valTinvar(:,1))
    hold on
    plot(Tistep(:,1),'--r')
    hold off
    
    nn = 0;
    validTin = zeros(DRsize,NumZ);
    validTin(1,:) = Tistep((nn*DRsize+1),:);
    for t = 1:DRsize-1
        validTin(t+1,:) = valA * validTin(t,:)' - valB.*Coolingrate(t+(nn*DRsize+1),:)' + valC * Tohist(t+(nn*DRsize+1)) + valD(:,t);
    end
    %% MSTB modeling
    a_MB = sum(inv(diag(Btild))*valA,'all')/sum(inv(diag(Btild)),'all');
    b_MB = 1/sum(inv(diag(Btild)),'all');

    figure(10)
    for i=1:5
        subplot(2,3,i)
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
    
    operate_time = 2+6*4:20*4;
    for week = [1:4, 7:11, 14:18]
        operate_time = [operate_time;96*week + [2+6*4:20*4]];
    end
    operate_time = reshape(operate_time',1,[]);
    
    Chillerpower = table2array(logTable(operate_time,19))/1e3; % [kW]
    Fanpower = table2array(logTable(operate_time,7))/1e3; % [kW]
    SumCoolingrate = sum(logTable{operate_time,[8,14:17]},2)/1e3; % [kW]
    Totalpower = Chillerpower+Fanpower;
    
    fig4 = figure(4);
    mdl = fitlm(SumCoolingrate,Chillerpower+Fanpower); % ,'intercept', false
    coeff_1 = mdl.Coefficients{2,1};
    coeff_2 = mdl.Coefficients{1,1};
    plot(SumCoolingrate,Chillerpower+Fanpower,'b.','MarkerSize',8)
    hold on
    plot(linspace(0,15,100),linspace(0,15,100)*coeff_1+coeff_2,'k') %
    hold off
    xlabel('Q(k) [kW]')
    ylabel('P(k) [kW]')
    xlim([0,15])
    ylim([0,10])
    set(gca,'FontSize',26)
    set(gcf,'position',[10,10,600,600])
    savefig(fig4,"5-zone P-Q")
    
    mtot = sum(logTable{operate_time,[2:6]},2);
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

    figure(2);
    plot(20:0.1:26, 20:0.1:26,'k','LineWidth',2)
    hold on
    plot(Tistep((nn*DRsize+1):(nn+1)*DRsize,:),validTin,'.','MarkerSize',12);
    save('5zonetemp.mat','Tistep','validTin')
    hold off
    xlim([23 26])
    ylim([23 26])
    ylabel(['Estimated Temp. [' char(176) 'C]'])
    xlabel(['Actual Temp. [' char(176) 'C]'])
    set(gca,'FontSize',22)
    set(gcf,'position',[10,10,600,600])
    
    mse = mean((Tistep((nn*DRsize+1):(nn+1)*DRsize,:)-validTin).^2,'all');
    rmse = sqrt(mse);                % Root Mean Squared Error
    
    cair = 1.03;
    folderName = strcat('5zoneID',int2str(bldg));
    fileName   = strcat('BuildingID',int2str(bldg),'.eio');
    fullPath   = fullfile(folderName, fileName);
    fid = fopen(fullPath, 'r');
    if fid == -1
        error('Error1: %s', fullPath);
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
                fprintf('Error2: %s\n', lastValueStr);
            end
        end
        tline = fgetl(fid);
    end
    m_low = m_high*0;
    qmax = cair * m_high * (24-13) *1.2;
    qmin = cair * m_low * (24-13) * 1.2;
    
    building_results(bldg).validTin = validTin;
    building_results(bldg).Tistep = Tistep((nn*DRsize+1):(nn+1)*DRsize,:);
    building_results(bldg).mse = mse;
    
    delete(strcat('5zoneID',int2str(bldg),'/*'));
end
nn=0;
fig22 = figure(22);
hold on
all_mse = [];
for bldg = 1:NumB
    plot(building_results(bldg).Tistep, building_results(bldg).validTin, '.k', 'MarkerSize', 16);
    all_mse = [all_mse; (building_results(bldg).Tistep - building_results(bldg).validTin).^2];
end

% Total RMSE Calculation
mse = mean(all_mse, 'all');
rmse = sqrt(mse);
fprintf('Total RMSE for all buildings: %f\n', rmse);
plot(20:0.1:26, 20:0.1:26,'r','LineWidth',2)
hold off
xlim([23 26])
ylim([23 26])
yticks([23:1:26])
xlabel(['Actual Temperature [' char(176) 'C]'])
ylabel(['Predicted Temperature [' char(176) 'C]'])
set(gca,'FontSize',24)
set(gcf,'position',[10,10,600,600])
box on;

cd ../
cd result
savefig(fig22,'5-zone Identification');
cd ../
