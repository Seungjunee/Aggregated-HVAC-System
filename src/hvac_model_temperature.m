classdef hvac_model_temperature
    properties
        ep
        ts
        Tset
        numZ
        delta
        coeff_1
        coeff_2
        Fancoeff_1
        Fancoeff_2
        Fancoeff_3
        A
        B
        Btild
        C
        D
        qbase
        qbase_true
        Qbase
        Pbase
        Pbase_true
        Pbmin
        Pbmax
        a
        a_MB
        bhat
        b_MB
        Ti
        soz
        soc
        Power
        logTable
        Massflowbase
        m_low
        m_high
    end
    
    methods
        function obj = hvac_model_temperature(ts,buildidx,folderpath,controlType)
            %   initialization
            %   ts : time step
            %   Ti0 : initial indoor temperature
            previousFolder = pwd;
            folderCleanup = onCleanup(@() cd(previousFolder));
            cd(folderpath)
            obj.ep = mlep;
            obj.ep.idfFile = strcat('Building',controlType,int2str(buildidx));
            obj.ep.epwFile = 'USA_FL_Miami.Intl.AP.722020_TMY3';
            obj.ep.outputDirName = strcat('BldgTempControl',int2str(buildidx));
            % Prepare data logging            
            coefficientData = load(strcat('coefficients',int2str(buildidx),'.mat'));
            baselineData = load(strcat('Baseline_info',int2str(buildidx),'.mat'));
            baselineFields = fieldnames(baselineData);
            for fieldIndex = 1:numel(baselineFields)
                fieldName = baselineFields{fieldIndex};
                coefficientData.(fieldName) = baselineData.(fieldName);
            end
            
            obj.ts = ts;    % time step
            obj.Tset = coefficientData.Tset;%   degC
            obj.numZ = length(obj.Tset);
            obj.delta = coefficientData.delta;%   degC
            
            obj.coeff_1 = coefficientData.coeff_1;
            obj.coeff_2 = coefficientData.coeff_2;
            obj.Fancoeff_1 = coefficientData.Fancoeff_1;
            obj.Fancoeff_2 = coefficientData.Fancoeff_2;
            obj.Fancoeff_3 = coefficientData.Fancoeff_3;
            
            obj.A = coefficientData.valA;
            obj.B = coefficientData.valB;
            obj.Btild = coefficientData.valB./obj.delta;
            obj.C = coefficientData.valC;
            obj.D = coefficientData.valD;
            
            obj.qbase = coefficientData.qbase;
            obj.qbase_true = coefficientData.qbase_true;
            obj.Qbase = coefficientData.Qbase;
            obj.Pbase = coefficientData.Pbase;
            obj.Pbase_true = coefficientData.Pbase_true;
            
            obj.Pbmin = coefficientData.Pbmin;
            obj.Pbmax = coefficientData.Pbmax;
            
            obj.a = coefficientData.a_MB; %   dynamics param a
            obj.a_MB = coefficientData.a_MB;
            obj.b_MB = coefficientData.b_MB;
            obj.bhat = obj.b_MB/obj.coeff_1;%   dynamics param b
            obj.ep.initialize
            obj.logTable = table('Size',[0, 1 + obj.ep.nOut],...
                'VariableTypes',repmat({'double'},1,1 + obj.ep.nOut),...
                'VariableNames',[{'Time'}; obj.ep.outputSigName]);
            obj.Massflowbase = coefficientData.Massflowbase;
            obj.m_low = coefficientData.m_low;
            obj.m_high = coefficientData.m_high;
            clear folderCleanup;
        end
        
        function obj = stop(obj)
            obj.ep.stop;
        end
        
        function soz = eval_soz(obj,Ti)
            soz = (obj.Tset-Ti)./obj.delta;
        end
        
        function Tindoor = eval_Tindoor(obj,soc)
            Tindoor = obj.Tset - obj.delta.*soc;
        end
        
        function Qbase = eval_Qbase(obj,To,Qp)
            Qbase = (To - obj.Tset)/obj.Rth + Qp;
        end
        
        function Pbase = eval_Pbase(obj,To,Qp)
            Pbase = (To - obj.Tset)*obj.k1/(obj.k2*obj.Rth) + ...
                (obj.k1*Qp + obj.k2*obj.l1 - obj.k1*obj.l2)/obj.k2;
        end
        
        function obj = eval_Pbconst(obj,Pbase)
            obj.Pbmin = Pbase - obj.Pmin;
            obj.Pbmax = obj.Pmax - Pbase;
        end
        
    end
    
end

