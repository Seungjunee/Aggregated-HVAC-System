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
            cd(folderpath)
            obj.ep = mlep;
            obj.ep.idfFile = strcat('Building',controlType,int2str(buildidx));
            obj.ep.epwFile = 'USA_FL_Miami.Intl.AP.722020_TMY3';
            obj.ep.outputDirName = strcat('BldgTempControl',int2str(buildidx));
            % Prepare data logging            
            load(strcat('coefficients',int2str(buildidx),'.mat'))
            load(strcat('Baseline_info',int2str(buildidx),'.mat'))
            
            obj.ts = ts;    % time step
            obj.Tset = Tset;%   degC
            obj.numZ = length(Tset); 
            obj.delta = delta;%   degC
            
            obj.coeff_1 = coeff_1;
            obj.coeff_2 = coeff_2;
            obj.Fancoeff_1 = Fancoeff_1;
            obj.Fancoeff_2 = Fancoeff_2;
            obj.Fancoeff_3 = Fancoeff_3;
            
            obj.A = valA;
            obj.B = valB;
            obj.Btild = valB./obj.delta;
            obj.C = valC;
            obj.D = valD;
            
            obj.qbase = qbase;
            obj.qbase_true = qbase_true;
            obj.Qbase = Qbase;
            obj.Pbase = Pbase;
            obj.Pbase_true = Pbase_true;
            
            obj.Pbmin = Pbmin;
            obj.Pbmax = Pbmax;
            
            obj.a = a_MB; %   dynamics param a
            obj.b_MB = b_MB;
            obj.bhat = b_MB/obj.coeff_1;%   dynamics param b
            obj.ep.initialize
            obj.logTable = table('Size',[0, 1 + obj.ep.nOut],...
                'VariableTypes',repmat({'double'},1,1 + obj.ep.nOut),...
                'VariableNames',[{'Time'}; obj.ep.outputSigName]);
            obj.Massflowbase = Massflowbase;
            obj.m_low = m_low;
            obj.m_high = m_high;
            cd ../
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

