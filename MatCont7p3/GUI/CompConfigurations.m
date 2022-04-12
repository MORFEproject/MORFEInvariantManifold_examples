classdef CompConfigurations
    %aka BranchManager
    
    
    properties
       %compconfigs = struct();
       conflist = {};
       initiallist = {};
       initials = struct();
       promiscuous = 0;
       
       OrbitConfs = {};
       OrbitConfsHom = {};
    end
    
    methods
        function obj = CompConfigurations()
            obj.OrbitConfs = {};
            
            obj.OrbitConfs{end+1} = SimConf('ode45', 4, true, 0);  %default refine: 4, index: 0
            obj.OrbitConfs{end+1} = SimConf('ode23', 1,true, 1); %default refine: 1
            obj.OrbitConfs{end+1} = SimConf('ode113', 1,true, 2); %default refine: 1
            obj.OrbitConfs{end+1} = SimConf_ode15s(1,true, 3); %default refine1 
            obj.OrbitConfs{end+1} = SimConf('ode23s', 1,true,  4); %default refine: 1
            obj.OrbitConfs{end+1} = SimConf('ode23t', 1,true, 5); %default refine: 1
            obj.OrbitConfs{end+1} = SimConf('ode23tb', 1,true, 6); %default refine: 1
            obj.OrbitConfs{end+1} = SimConf('ode78', 1,false, 7); %default refine: 1           
            obj.OrbitConfs{end+1} = SimConf('ode87', 1,false, 8); %default refine: 1    
            
 
            
            
            

            obj.OrbitConfsHom = {};
            for k = 1:length(obj.OrbitConfs); obj.OrbitConfsHom{end+1} = SimConf_ConnectHom('ConnectionSaddle' , 'Hom', obj.OrbitConfs{k}); end           
            for k = 1:length(obj.OrbitConfs); obj.OrbitConfsHom{end+1} = SimConf_ConnectHSN('ConnectionSaddleNode' , 'HSN', obj.OrbitConfs{k}); end           
            for k = 1:length(obj.OrbitConfs); obj.OrbitConfsHom{end+1} = SimConf_ConnectHet('ConnectionHet' , 'Het', obj.OrbitConfs{k}); end           

            conflist = [obj.OrbitConfs,obj.OrbitConfsHom];
            
            conflist{end+1} = ContConf_EP();
            conflist{end+1} = ContConf_EP_BP();
            conflist{end+1} = ContConf_EP_point('H', @init_H_EP);
            conflist{end+1} = ContConf_EP_point('LP', @init_LP_EP);
            conflist{end+1} = ContConf_EP_point('NE', @init_EP_EP); %neutral saddle equilibrium

            
            conflist{end+1} = ContConf_LP_point('BP', @init_BP_LP);
            conflist{end+1} = ContConf_LP_point('CP', @init_CP_LP);
            conflist{end+1} = ContConf_LP_point('ZH', @init_ZH_LP);
            conflist{end+1} = ContConf_LP_point('BT', @init_BT_LP);
            conflist{end+1} = ContConf_LP_LP(); 
            
            conflist{end+1} = ContConf_H('H', @init_H_H);
            conflist{end+1} = ContConf_H('GH', @init_GH_H);
            conflist{end+1} = ContConf_H('HH', @init_HH_H);
            conflist{end+1} = ContConf_H('ZH', @init_ZH_H);
            conflist{end+1} = ContConf_H('BT', @init_BT_H);
            conflist{end+1} = ContConf_H('NE', @init_H_H); %Neutral Saddle Equilibrium
            
            conflist{end+1} = ContConf_LC_H();
            conflist{end+1} = ContConf_LC_PD();
            conflist{end+1} = ContConf_LC_BPC();
            conflist{end+1} = ContConf_LC_LC(); %LC, BPC and NC (neutral saddle Eq.)
            
            conflist{end+1} = ContConf_LPC('LPC' , @init_LPC_LPC);
            conflist{end+1} = ContConf_LPC('BPC' , @init_BPC_LPC);
            conflist{end+1} = ContConf_LPC('CPC' , @init_CPC_LPC);
            conflist{end+1} = ContConf_LPC('LPNS', @init_LPNS_LPC);
            conflist{end+1} = ContConf_LPC('LPPD', @init_LPPD_LPC);
            conflist{end+1} = ContConf_LPC('R1'  , @init_R1_LPC);
            conflist{end+1} = ContConf_LPC('GPD' , @init_GPD_LPC);
            conflist{end+1} = ContConf_LPC_GH();
            
            conflist{end+1} = ContConf_NS('NS' , @init_NS_NS);
            conflist{end+1} = ContConf_NS('NC' , @init_NS_NS); %Neutral Saddle Cycle
            conflist{end+1} = ContConf_NS('LPNS' , @init_LPNS_NS);
            conflist{end+1} = ContConf_NS('PDNS' , @init_PDNS_NS);
            conflist{end+1} = ContConf_NS('CH' , @init_CH_NS);
            conflist{end+1} = ContConf_NS('R1' , @init_R1_NS);
            conflist{end+1} = ContConf_NS('R2' , @init_R2_NS);
            conflist{end+1} = ContConf_NS('R3' , @init_R3_NS);
            conflist{end+1} = ContConf_NS('R4' , @init_R4_NS);
            conflist{end+1} = ContConf_NS_HH();
            conflist{end+1} = ContConf_NS_ZH();
            
            conflist{end+1} = ContConf_PD('GPD', @init_GPD_PD);
            conflist{end+1} = ContConf_PD('LPPD', @init_LPPD_PD);
            conflist{end+1} = ContConf_PD('PDNS', @init_PDNS_PD);
            conflist{end+1} = ContConf_PD('R2', @init_R2_PD);
            conflist{end+1} = ContConf_PD('PD', @init_PD_PD);
            
            conflist{end+1} = ContConf_BP();
            conflist{end+1} = ContConf_BPC();
            
            conflist{end+1} = ContConf_Hom_Hom();
            conflist{end+1} = ContConf_Hom_LC();
            conflist{end+1} = ContConf_Hom_BT();
            %NCH -> Hom: disabled, init_NCH_Hom has unresolved problems.
%           conflist{end+1} = ContConf_Hom_NCH(); 
%           conflist{end+1} = ContConf_HSN_NCH();
            conflist{end+1} = ContConf_Hom_HTHom();
            conflist{end+1} = ContConf_HTHom();
            conflist{end+1} = ContConf_HSN_LC();
            conflist{end+1} = ContConf_HSN_HSN();
            conflist{end+1} = ContConf_HSN_HTHSN();
            conflist{end+1} = ContConf_HTHSN();
            conflist{end+1} = ContConf_HTHet();
            conflist{end+1} = ContConf_Het('Het', @init_Het_Het);
            conflist{end+1} = ContConf_Het('HTHet', @init_HTHet_Het);
            obj.conflist = conflist;

       
            
            
            obj.sanityCheck();
        end
        
        
        
        
        function sanityCheck(obj)
            labels = cell(1, length(obj.conflist));
            
            for k = 1:length(obj.conflist)
                conf = obj.conflist{k};
                label = conf.getLabel();
                
                assert(~any(strcmp(label, labels)), ['Duplicate: ', label]);
                labels{k} = label;
                
            end
   
            
        end
        
        function list = getConfigs(obj, settings)
            %names = fieldnames(obj.compconfigs);
            list = {};
            for index = 1:length(obj.conflist)
                %name = names{index};
                %conf = obj.compconfigs.(name);
                conf = obj.conflist{index};
                if conf.isAvailable(settings) || obj.promiscuous
                   list{end+1} = conf; 
                end
            end
            ip = settings.IP; label = ip.getLabel();
            
            priorcurve = DefaultValues.getPriority(label);
            if ~obj.promiscuous
                [~, ii] = sort(cellfun(@(x) getCurvePriorityValue(priorcurve, x.getSolutionLabel())*1e10 + x.getPrioritynumber(), list));
                list = list(ii);
            end
        end
        
        
        function s = createSettingTest(obj)
            s = CLSettings();
            s.installSetting('system', CLSystem('../Systems/adapt2.mat'));
            for k = 1:length(obj.conflist)
                cc = obj.conflist{k};
                cc.configureSettings(s);
            end
            s.revealAll();
            s.test();
        end
        
        function filteredlabels = getCurveLabels(obj)
            labels = cellfun(@(x) x.getSolutionLabel(),  obj.conflist, 'un', 0);
            filteredlabels = {};
            for k = 1:length(labels)
                if ~any(strcmp(labels{k}, filteredlabels))
                    filteredlabels{end+1} = labels{k};
                end
            end
            
            
        end
        
        
    end
    
    
    
    
    
end

function i = getCurvePriorityValue(list, label)
    result = find(strcmp(list, label));
    if isempty(result)
       i = length(list) + 1;
    else
       i = result; 
    end
end
