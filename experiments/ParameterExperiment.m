classdef ParameterExperiment < Experiment
    %CONTROLLERPERFORMANCEEXPERIMENT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = protected, GetAccess = public)
        % Pair of metricName, metricFun.
        Metrics = defineConvergenceExperimentsMetrics()
    end
    
    methods
        function obj = ParameterExperiment(caseName)
            args = struct;
            if nargin > 0
                args.caseName = caseName;
            end
            obj@Experiment('parameter',args);
            
        end
        
        % Programmatically read / generate additional parameters
        function setupAdditionalParameters(obj,setupArgs)
            
            % Read empirical load profiles
            loadProfiles = table2array(readtable(...
                sprintf('%scommon/load_profiles_2017.csv',...
                obj.GlobalDataFolder),'ReadVariableNames',true));
            loadProfiles = loadProfiles(:,2:end); % Remove hour index col
            
            obj.ExperimentParams.quantities = loadProfiles; % TxN load profiles
        end
        
        % Generate confounding variables, in the sense of multiple
        % 'control' groups across trials. These are generated for each
        % trial.
        function confoundingVars = generateConfoundingVariables(obj,trialInd)
            p = obj.ExperimentParams;
            confoundingVars = struct;
            
            t = randi(size(p.quantities,1)); % random time
            playerInds = randi(size(p.quantities,2),[1,2]); % 2 random players (with replacement)
            % Grab empirical price, quantity, and pv
            empiricalQuantity = p.quantities(t,playerInds);
            empiricalQuantity = max(empiricalQuantity,0.1); % Fit utility requires d > 0
            % Sample elasticities uniformly in interval
            elasticities = (p.maxElasticity-p.minElasticity)*rand(2,1) + p.minElasticity;
            % Constant price
            empiricalPrice = p.nominalPrice*ones(size(empiricalQuantity));
            % Fit utility
            u = fit_utility_N(elasticities,empiricalPrice,empiricalQuantity);
            
            % Solar power determined as random scale factor, uniformly
            % distributed in [0,maxSolarScaleFactor] times the load at that time
            confoundingVars.solarPower = p.maxSolarScaleFactor*rand(1,2).*empiricalQuantity;
            confoundingVars.playerInds = playerInds;
            confoundingVars.elasticities = elasticities;
            confoundingVars.utilityFunctions = u;
        end
        
        % Generate treatment variables. These are enumerated and compared
        % within each trial.
        function treatmentVars = generateTreatmentVariables(obj,confoundingVars)
            p = obj.ExperimentParams;
            treatmentVars = struct;
            gammaSample = p.minGamma:p.gammaStepSize:p.maxGamma;
            delta0Sample = p.minDelta0:p.delta0StepSize:p.maxDelta0;
            
            for i = 1:length(gammaSample)
                for j = 1:length(delta0Sample)
                    treatmentVars(i,j).gamma = gammaSample(i);
                    treatmentVars(i,j).delta0 = delta0Sample(j);
                end
            end
            treatmentVars = treatmentVars(:);
        end
        
        % Simulate a trial
        function results = simulateTreatment(obj,confoundingVars,treatmentVars)
            
            % In this experiement, we always have N=2 and T=1 with the
            % pi-player the second player
            u = confoundingVars.utilityFunctions;
            piInd = 2;
            p_s = confoundingVars.solarPower;
            delta_T = 1;
            S_max = zeros(1,2);
            P_b_max = zeros(1,2);
            s0 = zeros(1,2);
            delta0 = treatmentVars.delta0;
            gamma = treatmentVars.gamma;
            epsilon = obj.ExperimentParams.epsilon;
            maxIterations = obj.ExperimentParams.maxIterations;
            
            results = struct;
            
            c = struct; % Centralized solution results
            [c.price, c.q, c.d, c.p_b, c.s] = opt_centralized(...
                u,p_s,delta_T,P_b_max,S_max,s0);
            
            b = struct; % Bidding process results
            b.trajectories = struct; % Trajectory variables (3D arrays)
            [b.price, b.q, b.numIterations, b.trajectories.price,...
                b.trajectories.qFeas, b.trajectories.qReq,...
                b.trajectories.delta, b.d, b.priceQ, b.p_b, b.s] = ...
                bid_process(u,piInd,p_s,delta_T,P_b_max,S_max,s0,delta0,gamma,epsilon,maxIterations);
            
            results.centralized = c;
            results.bidding = b;
        end
        
    end
    
end