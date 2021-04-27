classdef ConvergenceExperiment < Experiment
    %CONTROLLERPERFORMANCEEXPERIMENT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = protected, GetAccess = public)
        % Pair of metricName, metricFun.
        Metrics = defineConvergenceExperimentsMetrics()
    end
    
    methods
        function obj = ConvergenceExperiment(caseName)
            %CONTROLLERPERFORMANCE Construct an instance of this class
            args = struct;
            if nargin > 0
                args.caseName = caseName;
            end
            obj@Experiment('convergence',args);
            
        end
        
        % Programmatically read / generate additional parameters
        function setupAdditionalParameters(obj,setupArgs)
            
            % Read sets for varying battery capacity, battery charge rate,
            % solar scaling, and time horizon
            t = readtable(sprintf('%sbatteryCapacitySet.csv',obj.CaseFolder));
            batteryCapacitySet = t{:,1};
            t = readtable(sprintf('%sbatteryChargeRateSet.csv',obj.CaseFolder));
            batteryChargeRateSet = t{:,1};
            t = readtable(sprintf('%ssolarScaleSet.csv',obj.CaseFolder));
            solarScaleSet = t{:,1};
            t = readtable(sprintf('%stimeHorizonSet.csv',obj.CaseFolder));
            timeHorizonSet = t{:,1};
            
            
            % Read empirical load and solar profiles
            loadProfiles = table2array(readtable(...
                sprintf('%scommon/load_profiles_2017.csv',...
                obj.GlobalDataFolder),'ReadVariableNames',true));
            loadProfiles = loadProfiles(:,2:end); % Remove hour index col
            solarProfiles = table2array(readtable(...
                sprintf('%scommon/solar_pv_profiles_2017.csv',...
                obj.GlobalDataFolder),'ReadVariableNames',true));
            solarProfiles = solarProfiles(:,2:end); % Remove hour index col

            
            obj.ExperimentParams.batteryCapacitySet = batteryCapacitySet;
            obj.ExperimentParams.batteryChargeRateSet = batteryChargeRateSet;
            obj.ExperimentParams.solarScaleSet = solarScaleSet;
            obj.ExperimentParams.timeHorizonSet = timeHorizonSet;
            
            obj.ExperimentParams.solarProfiles = solarProfiles; % TxN solar profiles
            obj.ExperimentParams.quantities = loadProfiles; % TxN load profiles
        end
        
        % Generate confounding variables, in the sense of multiple
        % 'control' groups across trials. These are generated for each
        % trial.
        function confoundingVars = generateConfoundingVariables(obj,trialInd)
            p = obj.ExperimentParams;
            confoundingVars = struct;
            
            % Select number of players and number of time steps randomly
            maxN = p.maxN;
            if ~any(p.solarProfiles(:))
                error('Must be some solar capacity somewhere');
            end
            T = p.timeHorizonSet(randi(length(p.timeHorizonSet)));
                N = randi([2,maxN]);
            % Select player and time blocks randomly.
            % Keep trying until at least one player has some solar.
            while true
                playerInds = randsample(size(p.quantities,2),N);
                startTime = randi(size(p.quantities,1)-T) + 1;
                solarShape = p.solarProfiles(startTime:startTime+T-1,playerInds);
                if any(solarShape(:))
                    break;
                end
            end
            
            % Grab nominal price, quantity, and pv
            empiricalPrice = p.nominalPrice*ones(T,N);
            empiricalQuantity = p.quantities(startTime:startTime+T-1,playerInds);
            empiricalQuantity = max(empiricalQuantity,0.1); % Fit utility requires d > 0
            % Sample elasticities uniformly in interval
            elasticities = (p.maxElasticity-p.minElasticity)*rand(N,1) + p.minElasticity;
            % Fit utility
            u = fit_utility_N(elasticities,empiricalPrice,empiricalQuantity);
            
            % Normalize the solar profile to sum to the empirical demand
            solarShape = solarShape / sum(solarShape(:)) * sum(empiricalQuantity(:));
            % Randomly distribute battery capacity summing to aggregate
            % mean load
            batteryAssignment = rand(1,N);
            batteryAssignment = batteryAssignment/ sum(batteryAssignment) * sum(mean(empiricalQuantity));
            
            % Randomly sample gamma uniformly in interval
            gamma = (p.maxGamma-p.minGamma)*rand(1) + p.minGamma;
            
            confoundingVars.piInd = randi(N);
            confoundingVars.solarShape = solarShape;
            confoundingVars.batteryAssignment = batteryAssignment;
            confoundingVars.batterySoC0 = rand(1,N);
            confoundingVars.playerInds = playerInds;
            confoundingVars.elasticities = elasticities;
            confoundingVars.utilityFunctions = u;
            confoundingVars.gamma = gamma;
        end
        
        % Generate treatment variables. These are enumerated and compared
        % within each trial.
        function treatmentVars = generateTreatmentVariables(obj,confoundingVars)
            batteryCapacitySet = obj.ExperimentParams.batteryCapacitySet;
            batteryChargeRateSet = obj.ExperimentParams.batteryChargeRateSet;
            solarScaleSet = obj.ExperimentParams.solarScaleSet;
            
            treatmentVars = struct;
            for i = 1:length(batteryCapacitySet)
                for j = 1:length(batteryChargeRateSet)
                    for k = 1:length(solarScaleSet)
                        treatmentVars(i,j,k).batteryCapacity = batteryCapacitySet(i);
                        treatmentVars(i,j,k).batteryChargeRate = batteryChargeRateSet(j);
                        treatmentVars(i,j,k).solarScale = solarScaleSet(k);
                    end
                end
            end
            treatmentVars = treatmentVars(:);
        end
        
        % Simulate a trial
        function results = simulateTreatment(obj,confoundingVars,treatmentVars)
            
            u = confoundingVars.utilityFunctions;
            piInd = confoundingVars.piInd;
            p_s = treatmentVars.solarScale*confoundingVars.solarShape;
            delta_T = obj.ExperimentParams.deltaT;
            S_max = treatmentVars.batteryCapacity*confoundingVars.batteryAssignment;
            P_b_max = treatmentVars.batteryChargeRate*S_max;
            s0 = confoundingVars.batterySoC0.*S_max;
            delta0 = obj.ExperimentParams.delta0; % Should be confounding?
            gamma = confoundingVars.gamma;
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