classdef StorageExperiment < Experiment
    %CONTROLLERPERFORMANCEEXPERIMENT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = protected, GetAccess = public)
        % Pair of metricName, metricFun.
        Metrics = defineStorageExperimentsMetrics()
    end
    
    methods
        function obj = StorageExperiment(caseName)
            %CONTROLLERPERFORMANCE Construct an instance of this class
            args = struct;
            if nargin > 0
                args.caseName = caseName;
            end
            obj@Experiment('storage',args);
            
        end
        
        % Programmatically read / generate additional parameters
        function setupAdditionalParameters(obj,setupArgs)
            
            % Read sets for varying battery capacity, battery charge rate,
            % solar scaling, and time horizon
            t = readtable(sprintf('%sbatteryCapacitySet.csv',obj.CaseFolder));
            batteryCapacitySet = t{:,1};
            t = readtable(sprintf('%sbatteryChargeRateSet.csv',obj.CaseFolder));
            batteryChargeRateSet = t{:,1};
            %t = readtable(sprintf('%ssolarScaleSet.csv',obj.CaseFolder));
            %solarScaleSet = t{:,1};
            t = readtable(sprintf('%stimeHorizonSet.csv',obj.CaseFolder));
            timeHorizonSet = t{:,1};
            
            % Read empirical solar, load and price profiles
            load_data_file = sprintf('%scommon/load_profiles_2017.csv',obj.GlobalDataFolder);
            pv_data_file = sprintf('%scommon/solar_pv_profiles_2017.csv',obj.GlobalDataFolder);
            t1 = table2array(readtable(load_data_file));
            t2 = table2array(readtable(pv_data_file));
            loadShape = t1(2:end, 2:end);
            solarShape = t2(2:end, 2:end);
            % Same profile each day
            prices = repmat([ones(10,1)*0.1; ones(6,1)*0.15; ones(5,1)*0.3; ones(3,1)*0.1],365,size(solarShape,2));
            prices = prices(1:size(solarShape, 1), 1:end);

            
            obj.ExperimentParams.batteryCapacitySet = batteryCapacitySet;
            obj.ExperimentParams.batteryChargeRateSet = batteryChargeRateSet;
            %obj.ExperimentParams.solarScaleSet = solarScaleSet;
            obj.ExperimentParams.timeHorizonSet = timeHorizonSet;
            
            obj.ExperimentParams.solarShape = solarShape; % TxN solar profiles
            obj.ExperimentParams.quantities = loadShape; % TxN load profiles
            obj.ExperimentParams.prices = prices; % TxN ToU price profiles
        end
        
        % Generate confounding variables, in the sense of multiple
        % 'control' groups across trials. These are generated for each
        % trial.
        function confoundingVars = generateConfoundingVariables(obj,trialInd)
            p = obj.ExperimentParams;
            confoundingVars = struct;
            
            % Select number of players randomly
            maxN = p.maxN;
            if ~any(p.solarShape(:))
                error('Must be some solar capacity');
            end
            % Select players randomly. Keep trying until at least one
            % player has some solar.
            while true
                N = randi([2,maxN]);
                % TODO: set start time randomly as well and have horizon
                % wrap around (wrap around can be solved by having data
                % span at least twice the max time
                max_T = max(p.timeHorizonSet);
                % Pick initial hour of the year
                T_init = randi([1, size(p.solarShape, 1) - (max_T + 1)]);
                % Pick random Horizon from Horizon Set
                T = p.timeHorizonSet(randi(length(p.timeHorizonSet)));
                % Define Range of Time
                T_range = T_init:(T_init + T - 1);
                playerInds = randsample(size(p.quantities,2),N);
                solarShape = p.solarShape(T_range,playerInds);
                if any(solarShape(:))
                    break;
                end
            end
            
            % Grab empirical price, quantity, and pv
            % Assume empirical price the same for everyone
            empiricalPrice = p.prices(T_range,playerInds);
            empiricalQuantity = p.quantities(T_range,playerInds);
            % Fit utility requires d > 0:
            empiricalQuantity = max(empiricalQuantity,0.1); 
            % Sample elasticities uniformly in interval
            elasticities = (p.maxElasticity-p.minElasticity)*rand(N,1) + p.minElasticity;
            % Fit utility
            u = fit_utility_N(elasticities,empiricalPrice,empiricalQuantity);
            
            % Normalize the solar profile to sum to the empirical demand
            %solarShape = solarShape / sum(solarShape(:));
            solarShape = solarShape / sum(empiricalQuantity(:));
            % Randomly distribute battery capacity summing to 1
            batteryAssignment = rand(1,N);
            batteryAssignment = batteryAssignment/sum(batteryAssignment);
            
            %Set Confounding Variables
            confoundingVars.piInd = randi(N);
            confoundingVars.solarShape = solarShape;
            confoundingVars.batteryAssignment = batteryAssignment;
            confoundingVars.batterySoC0 = 0.5 + 0.5*rand(1,N);
            confoundingVars.playerInds = playerInds;
            confoundingVars.elasticities = elasticities;
            confoundingVars.utilityFunctions = u;
        end
        
        % Generate treatment variables. These are enumerated and compared
        % within each trial.
        function treatmentVars = generateTreatmentVariables(obj,confoundingVars)
            batteryCapacitySet = obj.ExperimentParams.batteryCapacitySet;
            batteryChargeRateSet = obj.ExperimentParams.batteryChargeRateSet;
            %solarScaleSet = obj.ExperimentParams.solarScaleSet;
            
            treatmentVars = struct;
            for i = 1:length(batteryCapacitySet)
                for j = 1:length(batteryChargeRateSet)
                    treatmentVars(i,j).batteryCapacity = batteryCapacitySet(i);
                    treatmentVars(i,j).batteryChargeRate = batteryChargeRateSet(j);
                end
            end
            treatmentVars = treatmentVars(:);
        end
        
        % Simulate a trial
        function results = simulateTreatment(obj,confoundingVars,treatmentVars)
            
            u = confoundingVars.utilityFunctions;
            piInd = confoundingVars.piInd;
            p_s = confoundingVars.solarShape;
            delta_T = obj.ExperimentParams.deltaT;
            S_max = treatmentVars.batteryCapacity*confoundingVars.batteryAssignment;
            P_b_max = treatmentVars.batteryChargeRate*S_max;
            s0 = confoundingVars.batterySoC0.*S_max;
            delta0 = obj.ExperimentParams.delta0; 
            gamma = obj.ExperimentParams.gamma;
            epsilon = obj.ExperimentParams.epsilon;
            maxIterations = obj.ExperimentParams.maxIterations;
            
            results = struct;
            
            c = struct; % Centralized solution results
            [c.price, c.q, c.d, c.p_b, c.s, c.lambda_b, c.lambda_c] = opt_centralized(...
                u,p_s,delta_T,P_b_max,S_max,s0, 0);
            
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