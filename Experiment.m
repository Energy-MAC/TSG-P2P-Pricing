classdef Experiment < handle
    %EXPERIMENT Abstract class for defining experiments
    %   Detailed explanation goes here
    
    properties (SetAccess = protected, GetAccess = public)
        ExperimentParams
    end
    
    properties (SetAccess = private, GetAccess = public)
        CaseFolder
        ExperimentName
        GlobalDataFolder = 'data/';
        NTrials
        ResultsFolder
    end
    
    properties (Access = private)
        RandStreamTrials
        RandStreamExperiment
        tStruct = struct('confoundingVariables',[],'treatmentVariables',[],'outputs',[]); % Model for trial structure; could be another class
        completedTrialsCacheIsDirty = true
        completedTrials % Cache completed trials
    end
    
    properties (Abstract, SetAccess = protected, GetAccess = public)
        Metrics
    end
    
    properties (Dependent)
        CompletedTrials % completed trials, read cached or load from disk.
    end
    
    methods
        
        % Constructor: experimentName and globalDataFolder as path to data
        % directory are required. A caseName is optional to specify a
        % particular case of this experiment.
        function obj = Experiment(experimentName,args)
            if (nargin < 2)
                args = struct;
            end
            
            obj.ExperimentName = experimentName;
            
            % Set up directory structure
            if ~isfield(args,'caseName')
                obj.ResultsFolder = sprintf('%sexperiments/outputs/%s/',obj.GlobalDataFolder,obj.ExperimentName);
                obj.CaseFolder = sprintf('%sexperiments/inputs/%s/',obj.GlobalDataFolder,obj.ExperimentName);
            else
                caseName = args.caseName;
                obj.ResultsFolder = sprintf('%sexperiments/outputs/%s/%s/',obj.GlobalDataFolder,obj.ExperimentName,caseName);
                obj.CaseFolder = sprintf('%sexperiments/inputs/%s/%s/',obj.GlobalDataFolder,obj.ExperimentName,caseName);
            end
                
            if ~exist(obj.ResultsFolder,'dir')
                mkdir(obj.ResultsFolder);
            end
            if ~exist(sprintf('%strials/',obj.ResultsFolder),'dir')
                mkdir(sprintf('%strials/',obj.ResultsFolder));
            end
            
            % Get the number of trials to run
            p = readKeyValue(sprintf('%skey_values.csv',obj.CaseFolder));
            obj.NTrials = p.N_trials;
            
            % Set up the random number generation and streams for each
            % trial
            seed = p.seed;
            streams = RandStream.create('mrg32k3a','Seed',seed,'NumStreams',obj.NTrials+1,'CellOutput',true);
            obj.RandStreamExperiment = streams{1}; % First one is for experiment setup
            obj.RandStreamTrials = streams(2:end); % Second is for each trial
            
            % Remove those properties and set the rest as experiment
            % parameters
            p = rmfield(p,'N_trials');
            p = rmfield(p,'seed');
            obj.ExperimentParams = p;
            
            % Seed rng for setting up additional parameters
            RandStream.setGlobalStream(obj.RandStreamExperiment);
            if (isfield(args,'setupArgs'))
                setupArgs = args.setupArgs;
            else
                setupArgs = struct;
            end
            obj.setupAdditionalParameters(setupArgs);
        end
        
        % Run the experiment either with trials in parallel or series
        function trials = runExperiment(obj,order)
            if nargin < 2
                order = 'par';
            end
            switch lower(order)
                case 'par'
                    trials = obj.runTrialsPar(true);
                case 'ser'
                    trials = obj.runTrialsSerFrom();
                otherwise
                    error('Unrecognized execution order');
            end
        end
        
        function trials = runTrialsSerFrom(obj,startTrial,endTrial)
            if nargin < 2 || isempty(startTrial)
                startTrial = 1;
            end
            if nargin < 3 || isempty(endTrial)
                endTrial = obj.NTrials;
            end
            
            fprintf('Running %s in series starting from trial %i\r',obj.ExperimentName,startTrial);
            if (startTrial > 1)
                % Get trials that were previously saved
            end
            
            for i = startTrial:endTrial
%                try
                    t = obj.runTrial(i);
                    if (isempty(t.outputs))
                        warning('A completed trial should assign outputs')
                    end
                    obj.saveTrial(i,t);
 %               catch ME
  %                  warning('An error was encountered running trial %i. Skipping it. Error message was %s',i,ME.message);
   %             end
            end
            trials = obj.loadAllTrials();
        end
        
        function trials = runTrialsPar(obj,overwrite)
            % Get an array of trials struct
            fprintf('Running %s in parallel with %i trials. Overwriting any saved: %s\r',obj.ExperimentName,obj.NTrials,mat2str(overwrite));
            
            % Model for trial structure
            tStruct = obj.tStruct;
            
            % Functions to call within par loop
            loadTrialFun = @obj.loadTrial;
            saveTrialFun = @obj.saveTrial;
            runTrialFun = @obj.runTrial;
            
            % Run trials
            parfor i = 1:obj.NTrials
                if (overwrite)
                    % Create a fresh one
                    trial = tStruct;
                else
                    % Try and load
                    trial = feval(loadTrialFun,i);
                    if isempty(trial)
                        % Couldn't load a saved one so create fresh
                        trial = tStruct;
                    end
                end
                
                if (isempty(trial.outputs)) % Will be true for fresh trials
                    try
                        t = feval(runTrialFun,i);
                        if (isempty(t.outputs))
                            warning('A completed trial should assign outputs')
                        end
                        trial = t;
                        feval(saveTrialFun,i,trial);
                    catch ME
                        warning('An error was encountered running trial %i. Skipping it. Error message was %s',i,ME.message);
                    end
                else
                    fprintf('Trial %i already completed\n',i);
                end
            end
            fprintf('Completed all trials. Loading results...');
            trials = obj.loadAllTrials();
            fprintf('...done\n');
        end
        
        function trial = runTrial(obj,i)
            fprintf('Starting trial %i...\n',i);
            % Requires that obj.loadExperimentParameters has been called
            % first
            
            % Set the stream for this trial
            RandStream.setGlobalStream(obj.RandStreamTrials{i});
            
            tTrial = tic; % Time the trial
            trial = struct;
            fprintf('Generating confounding variables for trial %i...\n',i);
            tDataProcess = tic;
            trial.confoundingVariables = obj.generateConfoundingVariables(i); % Returns struct
            fprintf('...confounding variables generated for trial %i. Data processing took %g seconds\n',i,toc(tDataProcess));
            trial.treatmentVariables = obj.generateTreatmentVariables(trial.confoundingVariables); % Returns struct array
            stream = RandStream.getGlobalStream;
            streamState = stream.State;
            for j = 1:length(trial.treatmentVariables)
                stream.State = streamState; % Ensures same RNG for each trial
                fprintf('Simulating treatment %i of %i for trial %i...\n',j,length(trial.treatmentVariables),i);
                % Simulation
                tSim = tic; % Time the simulation with the treatment
                trial.outputs(j) = obj.simulateTreatment(trial.confoundingVariables,trial.treatmentVariables(j));
                fprintf('...done simulating treatment %i of %i for trial %i. Elapsed time: %g seconds\n',j,length(trial.treatmentVariables),i,toc(tSim));
            end
            fprintf('...trial %i completed. Elapsed time for trial: %g seconds\n',i,toc(tTrial));
        end
        
        function trial = loadTrial(obj,i)
            try
                t = load(sprintf('%strials/%i',obj.ResultsFolder,i),'trial');
                trial = t.trial;
            catch
                trial = [];
            end
        end
        
        function [trials,indCompleted] = loadAllTrials(obj)
            indCompleted = false(1,obj.NTrials);
            for i = 1:obj.NTrials
                t = obj.loadTrial(i);
                if (~isempty(t) && ~isempty(t.outputs))
                    trials(i) = t;
                    indCompleted(i) = true;
                else
                    trials(i) = obj.tStruct;
                end
            end
        end
        
        function [trials] = loadCompletedTrials(obj)
            [trials,indCompleted] = obj.loadAllTrials();
            trials = trials(indCompleted);
            obj.completedTrials = trials;
            obj.completedTrialsCacheIsDirty = false;
        end
        
        function saveTrial(obj,i,trial)
            save(sprintf('%strials/%i',obj.ResultsFolder,i),'trial');
            obj.completedTrialsCacheIsDirty = true;
        end
        
        function resultsDistStruct = computeResultsDistributionMetricMultiple(obj,metricNames)
            if nargin >= 2
                % Find indices of metricNames
                metricInds = [];
                for m = 1:length(metricNames)
                    metricInds = [metricInds find(strcmp(metricNames{m},obj.Metrics(:,1)),1)];
                end
            else
                % Do them all
                metricInds = 1:size(obj.Metrics,1);
            end
            trials = obj.CompletedTrials;
            
            % Infer the number of trials and treatments in each trial
            Ntrials = length(trials); % Number of trials
            if ~Ntrials
                % resultsDist = []; % Think this can be deleted unless
                % function needs to return an output. TODO: verify
                return;
            end
            Ntreats = length(trials(1).treatmentVariables); % Number of treatments; assumes all trials have the same treatments
            
            % Compute each metric
            for m = metricInds
                metricFun = obj.Metrics{m,2};
                % Preallocate a matrix for the results for each trial and
                % treatment. It will be Ntrials x Ntreats and the i,j element
                % will be a performance metric for that combination
                resultsDist = nan(Ntrials,Ntreats);
                % Compute the performance metric for each outcome
                for i = 1:Ntrials
                    if (length(trials(i).outputs) < Ntreats)
                        % Not all outputs were computed for each treatment for this trial
                        continue
                    end
                    for j = 1:Ntreats
                        % Compute the metric for that treatment
                        resultsDist(i,j) = metricFun(trials(i).outputs(j),obj.ExperimentParams,trials(i));
                    end
                end
                resultsDistStruct.(obj.Metrics{m,1}) = resultsDist;
            end
        end
        
        function resultsDist = computeResultsDistributionMetricSingle(obj,metricFun,passTrial)
            if (nargin < 3)
                passTrial = false; % Optional flag to request that the whole trial info (i.e. confounding and experiment params) be passed to compute the metric.
            end
            
            trials = obj.CompletedTrials;
            
            % Infer the number of trials and treatments in each trial
            Ntrials = length(trials); % Number of trials
            if ~Ntrials
                resultsDist = [];
                return;
            end
            Ntreats = length(trials(1).treatmentVariables); % Number of treatments; assumes all trials have the same treatments
            
            % Preallocate a matrix for the results for each trial and
            % treatment. It will be Ntrials x Ntreats and the i,j element
            % will be a performance metric for that combination
            resultsDist = nan(Ntrials,Ntreats);
            
            % Compute the performance metric for each outcome
            for i = 1:Ntrials
                if (length(trials(i).outputs) < Ntreats)
                    % Not all outputs were computed for each treatment for this trial
                    continue
                end
                for j = 1:Ntreats
                    % Compute the metric for that treatment
                    if (passTrial)
                        resultsDist(i,j) = metricFun(trials(i).outputs(j),obj.ExperimentParams,trials(i));
                    else
                        resultsDist(i,j) = metricFun(trials(i).outputs(j),obj.ExperimentParams);
                    end
                end
            end
        end
        
        function plotMetricConvergence(obj,metricName,params)
            % Plot the evolution of the coefficient of variation for each
            % treatment as the number of trials increase
            resultsDistStruct = computeResultsDistributionMetricMultiple(obj,{metricName});
            results = resultsDistStruct.(metricName);
            
            if (isfield(params,'treatmentInds'))
                treatmentInds = params.treatmentInds;
            else
                treatmentInds = 1:size(results,2);
            end
            if (isfield(params,'treatmentLabels'))
                treatmentLabels = params.treatmentLabels;
            else
                treatmentLabels = arrayfun(@(i) num2str(i),treatmentInds,'UniformOutput',false);
            end
            if (isfield(params,'metricLabel'))
                metricLabel = params.metricLabel;
            else
                metricLabel = metricName;
            end
            
            cov = nan(size(results,1)-1,length(treatmentInds));
            for i = 2:size(results,1)
                cov(i-1,:) = std(results(1:i,treatmentInds))./mean(results(1:i,treatmentInds));
            end
            plot(2:size(results,1),cov);
            if (length(treatmentLabels) > 1)
                legend(treatmentLabels{:});
            end
            ylabel('CoV');
            xlabel('Trials');
            title(sprintf('Metric variation: %s',metricLabel));
        end
        
        function [axs,f] = plotResultsDistributionTreatment(obj,metricName,treatmentLabels)
            % metricName: string defining the metric. Must be in
            % obj.Metrics.
            % treatmentLabels: (optional) cell array mapping the index of
            % the treatment to a label to be displayed in the subplot
            % title. If defined, must be of length Ntreatements.
            
            resultsDistStruct = computeResultsDistributionMetricMultiple(obj,{metricName});
            results = resultsDistStruct.(metricName);
            
            N = size(results,2); % Number of treatments
            treatmentLabelsDefined = nargin > 2 && length(treatmentLabels) == N;
            
            if (N > 12)
                warning('Only showing max of 12 treatments');
                N = 12;
            end
            
            if ~treatmentLabelsDefined
                treatmentLabels = arrayfun(@(i) sprintf('Treatment %i',i),1:N,'UniformOutput',false);
            end
            
            % Get number of rows
            if (N < 1)
                return;
            elseif (N <= 3)
                nRows = 1;
            elseif (N <= 8)
                nRows = 2;
            else
                nRows = 3;
            end
            nCols = ceil(N/nRows);
            
            % Use the same bin width for all
            if size(results,3) >= 60
                nBins = 20;
            else
                nBins = ceil(size(results,1)/3);
            end
            edges = linspace(min(results(:)),max(results(:)),nBins+1);
            
            f = figure('Name',metricName);
            yMax = 0;
            for j = 1:N
                axs(j) = subplot(nRows,nCols,j);
                histogram(results(:,j),edges,'Normalization','probability');
                title(treatmentLabels{j})
                t = get(axs(j),'YLim');
                yMax = max(yMax,t(2));
            end
            try
                sgtitle(sprintf('Distribution of %s for each treatment',metricName));
            end
            
            set(axs,'YLim',[0 yMax]); % Give all the same y axis
            
            f = figure('Name',sprintf('%s (overlay)',metricName));
            for j = 1:N
                histogram(results(:,j),edges,'Normalization','pdf');
                hold on
            end
            legend(treatmentLabels);
        end
        
        function plotResultsDistributionMetrics(obj,params)
        %function [ax,f] = plotRelativeResultsDistributionMetrics(obj,metricNames,baseInd,treatmentLabels,metricLabels)
            % Makes a stacked bar chart. Each group of bars corresponds to
            % a metric. Each bar in the group is the percent change in the
            % metric for that treatment relative to the baseline treatment.
            
            if (nargin < 2)
                params = struct;
            end
            
            if (isfield(params,'metricNames'))
                metricNames = params.metricNames;
                resultsDistStruct = computeResultsDistributionMetricMultiple(obj,metricNames);
                if length(fieldnames(resultsDistStruct)) < length(metricNames)
                    error('At least one of the metric names is not defined');
                end
            else
                % Use all metrics (infer from results field names)
                resultsDistStruct = computeResultsDistributionMetricMultiple(obj);
                metricNames = fieldnames(resultsDistStruct);
                if (length(metricNames)) < 1
                    warning('No metrics defined for experiment');
                    return;
                end
            end
            % N metric names
            N = length(metricNames);
            
            % Infer metric labels if not explicitly defined
            if (~isfield(params,'metricLabels') || length(params.metricLabels) ~= N)
                metricLabels = metricNames;
            else
                metricLabels = params.metricLabels;
            end
            
            
            % Infer treatments, assumes the same for all trials and metrics
            % M treatments
            if (isfield(params,'treatmentInd'))
                tInd = params.treatmentInd;
                M = length(tInd);
            else
                M = size(resultsDistStruct.(metricNames{1}),2);
                tInd = 1:M;
            end
            
            if (isfield(params,'treatmentLabels') && length(params.treatmentLabels) == M)
                treatmentLabels = params.treatmentLabels;
                treatmentsNumeric = isnumeric(treatmentLabels);
            else
                treatmentLabels = [];
                treatmentsNumeric = false;
            end
            
            
            if (treatmentsNumeric)
                % Treatments have a numeric interpretation. Sort them
                [treatmentLabels, tInd] = sort(treatmentLabels);
            end
            
            % Set normalization, default 'absolute'
            if (isfield(params,'normalization'))
                normalization = params.normalization;
            else
                normalization = 'absolute';
            end
            
            data = cell(1,N);
            switch lower(normalization)
                case 'absolute'
                    for i = 1:N
                        t = resultsDistStruct.(metricNames{i});
                        data{i} = t(:,tInd);
                    end
                    ylabelText = 'Metric value';
                case 'percentchangebytrial'
                    if (isfield(params,'baseIndex'))
                        baseInd = params.baseIndex;
                    else
                        baseInd = 1;
                    end
                    
                    baseInd2 = tInd == baseInd;
                    
                    % Set the y label to indicate percent change
                    ylabelText = sprintf('%% change relative to treatment: %s',treatmentLabels{baseInd2});
                    
                    % Remove the baseInd
                    tInd(baseInd2) = [];
                    treatmentLabels(baseInd2) = [];
                    % Function to compute relative change of each treatment to the
                    % reference.
                    relativeChange = @(x) (x-repmat(x(:,baseInd),1,size(x,2)))./abs(repmat(x(:,baseInd)+1e-10,1,size(x,2)))*100;
                    for i = 1:N
                        t = relativeChange(resultsDistStruct.(metricNames{i}));
                        data{i} = t(:,tInd);
                    end
                otherwise
                    error('Unrecognized normalization: %s',normalization);
            end
            
            if (treatmentsNumeric)
                % Make a line graph of the data
                plotMetricLine(data,treatmentLabels,metricLabels);
            else
                % Treatments are labeled categorical
                if (isfield(params,'barOpts'))
                    barOpts = params.barOpts;
                else
                    barOpts = struct;
                end
                plotGroupedBar(data,metricLabels,treatmentLabels,barOpts);
            end
            
            
            title(gca,'Metrics');
            ylabel(gca,ylabelText);
            xlabel(gca,'Treatment');
        end
        
        function y = get.CompletedTrials(self)
            if (self.completedTrialsCacheIsDirty)
                y = self.loadCompletedTrials();
            else
                y = self.completedTrials;
            end
        end
        
        function setupAdditionalParameters(obj,setupArgs)
        end
    end
    
    methods (Abstract)
        
        % Generate confounding variables
        generateConfoundingVariables(obj,trialInd)
        
        % Generate treatment variables
        generateTreatmentVariables(obj,confoundingVariables)
        
        % Simulate a trial
        simulateTreatment(obj,confoundingVariables,treatmentVariables)
    end
end

