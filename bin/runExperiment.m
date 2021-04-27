function expObj = runExperiment(experimentType,caseName,startTrial,endTrial)

disp('Starting...');

switch (lower(experimentType))
    case 'convergence'
        constructor = @ConvergenceExperiment;
    case 'parameter'
        constructor = @ParameterExperiment;
    otherwise
        error('Unsupported experiment type: ''%s''',experimentType);
end

if (nargin < 3)
    startTrial = '';
end
if (nargin < 4)
    endTrial = '';
end
        
if (ischar(startTrial))
    startTrial = str2double(startTrial);
end
if (ischar(endTrial))
    endTrial = str2double(endTrial);
end

if (isnan(startTrial))
    startTrial = 1;
end
if (isnan(endTrial))
    endTrial = [];
end

if nargin < 2 || isempty(caseName)
    fprintf('Running default case...\r');
    expObj = constructor();
else
    expObj = constructor(caseName);
    fprintf('Running case ''%s''...\r',caseName);
end

fprintf('Running trials %i to %i...\r',startTrial,endTrial);
expObj.runTrialsSerFrom(startTrial, endTrial);
disp('Finished!');
end

