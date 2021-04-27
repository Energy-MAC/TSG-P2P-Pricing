configParams = config();

addpath('code');
addpath('experiments');
addpath('bin');

% Import computational experiment superclass and utilities
addpath(configParams.experimentDir);
addpath(sprintf('%s/%s',configParams.experimentDir,'utils'));