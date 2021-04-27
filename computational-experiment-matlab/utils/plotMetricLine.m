function [f,ax] = plotMetricLine(data,treatments,metricLabels)
%PLOTGROUPEDBAR Plot a stacked bar chart of data with error bars
%   'data': A cell array of length(N). Each cell must contain a 2-D array
%   with M columns and O rows. The cell corresponds to a metric and will
%   be given a line. The column corresponds to a treatment is plotted as a
%   series of values against the 'treatment' input. Each row is trial.
%   Two lines will be plotted for each n in N: the mean over trials, and
%   the standard deviation over trials.

N = length(data);
if (N < 1)
    warning('No data provided');
    return
end
M = size(data{1},2);

y = nan(M,N);
err = nan(size(y));

% Compute mean and standard deviation
for i = 1:N
    y(:,i) = mean(data{i})';
    err(:,i) = std(data{i});
end

f = figure;
plot(treatments,y);
hold on;
ax = gca;
set(ax,'ColorOrderIndex',1);
plot(treatments,y+err,'--');
set(ax,'ColorOrderIndex',1);
plot(treatments,y-err,'--');
if (nargin >= 3 && length(metricLabels) > 1)
    legend(metricLabels{:})
end
end

