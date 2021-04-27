function plotGroupedBar(data,metricLabels,treatmentLabels,opts)
%PLOTGROUPEDBAR Plot a stacked bar chart of data with error bars
%   'data': A cell array of length(N). Each cell must contain a 2-D array
%   with M columns and O rows. The cell corresponds to a metric and will
%   be given a bar in each cluster of bars. The column corresponds to a
%   treatment, and the column values are the outcome for that treatment
%   across trials. Each row is trial. There will be M clusters of N bars.
%   The height of the bar is the mean of that metric over treatments for
%   that trial, while the error bar is the standard deviation.
%
%   'opts' is an optional struct of options.
%       opts.percentile: if a truthful scalar value, i.e. 0.05, the height of the
%           bar will be the median, and the error bars will be the 0.05 and
%           1 - 0.05 percentiles. Defaults to false TODO: let this take an array.

if (nargin < 4 || isempty(opts))
    opts = struct;
end
if (~isfield(opts,'percentile'))
    opts.percentile = false;
end

N = length(data);
if (N < 1)
    warning('No data provided');
    return
end
M = size(data{1},2);

y = nan(M,N);
if (opts.percentile)
    pos = nan(size(y));
    neg = nan(size(y));
    if (opts.percentile == true || opts.percentile == 1)
        opts.percentile = 0.05; % Defailt to 5 percent
    end
else
    err = nan(size(y));
end

% Compute mean and standard deviation
for i = 1:N
    if (opts.percentile)
        t = prctile(data{i},[opts.percentile,0.5,1-opts.percentile]*100);
        y(:,i) = t(2,:)'; % median value
        pos(:,i) = (t(3,:)-t(2,:))'; % distance between upper percentile and median
        neg(:,i) = (t(2,:)-t(1,:))'; % distance between lower percentile and median
    else
        y(:,i) = mean(data{i})';
        err(:,i) = std(data{i});
    end
end

b = bar(y);
hold on;

if (nargin >= 3)
    set(gca,'XTickLabel',treatmentLabels);
end

% Add error bars; thanks to https://www.mathworks.com/matlabcentral/answers/438514-adding-error-bars-to-a-grouped-bar-plot

nGroups = size(y,1); nBars = size(y,2);
groupwidth = min(0.8, nBars/(nBars + 1.5));
for i = 1:nBars
    x = (1:nGroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nBars);
    if (opts.percentile)
        errorbar(x, y(:,i), neg(:,i), pos(:,i), '.k');
    else
        errorbar(x, y(:,i), err(:,i),'.k');
    end
end

if (isfield(opts,'FaceColor'))
    for i = 1:length(b)
        b(i).FaceColor = opts.FaceColor(i,:);
    end
end

if (opts.percentile)
    barLabel = sprintf('%i%% to %i%%',[opts.percentile 1-opts.percentile]*100);
else
    barLabel = 'Std. Dev.';
end

legend(metricLabels{:},barLabel,'Location','NorthWest');

end

