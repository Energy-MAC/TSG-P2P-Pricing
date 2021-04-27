function metrics = defineConvergenceExperimentsMetrics()
metrics = {...
    'Convergence', @convergence; % 0-1 whether the bid process converged
    'NumIterations', @numIterations; % Number of iterations to convergence
    'Welfare', @welfareMetric; % Welfare / sum of all utilities for bidding
    'MeanConsumption', @consumption; % Consumption per user per time step
    'ConsumptionDiff', @consumptionDiff; % Consumption difference from centralized solution
    'WelfareDiff', @welfareDiff; % Absolute welfare difference from centralized solution
    'WelfareDiffPercent' @welfareDiffPercent;  % Percentage welfare difference from centralized solution
    'PriceDiff', @priceDiff; % Price difference from centralized solution
    };
end

function x = convergence(results,expParams,trial)
    if any(isnan(results.bidding.numIterations))
        x = false;
    else
        x = max(results.bidding.numIterations) < expParams.maxIterations;
    end
end

function x = numIterations(results,expParams,trial)
% Number of iterations if converged
    x = max(results.bidding.numIterations);
end

function x = welfareMetric(results,expParams,trial)
    u = trial.confoundingVariables.utilityFunctions;
    d = results.bidding.d;
    x = welfare(u,d);
end

function x = consumption(results,expParams,trial)
    x = mean(results.bidding.d(:));
end

function x = consumptionDiff(results,expParams,trial)
    d = results.bidding.d-results.centralized.d;
    x = norm(d(:))/norm(results.centralized.d(:));
end

function x = welfareDiff(results,expParams,trial)
    u = trial.confoundingVariables.utilityFunctions;
    w1 = welfare(u,results.bidding.d);
    w2 = welfare(u,results.centralized.d);
    x = w1-w2;
end

function x = welfareDiffPercent(results,expParams,trial)
    u = trial.confoundingVariables.utilityFunctions;
    w1 = welfare(u,results.bidding.d);
    w2 = welfare(u,results.centralized.d);
    x = (w1-w2)/w2*100;
end

function x = priceDiff(results,expParams,trial)
    p = results.bidding.price - results.centralized.price;
    x = norm(p(:))/norm(results.centralized.price(:));
end