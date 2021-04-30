doRun = true; % Set to true to run the experiment; leave false to just load saved results
% Either run the experiment or load saved results
if doRun
    expObj = runExperiment('storage','expC');
else
    expObj = StorageExperiment('expC');
end

% Compute the distribution across trials of each metric
r = expObj.computeResultsDistributionMetricMultiple();

% Reference these distributions as a struct fields. The row index is the
% trial* and the column index is the treatment/independent variable set
% For example, view the convergence by trial and treatment:
disp(r.Convergence);
disp('Trials where the algorithm didn''t converge for at least one independent variable:')
disp(find(any(~r.Convergence,2))); %Full convergence

%% Convergence Performance

iters = reshape(r.NumIterations, [60*20, 1]);
mean_iterations = mean(iters)
std_iterations = std(iters)
median_iterations = median(iters)


%% Battery parameters
num_iter = r.NumIterations;
size15 = [1,6,11,16];
size25 = size15 + 1;
size40 = size25 + 1;
size80 = size40 + 1;
size300 = size80 + 1;

data15 = reshape(num_iter(:, size15), numel(num_iter(:, size15)), 1);
data15(isnan(data15)) = [];

data25 = reshape(num_iter(:, size25), numel(num_iter(:, size25)), 1);
data25(isnan(data25)) = [];

data40 = reshape(num_iter(:, size40), numel(num_iter(:, size40)), 1);
data40(isnan(data40)) = [];

data80 = reshape(num_iter(:, size80), numel(num_iter(:, size80)), 1);
data80(isnan(data80)) = [];

data300 = reshape(num_iter(:, size300), numel(num_iter(:, size300)), 1);
data300(isnan(data300)) = [];

x = [data15;data25;data40;data80; data300];
g = [ones(size(data15)); 2*ones(size(data25)); 3*ones(size(data40)); 4*ones(size(data80)); 5*ones(size(data300))];
fig = figure();
boxplot(x,g, 'Symbol','')
xticklabels({'15','25','40','80','300'})
ylim([0, 450])
xlabel('Total Battery Capacity [kWh]')
ylabel('# Iterations')

%% Welfare Comparisons
T_24 = [];
T_1 = [];
T_12 = [];

results = expObj.CompletedTrials;
for i=1:60
    aux = size(results(i).confoundingVariables.solarShape);
    T = aux(1);
    if T == 24
        T_24(end+1) = i;
    elseif T == 12
        T_12(end+1) = i;
    else
        T_1(end+1) = i;
    end
end

welf_diff = r.WelfareDiff; 
%welf_diff = r.WelfareDiffPercent; %Uncomment for percentage differences

data1 = reshape(welf_diff(T_1, :), numel(welf_diff(T_1, :)), 1);
data1(isnan(data1)) = [];
data1 = abs(data1);
data12 = reshape(welf_diff(T_12, :), numel(welf_diff(T_12, :)), 1);
data12(isnan(data12)) = [];
data12 = abs(data12);
data24 = reshape(welf_diff(T_24, :), numel(welf_diff(T_24, :)), 1);
data24(isnan(data24)) = [];
data24 = abs(data24);


% Table Results
mean_w1 = mean(data1)
std_w1 = std(data1)
max_w1 = max(data1)

mean_w12 = mean(data12)
std_w12 = std(data12)
max_w12 = max(data12)

mean_w24 = mean(data24)
std_w24 = std(data24)
max_w24 = max(data24)

%% Analysis for special instance (Table II)

centr_sol = expObj.CompletedTrials(47).outputs(6).centralized;
bidd_sol = expObj.CompletedTrials(47).outputs(6).bidding;
u = expObj.CompletedTrials(47).confoundingVariables.utilityFunctions;
[y, uIndsCentr] = welfare(u, centr_sol.d);
[y, uIndsBid] = welfare(u, bidd_sol.d);

qInd = [1,2,4,5,6];

% Centralized Solution
centr_welf = uIndsCentr - centr_sol.price .* centr_sol.q;
centr_welf_qPlayers = centr_welf(:, qInd);
centr_welf_piPlayer = centr_welf(:, 3);
total_centr_welf_qPlayers = sum(centr_welf_qPlayers, 1);
total_centr_welf_piPlayer = sum(centr_welf_piPlayer);

% Bidding Solution
price_trade_qPlayers = bidd_sol.price .* bidd_sol.q(:, qInd);
bidd_welf_qPlayers = uIndsBid(:, qInd) - price_trade_qPlayers;
bidd_welf_piPlayer = uIndsBid(:, 3) + sum(price_trade_qPlayers, 2);
total_bidd_welf_qPlayers = sum(bidd_welf_qPlayers, 1);
total_bidd_welf_piPlayer = sum(bidd_welf_piPlayer);

% Totals
welf_centr = sum(total_centr_welf_qPlayers) + total_centr_welf_piPlayer;
welf_bidd =  sum(total_bidd_welf_qPlayers) + total_bidd_welf_piPlayer;

% Differences
Deltadiff = welf_centr - welf_bidd
Deltapercentage = Deltadiff/welf_centr * 100

% Plot
fig2 = figure();
hold on
plot(bidd_sol.price, 'Linewidth', 1.5)
plot(centr_sol.price, 'Linewidth', 1.5)
legend({'\pi_1', '\pi_2', '\pi_3', '\pi_4', '\pi_5', '\pi_{centr}'}, 'Location', 'northwest')
xlabel('Hour of the Day')
ylabel('Price [$/kWh]')
xlim([1,24])
ylim([0,16])
