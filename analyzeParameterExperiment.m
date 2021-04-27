doRun = true; % Set to true to run the experiment; leave false to just load saved results
caseName = '';
% Either run the experiment or load saved results
if doRun
    expObj = runExperiment('parameter');
else
    expObj = ParameterExperiment(caseName);
end

fprintf('Experiment loaded. %i of %i trials completed\n',...
    length(expObj.CompletedTrials), expObj.NTrials);

%%
% Compute the distribution across trials of the convergence and iterations
r = expObj.computeResultsDistributionMetricMultiple({'Convergence','NumIterations'});

%%
% Check overall convergence
convergenceRate = mean(r.Convergence(:));
fprintf('Convergence rate: %g\n',convergenceRate)

if any(~r.Convergence(:))
    disp('Trials where the algorithm didn''t converge for at least one independent variable:');
    disp(find(any(~r.Convergence,2)));
end

%%
% Build array of independent variables. First row is gamma, second is
% delta0
ivStructs = expObj.CompletedTrials(1).treatmentVariables;
ivs = nan(2,length(ivStructs));
for i = 1:length(ivStructs)
    ivs(1,i) = ivStructs(i).gamma;
    ivs(2,i) = ivStructs(i).delta0;
end

%%
% Build array for each trial and treatment variable. Columns are: [gamma,
% delta0, numIterations]
d = nan(numel(r.NumIterations),3);
k = size(r.NumIterations,1);
for i = 1:length(ivStructs)
    d((i-1)*k+1:i*k,1) = ivStructs(i).gamma;
    d((i-1)*k+1:i*k,2) = ivStructs(i).delta0;
    d((i-1)*k+1:i*k,3) = r.NumIterations(:,i);
end

%%
% Plot scatter number iterations against delta0, gamma
figure;
scatter3(d(:,1),d(:,2),d(:,3));
xlabel('gamma');
ylabel('delta0');

%%
% Plot surfaces for mean and max number iterations against delta0, gamma
% Build meshgrid
[Xq, Yq] = meshgrid(linspace(min(ivs(1,:)),max(ivs(1,:))),...
    linspace(min(ivs(2,:)),max(ivs(2,:))));
figure;
surf(Xq,Yq,griddata(ivs(1,:),ivs(2,:),mean(r.NumIterations,1),Xq,Yq),...
    'facecolor','r','facealpha',0.5);
hold on
surf(Xq,Yq,griddata(ivs(1,:),ivs(2,:),max(r.NumIterations,[],1),Xq,Yq),...
    'facecolor','b','facealpha',0.5);
xlabel('\gamma');
ylabel('\delta_0');
title('Number of iterations to convergence');
legend('Mean','Max');

%%
% Make 2-d plots
gammaVals = unique(ivs(1,:));
deltaVals = unique(ivs(2,:));
numBins = 4;

% w.r.t gamma
meanY = nan(length(gammaVals),numBins);
maxY = nan(size(meanY));
labels = {};
for i = 1:length(gammaVals)
    ind1 = ivs(1,:) == gammaVals(i);
    for j=1:numBins
        lb = deltaVals(1) + (deltaVals(end)-deltaVals(1))/numBins*(j-1);
        ub = deltaVals(1) + (deltaVals(end)-deltaVals(1))/numBins*j;
        if j < numBins
            ind2 = ivs(2,:) >= lb & ivs(2,:) < ub;
            labels{j} = sprintf('$%g \\leq \\delta_0 < %g$',lb,ub);
        else
            ind2 = ivs(2,:) >= lb & ivs(2,:) <= ub;
            labels{j} = sprintf('$%g \\leq \\delta_0 \\leq %g$',lb,ub);
        end
        meanY(i,j) = mean(mean(r.NumIterations(:,ind1 & ind2)));
        maxY(i,j) = max(max(r.NumIterations(:,ind1 & ind2)));
    end
end

f = figure;
set(f,'Position',get(f,'Position').*[1 1 1.1 0.7]);
setPageSize(f);
subplot(1,2,1);
semilogy(gammaVals,meanY,'LineWidth',1);
legend(labels,'Interpreter','latex','orientation','vertical','Location','best');
set(gca,'ColorOrderIndex',1)
hold on
semilogy(gammaVals,maxY,'--','LineWidth',0.75,'HandleVisibility','off')
%ylabel('No. Iterations')
%title('Mean and Max Iterations');
xlabel('\gamma');
shrinkWhitespace(gca);

% w.r.t delta_0
meanY = nan(length(deltaVals),numBins);
maxY = nan(size(meanY));
labels = {};
for i = 1:length(deltaVals)
    ind1 = ivs(2,:) == deltaVals(i);
    for j=1:numBins
        lb = gammaVals(1) + (gammaVals(end)-gammaVals(1))/numBins*(j-1);
        ub = gammaVals(1) + (gammaVals(end)-gammaVals(1))/numBins*j;
        if j < numBins
            ind2 = ivs(1,:) >= lb & ivs(1,:) < ub;
            labels{j} = sprintf('$%g \\leq \\gamma < %g$',lb,ub);
        else
            ind2 = ivs(1,:) >= lb & ivs(1,:) <= ub;
            labels{j} = sprintf('$%g \\leq \\gamma \\leq %g$',lb,ub);
        end
        meanY(i,j) = mean(mean(r.NumIterations(:,ind1 & ind2)));
        maxY(i,j) = max(max(r.NumIterations(:,ind1 & ind2)));
    end
end

subplot(1,2,2);
semilogy(deltaVals,meanY,'LineWidth',1);
legend(labels,'Interpreter','latex','orientation','vertical','Location','best');
set(gca,'ColorOrderIndex',1)
hold on
semilogy(deltaVals,maxY,'--','LineWidth',0.75,'HandleVisibility','off')
%ylabel('No. Iterations')
xlabel('\delta_0 (kWh)');
shrinkWhitespace(gca);


% Save in eps and pdf formats
%print -dpdf -painters figures/parameter
%saveas(f,'figures/parameter','epsc');
%saveas(f,'figures/parameter','png');




%% Helper functions
function shrinkWhitespace(ax)
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1) + -0.04;
bottom = ax.Position(2)+0.05;
ax_width = outerpos(3) - ti(1) - ti(3) + 0.04;
ax_height = ax.Position(4);
ax.Position = [left bottom ax_width ax_height];
end

function setPageSize(f)
u1 = get(f,'Units');
set(f,'Units','inches');
screenposition = get(f,'Position');
set(f,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',screenposition(3:4));
set(f,'Units',u1);
end