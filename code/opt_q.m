function [q, alpha, d, p_b, s] = opt_q(u,price,p_s,delta_T,P_b_max,S_max,s0,q_prev,delta,J0)
% Solve the optimization problem for the q player
% Arguments:
%   u: utility function handle, vectorized
%   price: Tx1 vector of prices, pi(t) is price at time t
%   p_s: Tx1 vector of own solar production
%   delta_T: scalar or Tx1 vector specifying the time interval in units
%   such that energy units = delta_T units * power units. A vector can be
%   given for a non-fixed interval
%   P_b_max: scalar or Tx1 vector of battery max charge/discharge rate
%   (max charge and discharge rate assumed to be equal)
%   S_max: scalar maximum storage capacity
%   s0: Initial stored energy in battery
%   q_prev: optional bidding parameter for step-limiting; previous quantity
%   delta: optional bidding paramter for step-limiting; max step length
%   J0: optional scalar, utility from no trade (passed in precomputed to avoid redundant computation)
% Returns:
%   q: Tx1 vector of trade quantities; this player receives q if positive
%   alpha: boolean, true if this quantity and price is better than or equal to J0
%   d: Tx1 vector of consumption
%   p_b: Tx1 vector of battery charge
%   s: Tx1 vector of energy stored
%
% All power variables are assumed to be average power over the interval and
% in the same units.

T = length(price);

% Include step-limiting constraint if information provided
if nargin < 8
    stepLimit = false;
else
    stepLimit = true;
end

% Variable deciding whether or not to solve a feasibility LP for initial
% point
findFeasibleInitialPoint = false;

% Decision variable x:
%   x(1:T): consumption at time t
%   x(T+1:2T): quantity at time t
% If has battery:
%   x(2T+1:3T): battery power at time t
%   x(3T+1:4T): battery energy at time t

% Linear inequality constraints A: none
% Linear equality constraints Aeq:
%   A(1:T,:) -> power balance at time t
%   A(T+1:2*T) -> conservation of energy at time t

dInd = 1:T;
qInd = T+1:2*T;
hasBattery = S_max > 0 && any(P_b_max) > 0;
if hasBattery
    xDim = 4*T;
    AeqDim = 2*T;
    p_bInd = 2*T+1:3*T;
    sInd = 3*T+1:4*T;
else
    xDim = 2*T;
    AeqDim = T;
end

Aeq = zeros(AeqDim,xDim);
beq = zeros(AeqDim,1);
lb = -inf(xDim,1);
ub = inf(xDim,1);
    
% Objective. It is the utility from consumption minus the price*quantity.
% Take negative for minimization.
f = @(x) -sum(u(x(dInd))) + price'*x(qInd);

% Power balance constraint
conInd = 1:T;
Aeq(conInd,dInd) = eye(T);
Aeq(conInd,qInd) = -eye(T);
if hasBattery
    Aeq(conInd,p_bInd) = -eye(T);
end
beq(conInd) = p_s;

% Energy conservation constraint
conInd = T+1:2*T;
if hasBattery
    Aeq(conInd,sInd) = eye(T)-diag(ones(T-1,1),-1);
    Aeq(conInd,p_bInd) = diag(delta_T.*ones(T,1));
    beq(conInd(1)) = s0;
end

% Set bounds:
% Demand must be positive
lb(dInd) = 0;
% Step limiting constraint, formulated as bounds
if stepLimit
    lb(qInd) = q_prev - delta;
    ub(qInd) = q_prev + delta;
end
if hasBattery
    % Battery charge limits
    lb(p_bInd) = -P_b_max;
    ub(p_bInd) = P_b_max;
    % Storage limits
    lb(sInd) = 0;
    ub(sInd) = S_max;
end

% Set up initial point
x0 = zeros(xDim,1);
x0(qInd) = q_prev;
if (findFeasibleInitialPoint)
    % To help with convergence, find a feasible point for initial guess by
    % taking q to be q_prev, and then solving a feasibility LP
    qIndLogical = false(xDim,1);
    qIndLogical(qInd) = true;
    x0(~qIndLogical) = linprog(zeros(sum(~qIndLogical),1),[],[],...
        Aeq(:,~qIndLogical),beq-Aeq(:,qIndLogical)*q_prev,...
        lb(~qIndLogical),ub(~qIndLogical),...
        optimset('Display','none'));
else
    x0(dInd) = max(p_s+x0(qInd),0);
    if hasBattery
        x0(sInd) = s0;
    end
end

opts = optimoptions('fmincon','Algorithm','sqp','Display','none',...
    'MaxFunctionEvaluations', 1e8, 'MaxIterations', 1e4);
% Solve problem
try
    [x,~,exitFlag] = fmincon(f,x0,[],[],Aeq,beq,lb,ub,[],opts);
    %disp(exitFlag)
catch err
    try
        opts = optimoptions('fmincon','Algorithm','interior-point','Display','none',...
            'MaxFunctionEvaluations', 1e8, 'MaxIterations', 1e4); %Try interior method if sqp fails
        [x,~,exitFlag] = fmincon(f,x0,[],[],Aeq,beq,lb,ub,[],opts);
        %disp(exitFlag)
    catch err2
        disp(err2)
    end
    disp(err)
end

if (exitFlag <= 0)
    % Try again starting with feasible initial point by solving an LP
    % holding q_prev fixed
    x0 = zeros(xDim,1);
    x0(qInd) = q_prev;
    qIndLogical = false(xDim,1); % 96x1
    qIndLogical(qInd) = true; % true between 25:48
    optionslin = optimoptions('linprog','Algorithm','interior-point', 'MaxIterations', 1e4);
    [aux, fvalaux, flagaux] = linprog(zeros(sum(~qIndLogical),1),[],[],...
        Aeq(:,~qIndLogical),beq-Aeq(:,qIndLogical)*q_prev,...
        lb(~qIndLogical),ub(~qIndLogical), optionslin);%,...
        %optimset('Display','none'));
    x0(~qIndLogical) = aux;
    
    [x,~,exitFlag] = fmincon(f,x0,[],[],Aeq,beq,lb,ub,[],opts);
    
    if (exitFlag <=0)
        % Try again with lower step tolerance
        opts.StepTolerance = 1e-12;
        [x,~,exitFlag] = fmincon(f,x0,[],[],Aeq,beq,lb,ub,[],opts);
        if (exitFlag <= 0)
            % Now, give up, re-run to display output
            opts.Display = 'iter-detailed';
            fmincon(f,x0,[],[],Aeq,beq,lb,ub,[],opts);
            error('Optimization did not converge. Flag: %g',exitFlag);
        end
    end
end

[~,~,~,d0] = opt_pi(u,-q_prev,zeros(T,1),p_s,delta_T,P_b_max,S_max,s0,zeros(T,1));
J = sum(u(d0)) - price'*q_prev;
if isnan(J0)
    alpha = nan;
else
    alpha = J - J0 >= -1e-12;
end
% Save output variables
q = x(qInd);
d = x(dInd);
if hasBattery
    p_b = x(p_bInd);
    s = x(sInd);
else
    p_b = zeros(T,1);
    s = zeros(T,1);
end

end

