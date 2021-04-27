function [price, q, alpha, d, p_b, s] = opt_pi(u,q,qPromised,p_s,delta_T,P_b_max,S_max,s0,qFeas,J0,pricePromised)
% Solve the optimization problem for the pi player
% Arguments:
%   u: own utility function, 1xT vectorized
%   q: TxN vector of requested quantity, q(t,n) is the quantity requested by player n
%   at time t
%   qPromised: TxX array of committed quantities (X is number of settled
%   q-players)
%   p_s: Tx1 vector of own solar production
%   delta_T: scalar or Tx1 vector specifying the time interval in units
%   such that energy units = delta_T units * power units. A vector can be
%   given for a non-fixed interval
%   P_b_max: scalar or Tx1 vector of battery max charge/discharge rate
%   (max charge and discharge rate assumed to be equal)
%   S_max: scalar maximum storage capacity
%   s0: Initial stored energy in battery
%   qFeas: TxN array of a known feasible quantity
%   J0: scalar, own welfare at zero q, used for computing alpha. Defaults
%   to nan
%   pricePromised: TxX array of committed prices
% Returns:
%   price: Tx1 vector of prices
%   q: TxN vector of trade quantities; this player GIVES q if positive
%   alpha: boolean, true if this quantity and price is better than J0
%   d: Tx1 vector of consumption
%   p_b: Tx1 vector of battery charge
%   s: Tx1 vector of energy stored
%
% All power variables are assumed to be average power over the interval and
% in the same units.
%
% Note, p_s can be modified to not be strictly solar (renewable)
% production, but also committed quantities in other trades

if nargin < 10
    J0 = nan;
    pricePromised = nan;
end

[T, ~] = size(p_s);

if isempty(qPromised)
    qPromised = zeros(T,1);
end
if isempty(pricePromised)
    pricePromised = zeros(T,1);
end

hasBattery = S_max > 0 && any(P_b_max) > 0;

% First, solve a problem to ensure feasibility: shrink all quantities
% uniformly, minimize the shrinking such that the main problem is feasibly

beta = beta_projection(q,qPromised,p_s,delta_T,P_b_max,S_max,s0,qFeas);
% Project q if necessary back towards feasible point
q = beta*qFeas + q*(1-beta);

% Decision variable x:
%   x(1:T): consumption at time t
% If has battery:
%   x(T+1:2T): battery power at time t
%   x(2T+1:3T): battery energy at time t

% Linear inequality constraints A: none
% Linear equality constraints Aeq:
%   Aeq(1:T,:) -> power balance at time t
%   Aeq(T+1:2*T) -> conservation of energy at time t

dInd = 1:T;
if hasBattery
    xDim = 3*T;
    AeqDim = 2*T;
    p_bInd = T+1:2*T;
    sInd = 2*T+1:3*T;
else
    xDim = T;
    AeqDim = T;
end

Aeq = zeros(AeqDim,xDim);
beq = zeros(AeqDim,1);
lb = -inf(xDim,1);
ub = inf(xDim,1);
nonlcon = []; % Will get set lower if applicable
    
% Objective
f = @(x) -sum(u(x(dInd)));

% Power balance constraint
conInd = 1:T;
Aeq(conInd,dInd) = eye(T);
if hasBattery
    Aeq(conInd,p_bInd) = -eye(T);
end
beq(conInd) = p_s-sum(qPromised,2)-sum(q,2);

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
if hasBattery
    % Battery charge limits
    lb(p_bInd) = -P_b_max;
    ub(p_bInd) = P_b_max;
    % Storage limits
    lb(sInd) = 0;
    ub(sInd) = S_max;
end

% Feasible x0
x0 = zeros(xDim,1);
x0(dInd) = p_s-sum(qPromised,2);
if hasBattery
    x0(sInd) = s0;
end
opts = optimoptions('fmincon','Algorithm','sqp','Display','none',...
    'MaxFunctionEvaluations', 1e8, 'MaxIterations', 1e4);  % sqp to get better accuracy on lagrange multipliers
opts.StepTolerance = opts.ConstraintTolerance/100; % prevent early termination
[x,fval,exitflag,~,lambda] = fmincon(f,x0,[],[],Aeq,beq,lb,ub,nonlcon,opts);

if (exitflag <= 0)
    error('Optimization did not converge. Flag: %g',exitflag);
end

price = abs(lambda.eqlin);
price = price(1:T); %Only takes the power balance constraints
if isnan(J0)
    alpha = nan;
else
    J = -fval + price'*sum(q,2) + pricePromised(:)'*qPromised(:);
    alpha = J - J0 >= -1e-6;
end
d = x(dInd);
if hasBattery
    p_b = x(p_bInd);
    s = x(sInd);
else
    p_b = zeros(T,1);
    s = zeros(T,1);
end
end

function beta = beta_projection(q,qPromised,p_s,delta_T,P_b_max,S_max,s0,qFeas)

T = size(q,1);

hasBattery = S_max > 0 && any(P_b_max) > 0;
% Solve a problem to ensure feasibility: shrink all quantities
% uniformly, minmize the shrinking such that the main problem is feasibly
% Decision variable x:
%   x(1) = beta
% If battery:
%   x(2:T+1) = battery power
%   x(T+2:2*T+1) = battery energy

betaInd = 1;
ADim = T;
if hasBattery
    xDim = 1+2*T;
    AeqDim = T;
    p_bInd = 2:T+1;
    sInd = T+2:2*T+1;
else
    xDim = 1;
    AeqDim = 0;
end

% Linear inequality constraints A:
%   A(1:T,:) -> quantity exported must be less than solar plus power from
%   battery
% Linear equality constraints Aeq:
%   Aeq(1:T,:) -> conservation of energy at time t
A = zeros(ADim,xDim);
b = zeros(ADim,1);
Aeq = zeros(AeqDim,xDim);
beq = zeros(AeqDim,1);
lb = -inf(xDim,1);
ub = inf(xDim,1);

% Objective: minimize beta
f = zeros(xDim,1);
f(1) = 1;
% Linear inequality power balance constraint
conInd = 1:T;
A(conInd,betaInd) = sum(qFeas,2)-sum(q,2);
if hasBattery
    A(conInd,p_bInd) = -eye(T);
end
b(conInd) = p_s-sum(qPromised,2)-sum(q,2);

% Energy conservation constraint
if hasBattery
    conInd = 1:T;
    Aeq(conInd,sInd) = eye(T)-diag(ones(T-1,1),-1);
    Aeq(conInd,p_bInd) = diag(delta_T.*ones(T,1));
    beq(conInd(1)) = s0;
end

% Bounds for beta: 0 <= beta <= 1
lb(1) = 0;
ub(1) = 1;

% Bounds for battery    
if hasBattery
    % Battery charge limits
    lb(p_bInd) = -P_b_max;
    ub(p_bInd) = P_b_max;
    % Storage limits
    lb(sInd) = 0;
    ub(sInd) = S_max;
end
opts = optimoptions('linprog','Display','off',...
    'ConstraintTolerance',1e-9,'OptimalityTolerance',1e-10);
while (opts.ConstraintTolerance <= 1e-6)
    [x,~,exitFlag] = linprog(f,A,b,Aeq,beq,lb,ub,opts);
    if exitFlag > 0
        % Success
        break;
    end
    % Try again with relaxed constraint tolerance
    opts.ConstraintTolerance = opts.ConstraintTolerance*10;
end
if (exitFlag <= 0)
    error('Optimization did not converge. Flag: %g',exitFlag);
end
beta = x(1);

% Force to 0 if close
if abs(beta) < 1e-9
    beta = 0;
end

end
