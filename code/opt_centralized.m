function [price, q, d, p_b, s, lambda_b, lambda_c] = opt_centralized(u,p_s,delta_T,P_b_max,S_max,s0,marg_util)
    % Solve the centralized optimization problem for all players
    % Arguments:
    %   u: cell array of length N that contains utility functions
    %   p_s: TxN array of each player's own solar production
    %   delta_T: scalar or Tx1 vector specifying the time interval in units
    %   such that energy units = delta_T units * power units. A vector can be
    %   given for a non-fixed interval
    %   P_b_max: 1xN or TxN array of battery max charge/discharge rate for
    %   player n at time t (max charge and discharge rate assumed to be equal)
    %   S_max: 1xN array of maximum storage capacity for each player N
    %   s0: 1xN array of initial stored energy in each battery
    %   marg_util: (optional) cell array of length N that contains marginal
    %       utility functions
    % Returns:
    %   price: Tx1 vector of prices
    %   q: TxN vector of trade quantities; this player gives q if positive
    %   d: Tx1 vector of consumption
    %   p_b: Tx1 vector of battery charge
    %   s: Tx1 vector of energy stored
    %   lambda_b: TxN vector of optimal dual variables (upper - lower) on the 
    %       battery power constraint
    %   lambda_c: TxN vector of optimal dual variables (upper - lower) on the
    %       battery energy constraint
    %
    % All power variables are assumed to be average power over the interval and
    % in the same units.
    %
    % Note, p_s can be modified to not be strictly solar (renewable)
    % production, but also committed quantities in other trades

    [T, N] = size(p_s);
    NT = N*T;

    if nargin < 8 || ~iscell(marg_util)
        marg_util = {};
    end

    % Decision variable x:
    %   x(1:N*T): consumption by each player at time t
    %       x(1:T) consumption by player 1 at time t
    %       x(T+1:2*T) consumption by player 2 at time t
    %       ...
    %   x(N*T+1:2*N*T): quantity imported by player n at time t
    %   x(2*N*T+1,3*N*T): battery power at time t
    %   x(3*N*T+1,4*N*T): battery power at time t
    % We can optionally drop the "s" variable for storage and compute it
    % directly from p_b
    sub_s = true;
    sub_q = true;

    xDim = 0;
    i = 1;
    inc = NT;
    dInd = i:i+inc-1;
    xDim = xDim + length(dInd);
    i = i+inc;
    if ~sub_q
        inc = NT;
        qInd = i:i+inc+1;
        xDim = xDim + length(qInd);
        i = i+inc;
    end
    p_bInd = i:i+inc-1;
    xDim = xDim + length(p_bInd);
    i = i+inc;
    if ~sub_s
        sInd = i:i+inc-1;
        xDim = xDim + length(sInd);
    end

    AeqDim = T + (~sub_s)*NT + (~sub_q)*NT;

    Aeq = zeros(AeqDim,xDim);
    beq = zeros(AeqDim,1);
    lb = -inf(xDim,1);
    ub = inf(xDim,1);
    nonlcon = [];
    A = [];
    b = [];
        
    % Objective
    if ~isempty(marg_util)
        f = @(x) neg_welfare_with_grad(x,dInd,u,T,N,marg_util);
    else
        f = @(x) -welfare(u,reshape(x(dInd),T,N));
    end

    % Build constraint matrix
    % Keep track of constraint index
    i = 1;

    % Sum of injections (quantities exchanged) is zero
    inc = T;
    Ai = zeros(inc,xDim);
    if sub_q
        % Sum of all d equals sum of battery power + sum of solar at each time
        Ai(:,dInd) = kron(ones(1,N),eye(T));
        Ai(:,p_bInd) = kron(ones(1,N),-eye(T));
        bi = sum(p_s,2);
    else
        % Sum of all injections is zero
        Ai(:,qInd) = kron(ones(1,N),eye(T));
        bi = zeros(inc,1);
    end
    Aeq(i:i+inc-1,:) = Ai;
    beq(i:i+inc-1,:) = bi;
    i = i+inc;

    % Power balance constraint for each player
    if ~sub_q
        % If we are substituting q, then this constraint relating d = q + p_b +
        % p_s is not relevant
        inc = NT;
        Ai = zeros(inc,xDim);
        Ai(:,dInd) = eye(NT);
        Ai(:,qInd) = -eye(NT);
        Ai(:,p_bInd) = -eye(NT);
        bi = p_s(:);

        Aeq(i:i+inc-1,:) = Ai;
        beq(i:i+inc-1,:) = bi;
        i = i+inc;
    end

    % Energy conservation constraint. Can be written one of two ways
    if ~sub_s
        % Write an equality constraint relating s_t, s_t+1 and p_b_t
        % We will write the 0 <= s_t <= S_max as a bound
        inc = NT;
        Ai = zeros(inc,xDim);
        Ai(:,sInd) = kron(eye(N),eye(T)-diag(ones(T-1,1),-1));
        Ai(:,p_bInd) = kron(eye(N),diag(delta_T.*ones(T,1)));
        bi = zeros(NT,1);
        bi(1:T:end) = s0;

        Aeq(i:i+inc-1,:) = Ai;
        beq(i:i+inc-1,:) = bi;
        i = i+inc;
    else
        % We are substituting for s, so write the inequality constraint
        % 0 <= s_0 - A*p_b <= S_max where A is a lower triangular matrix for
        % each agent N
        A = zeros(2*NT,xDim);
        b = zeros(2*NT,1);
        lambda_c_upper_ind = 1:NT;
        lambda_c_lower_ind = NT+1:2*NT;
        A(lambda_c_upper_ind,p_bInd) = -kron(eye(N),tril(ones(T)));
        A(lambda_c_lower_ind,p_bInd) = kron(eye(N),tril(ones(T)));
        b(lambda_c_upper_ind) = kron(S_max-s0,ones(T,1));
        b(lambda_c_lower_ind) = kron(s0,ones(T,1));
    end

    % Set bounds:
    % Demand must be positive
    lb(dInd) = 0;
    % Battery charge limits
    if size(P_b_max,1) == 1
        P_b_max = repmat(P_b_max,T,1);
    end
    lb(p_bInd) = -P_b_max(:);
    ub(p_bInd) = P_b_max(:);
    % Storage limits
    if ~sub_s
        ub(sInd) = kron(S_max',ones(T,1));
    end

    x0 = zeros(xDim,1);
    try
        opts = optimoptions('fmincon','Algorithm','interior-point','Display','iter', 'MaxFunctionEvaluations', 1e8,'StepTolerance',1e-20); % sqp to get better accuracy on lagrange multipliers
        opts.OptimalityTolerance = 1e-5;
        opts.MaxIterations = 2000;
        if ~isempty(marg_util)
            opts.SpecifyConstraintGradient = true;
            opts.CheckGradients = true;
        end
        [x,~,exitflag,~,lambda] = fmincon(f,x0,sparse(A),b,sparse(Aeq),beq,lb,ub,nonlcon,opts); %sqp

        if (exitflag <= 0)
            error('Optimization did not converge. Flag: %g',exitflag);
        end
        if (exitflag > 1)
            warning('Optimization may not have converged to optimal point. Flag: %g', exitflag);
        end
    catch err
        error(err);
    end

    d = reshape(x(dInd),T,N);
    p_b = reshape(x(p_bInd),T,N);
    if sub_q
        q = d - p_b - p_s;
    else
        q = reshape(x(qInd),T,N);
    end
    if sub_s
        s = A(1:NT,p_bInd)*x(p_bInd)+b(NT+1:2*NT);
    else
        s = x(sInd);
    end
    s = reshape(s,T,N);

    price = abs(lambda.eqlin(1:T));
    lambda_b = reshape(lambda.upper(p_bInd) - lambda.lower(p_bInd),T,N);
    if sub_s
        lambda_c = lambda.ineqlin(lambda_c_upper_ind) - lambda.ineqlin(lambda_c_lower_ind);
    else
        lambda_c = lambda.upper(sInd) - lambda.lower(sInd);
    end
    lambda_c = reshape(lambda_c,T,N);

end

function [fval, grad_val] = neg_welfare_with_grad(x,dInd,u,T,N,marg_util)
    [fval, ~, marg_util_vals] = welfare(u,reshape(x(dInd),T,N),marg_util);
    fval = -fval;
    grad_val = zeros(size(x));
    grad_val(dInd) = -marg_util_vals(:);
end