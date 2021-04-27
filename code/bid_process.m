function [priceFinal,qFinal,numIterationsQ,price,qOffer,qReq,delta,dFinal,priceFinalQ,pBFinal,sFinal] = ...
    bid_process(u,piInd,p_s,delta_T,P_b_max,S_max,s0,delta0,gamma,epsilon,maxIterations)
% Run the bid process for N players. There is 1 pi-player who sets price
% and feasible quantity, and Nq = N-1 q-players who request quantities
% Arguments:
%   u: cell array vector of length N of utility function handles for each
%   player, each is vectorized w.r.t time
%   piInd: numeric index of the pi-player
%   p_s: TxN array of each player's own solar production
%   delta_T: scalar or Tx1 vector specifying the time interval in units
%   such that energy units = delta_T units * power units. A vector can be
%   given for a non-fixed interval
%   P_b_max: 1xN or TxN array of battery max charge/discharge rate for
%   player n at time t (max charge and discharge rate assumed to be equal)
%   S_max: 1xN array of maximum storage capacity for each player N
%   s0: 1xN array of initial stored energy in each battery
%   delta0: scalar, initial value of bidding paramter for step-limiting; max step length
%   epsilon: scalar, tolerance for convergence
%   maxIterations: scalar, maximum number of iterations until termination
% Returns:
%   priceFinal: TxNq array of agreed price at time t for the n'th q player
%   qFinal: TxNq array of agreed quantity at time t for n'th q player
%   numIterationsQ: Nqx1 vector, number of iterations to convergence for each q player
%   price: IxT array of the time-of-use price at each iteration
%   qFeas: IxTxN array of feasible quantity returned at each iteration by
%   the pi-player; settled quantities are copied forward
%   qReq: Same as qFeas, but quantity requested by q-player
%   delta: IxTxN array of deltas
%
% All power variables are assumed to be average power over the interval and
% in the same units.
%   pi_Ind: index of the pi player

[T, N] = size(p_s);

Nq = N-1; % Only 1 pi-player allowed at this point
qInd = true(N,1);
qInd(piInd) = false;
qIndNum = find(qInd); % This maps the q-player's index in the set of q players to the set of all players

% Initialize and pre-allocated variables
qReq = zeros(1,T,Nq); % Trajectory of quantities requested by q players
qOffer = zeros(1,T,Nq); % Trajectory of offered quantities returned by pi players
qFeas = zeros(1,T,Nq); % Trajectory of quantities known to be feasible
alpha = false(1,N); % Trajectory of alpha (alpha(i,n) = true means that player n would except the offer at iteration i
price = nan(1,T); % Trajectory of prices at each iteration
numIterationsQ = nan(Nq,1); % Number of iterations until termination for each q player
priceFinal = nan(T,Nq); % Final price agreed on for each q player at each time
qFinal = zeros(T,Nq); % Final quantity agreed on for each player at each time. Augments to include pi player at end
dFinal = nan(T,N); % Final consumption agreed on for each player at each time.
delta = delta0*ones(1,T,Nq);
priceFinalQ = nan(T,Nq); % Final shadow price of the q players on agreed upon solution
pBFinal = nan(T,N); % Final battery charge power of each player at each time. 
sFinal = nan(T,N); % Final battery state of charge of each player at each time

% Compute the baseline utility for each player at zero quantity. This is
% used for each player to determine whether a trade is better than doing
% nothing
J0 = zeros(1,N);
for n = 1:length(J0)
    [~,~,~,d0] = opt_pi(...
        u{n},zeros(T,1),zeros(T,1),p_s(:,n),delta_T,P_b_max(:,n),S_max(n),s0(n),...
        zeros(T,1));
    J0(n) = sum(u{n}(d0));
end

minDelta = 1e-7; % This prevents delta from getting so small that it causes
% numerical feasibility issues in solving the q player's problem
osc_tol = 1e-6; % This is a tolerance that relaxes the inequality constraints for oscillation

% Start iteration loop
i = 1;
while (true)
%     % Comment out/in and place breakpoint to intercept the loop while debugging.
%     if (i > 21)
%         disp('Debugging: many iterations');
%     end

    % Copy forward delta and qReq for next iteration; will get overwritten
    % for q players that haven't settled
    delta(i+1,:,:) = delta(i,:,:);
    qReq(i+1,:,:) = qReq(i,:,:);
    qExitedInd = ~isnan(numIterationsQ); % Logical index in set of q players for which have exited
    
    % Run pi-player optimization and save outputs. price, qOffer, and alpha
    % are tracked as trajectories and used during the iteration. The
    % remaining outputs are final values and overwritten each iteration
    [price(i,:), qOffer(i,:,~qExitedInd), alpha(i,piInd), dFinal(:,piInd), ...
        pBFinal(:,piInd),sFinal(:,piInd)] = opt_pi(...
        u{piInd},...% pi-player's utiity function
        reshape(qReq(i,:,~qExitedInd),T,sum(~qExitedInd)),...% Quantity requested by all q-players that haven't exited, shaped appropriately
        qFinal(:,qExitedInd),...% Settled / promised quantity to q-players who have exited
        p_s(:,piInd),delta_T,P_b_max(:,piInd),S_max(piInd),s0(piInd),...% Standard parameters for the pi-player
        reshape(qFeas(i,:,~qExitedInd),T,sum(~qExitedInd)),...% A known feasible point, used in the beta projection
        J0(piInd),...% pi-player's utility for no trade, used to compute alpha
        priceFinal(:,qExitedInd));% price promised to settled q-players, used to compute alpha
    
    % Solve each q-player's problem given price and qOffer to record if
    % they want to exit
    requestExit = false(Nq,1);
    for n = 1:Nq
        if ~isnan(numIterationsQ(n))
            % Player has already settled, so skip
            continue
        end
        ind = qIndNum(n); % n is the index in the set of q players; ind is the index in the set of all players
        qPrev = reshape(qOffer(i,:,n),T,1);
        [qNext, alpha(i,ind)] = opt_q(...
            u{ind},price(i,:)',p_s(:,ind),delta_T,P_b_max(:,ind),...
            S_max(ind),s0(ind),qPrev,reshape(delta(i,:,n),T,1),J0(ind));
        qReq(i+1,:,n) = qNext;
        
        % Check tolerance
        if norm(qNext-qPrev,1) <= max(gamma*epsilon*norm(qPrev,1),T*minDelta+1e-12)
            % Request to exit
            requestExit(n) = true;
            continue; % Do not shrink if an exit is requested
        end
        
        % Check oscillation at each time
        if (i > 1)
            qPrev2 = reshape(qOffer(i-1,:,n),T,1);
        else
            qPrev2 = zeros(T,1);
        end
        oscillating = ~(...
            (qNext - qPrev > -osc_tol & qPrev - qPrev2 > osc_tol) | ...
            (qNext - qPrev < osc_tol & qPrev - qPrev2 < -osc_tol) ...
            );
        % Shrink delta where oscillating
        delta(i+1,oscillating,n) = max(delta(i,oscillating,n)*gamma,minDelta);
    end
    
    % Resolve the exits
    % Line below can break convergence; it's not guaranteed we will hit a
    % point where alpha for the pi-player is true with more than one
    % player. This can happen, e.g. with 2 q-players (QA and QB) and one pi
    % player (P). Say QA buys a positive q at a low price, and P is buying
    % from QB to cover it, and QA exits. If it turns out at this price QB
    % wants to sell less, then QB will bid less and the price goes up. Now
    % P has to buy at a higher price to cover what they promised to QB. It
    % is possible that this loss can outweigh the gain that P gets from
    % other trades. However, the smaller epsilon is, the less likely this
    % is to happen.
    %if all(alpha(i,[piInd;qIndNum(~qExitedInd)])) && any(requestExit)
    % Allow exit if all Q players who have not exited are satisfied
    if all(alpha(i,qIndNum(~qExitedInd))) && any(requestExit)
        % The last offer was a feasible trade, so store at as most recent.
        % This is the point the pi-player uses as the reference in the beta
        % shrinking process to find a feasible point.
        qFeas(i+1,:,:) = qOffer(i,:,:);
        % Allow all q players who requested an exit to exit, because every
        % player would prefer this allocation to no trade
        for n = find(requestExit)'
            % Exit, settle on the quantity and price returned by the pi
            % player
            ind = qIndNum(n);
            qFinal(:,n) = reshape(qOffer(i,:,n),T,1);
            priceFinal(:,n) = price(i,:);
            numIterationsQ(n) = i;
            % Compute own shadow price, consumption, and battery given q to
            % record as outputs. This is equivalent to applying the pi
            % problem to the q player (with the opposite sign of the final
            % q).
            [priceFinalQ(:,n), qTemp, ~, dFinal(:,ind), pBFinal(:,ind),...
                sFinal(:,ind)] = opt_pi(...
                u{ind},-qFinal(:,n),zeros(T,1),p_s(:,ind),delta_T,P_b_max(:,ind),S_max(ind),s0(ind),zeros(T,1));
            % Double check that this solution returns the same q (i.e.
            % there was no beta shrinking on the solution to get
            % feasibility). This should always be the case barring a bug or
            % numerical issue because the final q was an offer from the
            % pi-player. Every offer from the pi-player lies on the line
            % connecting two feasible points, and is therefore convex
            % because the q-player's constraints are convex (linear). Note,
            % the two feasible points are the quanitity requested by the q
            % player, and a known feasible point.
            if any(abs(qFinal(:,n)+qTemp) > 1e-3)
                error('Something is wrong. Appears the final q is infeasible for the q player.');
            end
        end
    else
        qFeas(i+1,:,:) = qFeas(i,:,:);
    end
    
    % Exit the process if all q players have exited or the maximum
    % iterations are hit
    if (~any(isnan(numIterationsQ))) || i >= maxIterations
        fprintf('Exiting bid process after %i iterations\n',i);
        break;
    end
    qOffer(i+1,:,:) = qOffer(i,:,:);
    i = i+1;
end

% Remove extraneous values of q from the trajectory
qReq = qReq(2:end,:,:); % First value is trivially zero
qFinal = [qFinal(:,1:piInd-1) -sum(qFinal,2) qFinal(:,piInd:end)];

end