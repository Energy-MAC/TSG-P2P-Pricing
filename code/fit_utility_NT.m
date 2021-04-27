function [u, f, g] = fit_utility_NT(elasticity,empiricalPrice,empiricalQuantity)
% Fit utility functions for N players over T time steps. The fit computes
% NT utility functions by centering a constant elasticity function around
% each empirical price / quantity pair using the given elasticity for each
% player.
%
% WARNING: Better to use fit_utility_N for performance, which vectorizes
% the utility functions over T
%
% Arguments:
%   elasticity: vector of length N, elasticity for each player
%   empiricalPrice: NxT array of empirical prices
%   empiricalQuantities: NxT array of empirical quantities

% This could be an argument, used to offset constant elasticity
d0 = 1e-2;
[T, N] = size(empiricalPrice);
u = cell(T,N); % Utility function
f = cell(T,N); % Maringal utility function
g = cell(T,N); % Inverse marginal utility / demand function
for n = 1:N
    for t = 1:T
        [u{t,n}, f{t,n}, g{t,n}] = fit_utility(elasticity(n), empiricalPrice(t,n),empiricalQuantity(t), d0);
    end
end
end

