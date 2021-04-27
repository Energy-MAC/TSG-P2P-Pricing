function [u, f, g] = fit_utility_N(elasticity,empiricalPrice,empiricalQuantity)
% Fit utility functions for N players over T time steps. The fit computes
% NT utility functions by centering a constant elasticity function around
% each empirical price / quantity pair using the given elasticity for each
% player. The return values are 1xN cell arrays of function handles. The 
% first return value, u, is such that each function handle u(n) accepts a
% TxM array 'd' of consumption values, where d(t,m) is some consumption at
% time t. u{n}(d) returns a TxM array of utilities. The additional
% arguments f and g behave the same.
%
% Arguments:
%   elasticity: vector of length N, elasticity for each player
%   empiricalPrice: NxT array of empirical prices
%   empiricalQuantities: NxT array of empirical quantities

% This could be an argument, used to offset constant elasticity
d0 = 1e-2;
[~, N] = size(empiricalPrice);
u = cell(1,N); % Utility function
f = cell(1,N); % Maringal utility function
g = cell(1,N); % Inverse marginal utility / demand function
for n = 1:N
    [u{n}, f{n}, g{n}] = fit_utility(elasticity(n), empiricalPrice(:,n),empiricalQuantity(:,n), d0);
end
end

