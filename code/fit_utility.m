function [u, f, g] = fit_utility(elasticity,price,consumption, d0)
% Returns a quasi-constant elasticity utility function to an empirical price and a consumption
% pair given an elasticity. The returned utility function will have
% elasticity equal to 'elasticity' at consumption equal to 'consumption'
% and marginal utility 'f' equal to 'price' at consumption, i.e.
% f(consumption) = price.

% Arguments
%   elasticity: scalar, elasticity value
%   price: Tx1 vector of observed prices
%   consumption: Tx1 vector of observed quantities corresponding to price
%   d0: strictly positive parameter to ensure that marginal utility is
%   finite at 0.

% d0 > 0 is used so that the marginal utility is finite at 0.
% If d0 = 0 with constant elasticity then the utility function is negative
% infinity at zero, u(0) -> -Inf, so d0 must be strictly positive

% The elasticity is given as a function of d by elasticity(d) = r(1+d0/d)

% The closer d0 is to zero, the closer to constant elasticity the demand
% function is, and the larger the marginal utility will be at 0.

if nargin < 4
    d0 = 1e-2;
end
if d0 <= 0
    error('d0 must be strictly greater than zero');
end

if any(consumption <= 0) || any(price <= 0)
    error('Consumption and price must be strictly greater than zero for fit');
end

[T,x] = size(price);
[T2,x2] = size(consumption);
if ~(x == 1 && x2 == 1 && T == T2)
    error('Both price and consumption must be Tx1 arrays');
end

% Compute r so that the elasticity at 'consumption' is equal to
% 'elasticity'
r = elasticity./(1+d0./consumption);

% Pre-compute some constants for faster evaluation of the utility function
a = r./(r+1)./(consumption+d0).^(1./r).*price;
b = (1./r+1);
c = d0.^(1./r+1);

% Utility function
u = @(d) a.*((d+d0).^b-c);
if nargout > 1
    % Marginal utility
    f = @(d) ((d+d0)./(consumption+d0)).^(1./r).*price;
end
if nargout > 2
    % Demand function
    g = @(pi) (pi./price).^r.*consumption;
end
end