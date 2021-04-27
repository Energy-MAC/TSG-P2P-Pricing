function [y, u_vals, marg_u_vals] = welfare(u,d,marg_util)
% Compute the welfare given consumption d and utilities u. Assumes u is in
% the format that u{n} is a vectorized utility function such that u{n}(d)
% returns a vector of the same size d.
    if nargin < 3 || ~iscell(marg_util)
        marg_util = {};
    end
    
    [T,N] = size(d);
    u_vals = nan(T,N);
    for n = 1:N
        u_vals(:,n) = u{n}(d(:,n));
    end
    y = sum(u_vals(:));
    
    if nargout > 2
        if isempty(marg_util) && T*N > 1
            error('Must provide marginal utility functions if marginal value is requested');
        else
            marg_u_vals = nan(T,N);
            for n = 1:N
                marg_u_vals(:,n) = marg_util{n}(d(:,n));
            end
        end
    end
end