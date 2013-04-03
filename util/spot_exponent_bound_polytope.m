function [ptope] = spot_exponent_bound_polytope(pow)
%
%   [ptope] = spot_exponent_bound_polytope(pow)
%
%   pow   -- matrix of nonnegative integers.
%   ptope -- matrix of nonnegative integers such that
%            ZZ^n intersect conv(pow(:,1),pow(:,2),...,pow(:,end))
%            belongs to ptope.
%
    options.sos.newton = 1;                                             
    options.sos.csp = 0;
    options.verbose = 1;
    options.sos.inconsistent = 0;
    [C,csclasses] = corrsparsity(pow,options);
    mpow = monomialgeneration(pow,csclasses);

    ptope = mpow{1}; 
end