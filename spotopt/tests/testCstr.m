function [] = testCstr(prorig,solver,logger)
%
%  This is an automated testing script.
%
%  [] = testCstr(pr,solver)
%
%  Tests all binary combinations of constraints in the 
%
    
    if nargin < 3, logger = @disp; end
    
    listing = what('testCstr');
    mFiles = listing.m;
    
    N = length(mFiles);
    
    for i = 1:N
        fns{i} = str2func(mFiles{i}(1:length(mFiles{i})-2));
    end
    
    obj  = msspoly(zeros(N,1));
    zero = msspoly(zeros(N,1));
    tol  = zeros(N,1);
    
    pr = prorig;
    for i = 1:N
        [pr,obj(i),zero(i),tol(i)] = fns{i}(pr);
    end
    
    sol = pr.minimize(solver,sum(obj));
    zero = double(sol.eval(zero));
    J = find(tol <= abs(double(sol.eval(zero))));
    
    if length(J) > 0,
        logger(sprintf('Combined test failed:'));
    end
    for j = 1:length(J)
        logger(sprintf('\tTest %s: combined (%d > %d).',...
                       mFiles{J(j)},...
                       zero(J(j)),tol(J(j))));
    end
    
    for i = 1:N
        pr = prorig;
        [pr,obj,zero,tol] = fns{i}(pr);
        sol = pr.minimize(solver,obj);
        zero = double(sol.eval(zero));
        
        if abs(zero) > tol
            logger(sprintf('\tIndependent test %s failed:  (%d > %d).',...
                           mFiles{J(j)},zero,tol));
        end
    end
end