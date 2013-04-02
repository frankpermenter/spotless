function [err,answer,pinf,dinf,derr] = runFeasibleTest(solver,test,varargin)
    [pr,objective,answer] = test(varargin{:});
    
    [dl,dobj] = toDual(pr,objective);
    
    sol = solver.minimize(pr,objective);
    dsol = solver.minimize(dl,-dobj);
    
    pinf = sol.primalInfeasible;
    dinf = sol.dualInfeasible;
    
    if pinf || dinf
        err = Inf;
    else
        err = [answer-double(sol.eval(objective)) ...
               answer-double(dsol.eval(dobj))];
    end
end
