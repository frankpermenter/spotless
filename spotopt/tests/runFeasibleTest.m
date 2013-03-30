function [err,answer,pinf,dinf] = runFeasibleTest(solver,test,varargin)
    [pr,objective,answer] = test(varargin{:});

    sol = solver.minimize(pr,objective);
    
    pinf = sol.primalInfeasible;
    dinf = sol.dualInfeasible;
    
    if pinf || dinf
        err = Inf;
    else
        err = answer-double(sol.eval(objective));
    end
end
