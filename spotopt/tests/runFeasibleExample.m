function [err,answer,pinf,dinf] = runFeasibleExample(solver,test,varargin)
    [pr,objective,answer] = test(varargin{:});

    sol = pr.minimize(solver,objective);

    pinf = sol.primalInfeasible;
    dinf = sol.dualInfeasible;
    
    if pinf || dinf
        err = Inf;
    else
        err = answer-double(sol.eval(objective));
    end
end
