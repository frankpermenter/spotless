classdef (Abstract) spotsdpsolver 
    methods (Abstract)
        %inst = isInstalled(solvr);
        can = canSolve(solvr,sdp);
        sol = optimize(solvr,sdp,objective,options);
    end
end