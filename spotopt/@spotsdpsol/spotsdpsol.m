classdef spotsdpsol
    properties 
        program = [];
        user_variables = [];
        primalSolution = [];
        G = [];
        h = [];
        info = [];
    end
    
    methods
        function f = primalInfeasible(sol)
            f = sol.info.pinf;
        end
        function f = dualInfeasible(sol)
            f = sol.info.dinf;
        end
        
        function sol = spotsdpsol(prog,info,uservars,psol,G,h)
            if nargin < 5, G = speye(length(psol)); end
            if nargin < 6, h = sparse(length(psol),1); end

            sol.program = prog;
            sol.info = info;
            sol.user_variables = uservars;
            sol.primalSolution = psol;
            sol.G = G;
            sol.h = h;

        end
        
        function ev = eval(sol,exp)
            if sol.primalInfeasible
                error('Primal Infeasible, cannot evalutate primal variables.');
            end

            ev = subs(exp,[sol.user_variables],...
                      [sol.G*sol.primalSolution+sol.h]);
        end
    end
end