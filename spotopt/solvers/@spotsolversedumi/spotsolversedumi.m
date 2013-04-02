classdef spotsolversedumi < spotsdpsolver
    properties
        options = struct();
    end
    
    
    methods 
        
        function solvr = spotsolversedumi(options)
            ;
            % TODO: check that all options are legal.
            if nargin > 0
                solvr.options = options;
            else
                solvr.options = struct('fid',0);
            end
        end
        
        
        function can = canSolve(solvr,sdp)
            can = 1;
        end
        

        
        function sol = minimizePrimalForm(solver,pr,objective)
            if nargin < 2, objective = 0; end
            
            [A,b,c,K,G,h,varNo] = spotsdpsolver.sdpToPrimalSedumiFormat(pr,objective);

            [x,y,info] = sedumi(A,b,c,K,solver.options);

            if info.pinf, 
                primalSol = NaN*ones(size(varNo,2),1);
            else
                primalSol = x(varNo);
            end

            primalSol = G*primalSol + h;
            
            sol = spotsdpsol(pr,info,pr.variables,primalSol);
       end
       
       function sol = minimizeDualForm(solver,pr,objective)

           [A,b,c,K,G,h] = spotsdpsolver.sdpToDualSedumiFormat(pr,objective);
           
           [x,y,info] = sedumi(A,b,c,K,solver.options);


           if info.dinf || info.pinf,
               primalSol = NaN*ones(size(pr.variables,1),1);
           else
               primalSol = y;
           end
           
           sol = spotsdpsol(pr,info,pr.variables,G*primalSol+h);
       end
    end
end
