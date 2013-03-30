classdef spotsolversdpnal < spotsdpsolver
    properties
        options = struct();
    end
    
    methods 
        
        function solvr = spotsolversdpnal(options)
            ;
            % TODO: check that all options are legal.
            if nargin > 0
                solvr.options = options;
            else
                solvr.options = struct();
            end
        end
        
        function can = canSolve(solvr,sdp)
            can = 1;
        end
        
       function sol = minimize(solver,pr,objective)
            error('Optimization in primal form not supported.');
       end
       
       function sol = minimizeDualForm(solver,pr,objective)
           [A,b,c,K,G,h] = spotsdpsolver.sdpToDualSedumiFormat(pr,objective);
           
           [blk,Asdpt,Csdpt,bsdpt] = read_sedumi_sdpnal_0(A,full(b),c,K);

           [obj,X,y,Z,infosdpt] = sdpnal(blk,Asdpt,Csdpt,bsdpt);
           info = struct('pinf',infosdpt.termcode == 1,...
                         'dinf',infosdpt.termcode == 2);
               
           if info.dinf || info.pinf,
               primalSol = NaN*ones(size(pr.variables,1),1);
           else
               primalSol = y;
           end
           
           sol = spotsdpsol(pr,info,pr.variables,G*primalSol+h);
       end
    end
end
