classdef spotprog
    properties 
        longname = '';
        sdp=[];         % Current SDP Problem Description
        
        numVar    = 0;  % Number of variables.
        
        % Mapping taking sdp variables into user facing variables.
        variableExpr = [];
        preProc = {};
        log = {};
    end
    
    methods ( Access = private )        
        function nm = varName(prog)
            nm = [prog.sdp.name 'var'];
        end
        
        function [prog,v] = newVariables(prog,var_sdp)
            prog.variableExpr = [ prog.variableExpr ; var_sdp ];
            v = msspoly(prog.varName,[length(var_sdp) prog.numVar]);
            prog.numVar = prog.numVar + length(var_sdp);
        end
                
        function sdpVar = toSDPVariables(pr,var)
            sdpVar = subs(var,pr.variables,pr.variableExpr);
        end
        
    end

    methods
        function prog = spotprog(shortname,longname)
            if nargin >= 1,
                prog.sdp = spotsdp(shortname);
            else
                prog.sdp = spotsdp;
            end
            
            if nargin >= 2,
                prog.longname = longname;
            else
                prog.longname = 'default';
            end
        end
        
        function v = variables(pr)
            v = msspoly(pr.varName,pr.numVar);
        end
        
        function pr = withEqs(pr,e)
            pr.sdp = pr.sdp.withEqs(pr.toSDPVariables(e));
        end
        function pr = withPos(pr,e)
            pr.sdp = pr.sdp.withPos(pr.toSDPVariables(e));
        end
        function pr = withLor(pr,e)
            pr.sdp = pr.sdp.withLor(pr.toSDPVariables(e));
        end
        function pr = withRLor(pr,e)
            pr.sdp = pr.sdp.withRLor(pr.toSDPVariables(e));
        end
        function pr = withPSD(pr,e)
            pr.sdp = pr.sdp.withPSD(pr.toSDPVariables(e));
        end
        function pr = withBlkPSD(pr,e)
            pr.sdp = pr.sdp.withBlkPSD(pr.toSDPVariables(e));
        end
        
        
        function [pr,v] = newFree(pr,dim)
            [pr.sdp,s] = pr.sdp.newFree(dim);
            [pr,v] = newVariables(pr,s);
        end
        
        function [pr,v] = newPos(pr,dim)
            [pr.sdp,s] = pr.sdp.newPos(dim);
            [pr,v] = newVariables(pr,s);
        end

        function [pr,v] = newLor(pr,dim)
            [pr.sdp,s] = pr.sdp.newLor(dim);
            [pr,v] = newVariables(pr,s(:));
            v = reshape(v,size(s));
        end
        
        function [pr,v] = newRLor(pr,dim)
            [pr.sdp,s] = pr.sdp.newRLor(dim);
            [pr,v] = newVariables(pr,s(:));
            v = reshape(v,size(s));
        end
        
        function [pr,v] = newPSD(pr,dim)
            [pr.sdp,s] = pr.sdp.newPSD(dim);
            [pr,v] = newVariables(pr,mss_s2v(s));
            v = mss_v2s(v);
        end
        
        
        function [pr,v] = newBlkPSD(pr,dim)
            [pr.sdp,s] = pr.sdp.newBlkPSD(dim);
            [pr,v] = newVariables(pr,s(:));
            v = reshape(v,size(s));
        end
        
        function pr = setPreProc(pr,preproc)
        %
        %   pr = addPreProc(pr,preproc)
        %
        %   pr      -- a spotprog
        %   preproc -- a SpotProgPreProc object.
        %
            pr.preProc{end+1} = preproc;
        end
        
        
        function [sdp,A,b,obj] = prepare(pr,objective)
            objective = subs(objective,pr.variables,pr.variableExpr);
            
            [A,b] = spot_decomp_linear(pr.variableExpr,pr.sdp.variables);
            b = -b;
            
            %  A*(Ai*z+bi)+b
            sdp = pr.sdp;
            for i = 1:length(pr.preProc)
                [sdp,Ai,bi,log{i}] = pr.preProc{i}.preProcess(sdp);
                A = A*Ai;
                b = b+A*bi;
            end
            
            obj = subs(objective,pr.variables,b+A*sdp.variables);
        end
        
        function sol = minimize(pr,solver,objective)
        %
        %  sol = pr.minimize(objective,solver)
        %  [sol,err] = pr.minimize(objective,solver)
        %
        %  objective -- 1-by-1 msspoly linear in pr.variables
        %  solver    -- spotsolver
        %
        %  sol -- If the given solver can handle the problem type, sol is
        %         a spotsol.  If the program cannot be handled 
        %         err = 1 for two output arguments and an error is
        %         thrown o.w.
            ;
            % First solve the SDP
            [sdp,A,b,obj] = prepare(pr,objective);
        
            int_sol = solver.minimize(sdp,obj);
            
            sol = spotsdpsol(pr,int_sol.info,pr.variables,...
                             int_sol.primalSolution,A,b);
        end
        
    end

end
