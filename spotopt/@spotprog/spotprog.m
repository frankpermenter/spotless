classdef spotprog
    properties 
        longname = '';
        sdp=[];         % Current SDP Problem Description
        
        numVar    = 0;  % Number of variables.
        
        % Mapping taking sdp variables into user facing variables.
        variableExpr = [];
        
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
        
        
        function pr = with(pr,cstr,varargin)
        %
        %   [pr,x] = pr.new(cstr,expr1,expr2,...)
        %
        %   cstr  -- spotcstr object.
        %   exprI -- msspoly depending on pr.variables.
        %   
        %
            if nargin ~= 3,
                error('Incorrect number of arugments.');
            end 
            
            for i = 1 :length(varargin)
                varargin{i} = pr.toSDPVariables(varargin{i});
            end
            
            if isa(cstr,'char')
                switch cstr
                  case 'pos',
                    pr.sdp = pr.sdp.withPos(varargin{1});
                  case 'lor',
                    pr.sdp = pr.sdp.withLor(varargin{1});
                  case 'rlor',
                    pr.sdp = pr.sdp.withRLor(varargin{1});
                  case 'psd',
                    pr.sdp = pr.sdp.withPSD(varargin{1});
                  case 'blkpsd',
                    pr.sdp = pr.sdp.withBlkPSD(varargin{1});
                  otherwise,
                    error('Constraint must be string: pos,lor,rlor,psd,blkpsd.');
                end
            else
                error('Constraint must be string.');
            end
        end
        
        
        function [prog,v] = new(prog,type,dim)
            if nargin ~= 3,
                error('Incorrect number of arugments.');
            end 
            
            if ~isa(type,'char') 
                error('Type must be string.');
            elseif ~ismember(type,{'free','pos','lor','rlor','psd','blkpsd'})
                error('Type must be string: free,pos,lor,rlor,psd,blkpsd.');
            end
            
            
            % Now standardize dimensions
            if ismember(type,{'free','pos'})
                if spot_hasSize(dim,[1 1])
                    dim = [ dim 1];
                end
            end
            
            if strcmp(type,'psd') && ~spot_hasSize(dim,[1 1])
                error(['Dimension for this type must be 1-by-1.']);
            elseif ~spot_hasSize(dim,[1 2])
                error(['Dimension for this type must be 1-by-2.']);
            end
            
            if ~spot_isIntGE(dim,1), 
                error(['Dimensions must be positive integers.']);
            elseif strcmp(type,'rlor') & dim(1) < 2,
                error(['First dimension for type rlor must be >= 2.']);
            end
                        
            switch type
              case 'free',
                [prog.sdp,l] = prog.sdp.newFree(dim);
              case 'pos',
                [prog.sdp,l] = prog.sdp.newPos(dim);
              case 'lor',
                [prog.sdp,l] = prog.sdp.newLor(dim);
              case 'rlor',
                [prog.sdp,l] = prog.sdp.newRLor(dim);
              case 'psd',
                [prog.sdp,Q] = newPSD(spotsdp.psdDimToNo(dim));
                l = mss_s2v(Q);
              case 'blkpsd',
                [prog.sdp,l] = prog.sdp.newBlkPSD(dim);
            end
            
            [prog,v] = newVariables(prog,l);
            
            if strcmp(v,'psd')
                v = mss_v2s(v);
            end
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
            objective = subs(objective,pr.variables,pr.variableExpr);
        
            int_sol = solver.minimize(pr.sdp,objective);
            
            [A,b] = spot_decomp_linear(pr.variableExpr,pr.sdp.variables);

            sol = spotsdpsol(pr,int_sol.info,pr.variables,...
                             int_sol.primalSolution,A,b);
        end
        
    end

end
