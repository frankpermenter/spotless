classdef spotprog
    properties 
        longname = '';
        sdp=[];         % Current SDP Problem Description
        
        numVar    = 0;  % Number of variables.
        
        % Mapping taking sdp variables into user facing variables.
        var_expr = [];
        
    end
    
    methods (Static)

    end
    
    methods ( Access = private )        
        function nm = varName(prob)
            nm = [prob.sdp.name 'var'];
        end
        
        function [prob,v] = newVariables(prob,var_sdp)
            prob.var_expr = [ prob.var_expr ; var_sdp ];
            v = msspoly(prob.varName,[length(var_sdp) prob.numVar]);
            prob.numVar = prob.numVar + length(var_sdp);
        end
                
        function sdpVar = toSDPVariables(prob,var)
            sdpVar = subs(var,pr.variables,pr.variableExpr);
        end
        
    end

    methods
        function prob = spotprog(shortname,longname)
            if nargin >= 1,
                prob.sdp = spotsdp(shortname);
            else
                prob.sdp = spotsdp;
            end
            
            if nargin >= 2,
                prob.longname = longname;
            else
                prob.longname = 'default';
            end
        end
        
        function v = variables(pr)
            v = msspoly(prob.varName,prob.numVar);
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

        function cstr = newFree(dim)
            if ~spot_hasSize(dim,[1 1]``) || ~spot_isIntGE(dim,1)
                error('Dimension must be scalar positive integer.');
            end
            [prob.sdp,l] = newLor(prob.sdp,dim);
            [prob,v] = newVariables(prob,l(:));
            
            v = reshape(v,size(l));
        end

        function cstr = newPos(dim)
            if ~spot_hasSize(dim,[1 1]``) || ~spot_isIntGE(dim,1)
                error('Dimension must be scalar positive integer.');
            end
            [prob.sdp,l] = newLor(prob.sdp,dim);
            [prob,v] = newVariables(prob,l(:));
            
            v = reshape(v,size(l));
        end
        
        function [prob,l] = newLor(prob,dim)
            if spot_hasSize(dim,[1 1]), dim = [dim 1]; end
            if ~spot_hasSize(dim,[1 2]) || ~spot_isIntGE(dim,1)
                error('Dimension must be 1x2 positive integer.');
            end
            
            [prob.sdp,l] = newLor(prob.sdp,dim);
            [prob,v] = newVariables(prob,l(:));
            
            v = reshape(v,size(l));
        end
        
        function [prob,l] = newLor(prob,dim)
            if spot_hasSize(dim,[1 1]), dim = [dim 1]; end
            if ~spot_hasSize(dim,[1 2]) || ~spot_isIntGE(dim,1)
                error('Dimension must be 1x2 positive integer.');
            end
            
            [prob.sdp,l] = newLor(prob.sdp,dim);
            [prob,v] = newVariables(prob,l(:));
            
            v = reshape(v,size(l));
        end
        
        function [prob,V] = newPSD(prob,dim)
            if ~spot_hasSize(dim,[1 1]) || ~spot_isIntGE(dim,1)
                error('Dimension must be scalar positive integer.');
            end
            
            n = spotsdp.psdDimToNo(dim);
            
            [prob.sdp,Q] = newPSD(prob.sdp,dim);
            
            q = mss_s2v(Q);
            
            [prob,v] = newVariables(prob,q);
            
            V = mss_v2s(v);
        end
        
         function [prob,Vs] = newBlkPSD(prob,dim)
            if ~spot_hasSize(dim,[1 2]) || ~spot_isIntGE(dim,1)
                error('Dimension must be 1x2 positive integer.');
            end
            
            n = spotsdp.psdDimToNo(dim(1));
            
            [prob.sdp,Qs] = newBlkPSD(prob.sdp,dim);
            
            [prob,Vs] = newVariables(prob,Qs(:));
            
            Vs = reshape(Vs,size(Qs));
        end
        
        
        function sol = minimize(objective,solver)
        %
        %  sol = pr.minimize(objective,solver)
        %  [sol,err] = pr.minimize(objective,solver)
        %
        %  objective -- 1-by-1 msspoly linear in pr.variables
        %  solver    -- spotsolver
        %
        %  sol -- If the given solver can handle the problem type, sol is
        %         a spotsol.  If the problem cannot be handled 
        %         err = 1 for two output arguments and an error is
        %         thrown o.w.
        end
        
    end

end
