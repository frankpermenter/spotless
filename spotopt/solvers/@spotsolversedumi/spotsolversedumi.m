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
        
        function sol = optimize(solvr,pr,objective)
            if nargin < 2, objective = 0; end
            if nargin < 3, options = struct('fid',0); end

            objective = msspoly(objective);
            if ~pr.realLinearInDec(objective)
                error('Objective must be real and linear in dec. variables.');
            end
            
            %  First, construct structure with counts of SeDuMi
            %  variables.
            K = struct();
            K.f = pr.freeNum;
            K.l = pr.posNum;
            K.q = pr.lorDim;
            K.r = pr.rlorDim;
            K.s = pr.psdDim;
            
            KvarCnt = K.f+K.l+sum(K.q)+sum(K.r)+sum(K.s.^2);
            
            
            v = [ pr.freeVariables
                  pr.posVariables
                  pr.lorVariables
                  pr.rlorVariables ];
            
            vpsd = pr.psdVariables;
            
            vall = [v;vpsd];
            
            % Assign column numbers to v.
            psdVarNo = zeros(1,length(vpsd));
            psdVarNoSymm = zeros(1,length(vpsd));
            psdVarOff = 0;    % Progress in variables, storing
                              % upper triangle.
            psdRedVarOff = 0; % Progress in variables, storing
                              % entire matrix.
            for i = 1:length(pr.psdDim)
                n = pr.psdDim(i);
                m = n*(n+1)/2;
                psdVarNo(psdVarOff + (1:m)) = psdRedVarOff+mss_s2v(reshape(1:n^2,n,n));
                psdVarNoSymm(psdVarOff + (1:m)) = psdRedVarOff+mss_s2v(reshape(1:n^2,n,n)');
                psdVarOff = psdVarOff + m;
                psdRedVarOff = psdRedVarOff + n^2;
            end
            
            varNo = [ 1:length(v) length(v)+psdVarNo];
            varNoSymm = [ 1:length(v) length(v)+psdVarNoSymm];

            [b,A] = spotsolversedumi.linearToSedumi(pr.equations);
            [~,c] = spotsolversedumi.linearToSedumi(objective);
            
            dualObjective = b'*pr.dualVariables;

            [x,y,info] = sedumi(A,b,c,K,options);
            
            if info.pinf, 
                primalSol = NaN*ones(size(length(varNo),1));
            else
                primalSol = x(varNo);
            end
            
            if info.dinf
                dualSol = NaN*ones(size(y));
                dualSlack = NaN*ones(size(primalSol));
            else
                dualSol = y;
                z = c'-A'*y;
                dualSlack = (z(varNo)+z(varNoSymm))/2;
            end

            sol = spotsdpsol(pr,info,...
                             objective,dualObjective,...
                             primalSol,dualSol,dualSlack);
        end
    end
    
    methods (Static, private)
        function [As,bs] = linearToSedumi(lin,vall,varNo,KvarCnt)
            [A,bs] = spotsdpsolver.decompLinear(lin,vall);
            [i,j,s] = find(A);
            As = sparse(i,varNo(j),s,size(A,1),KvarCnt);
        end
    end
end
