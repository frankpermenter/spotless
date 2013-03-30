classdef  spotsdpsolver 
    methods (Abstract)
        %inst = isInstalled(solvr);
        can = canSolve(solvr,sdp);
        sol = minimize(solvr,sdp,objective);
    end
    
    methods (Static)       
        function [As,bs] = linearToSedumi(lin,vall,varNo,KvarCnt)
            [A,bs] = spot_decomp_linear(lin,vall);
            [i,j,s] = find(A);
            As = sparse(i,varNo(j),s,size(A,1),KvarCnt);
        end
        
        function [psdVarNo,psdVarNoSymm] = upperTriToFullVarNo(npsd,psdDim)
        % Assign column numbers to v.
            psdVarNo = zeros(1,npsd);
            psdVarNoSymm = zeros(1,npsd);
            psdVarOff = 0;    % Progress in variables, storing
                              % upper triangle.
            psdRedVarOff = 0; % Progress in variables, storing
                              % entire matrix.
            for i = 1:length(psdDim)
                n = psdDim(i);
                m = n*(n+1)/2;
                psdVarNo(psdVarOff + (1:m)) = psdRedVarOff+mss_s2v(reshape(1:n^2,n,n));
                psdVarNoSymm(psdVarOff + (1:m)) = psdRedVarOff+mss_s2v(reshape(1:n^2,n,n)');
                psdVarOff = psdVarOff + m;
                psdRedVarOff = psdRedVarOff + n^2;
            end
        end
       
        function [A,b,c,K,G,h,varNo] = sdpToPrimalSedumiFormat(pr,objective)
            objective = msspoly(objective);
            if ~realLinearInDec(pr,objective)
                error('Objective must be real and linear in dec. variables.');
            end
            
            user_variables = pr.variables;
            
            [pr,G,h] = pr.standardPrimalWithFree();

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
 
            [psdVarNo] = spotsdpsolver.upperTriToFullVarNo(length(vpsd),pr.psdDim);
            
            varNo = [ 1:length(v) length(v)+psdVarNo];
            
           [A,b] = spotsdpsolver.linearToSedumi(pr.equations,vall,varNo,KvarCnt);
           [c,~] = spotsdpsolver.linearToSedumi(objective,vall,varNo,KvarCnt);
        end
        
        function [A,b,c,K,G,h] = sdpToDualSedumiFormat(pr,objective)
            user_variables = pr.variables;
           
            objective = msspoly(objective);
            if ~realLinearInDec(pr,objective)
                error('Objective must be real and linear in dec. variables.');
            end
            
            if ~spotsdp.isStandardDualForm(pr)
                error(['Right now the program must be in standard dual ' ...
                       'format.']);
                
                [pr,G,h] = pr.standardDual();
            else
                G = speye(size(user_variables,1));
                h = sparse([],[],[],size(user_variables,1),size(user_variables,2));
            end
           
            objective = subs(objective,user_variables,G*pr.variables+h);
           
            % minimize b'y
            %          c-A'y in K.
            %
            % (1) Decide on the appropriate cone sizes.
            % (2) Construct A, c.
            
            %  First, construct structure with counts of SeDuMi
            %  variables.
            K = struct();
            K.f = 0;
            K.r = 0;
            
            pos = pr.posCnst;
            K.l = length(pr.posCnst);
            
            K.q = [];
            lor = [];
            for i = 1:length(pr.lorCnst)
                lor = [ lor ; pr.lorCnst{i}(:) ];
                K.q = [ K.q size(pr.lorCnst{i},1)*ones(1,size(pr.lorCnst{i},2))];
            end
            
            v = [ pos ; lor ];
            
            K.s = [];
            vpsd = [];
            for i = 1:length(pr.psdCnst)
                vpsd = [ vpsd ; pr.psdCnst{i}(:) ];
                K.s = [ K.s spotsdp.psdNoToDim(size(pr.psdCnst{i},1))*ones(1,size(pr.psdCnst{i},2))];
            end
            
            KvarCnt = K.f+K.l+sum(K.q)+sum(K.r)+sum(K.s.^2);
            
            [psdVarNo,psdVarNoSymm] = spotsdpsolver.upperTriToFullVarNo(length(vpsd),K.s);
            
            varNo = [ 1:length(v) length(v)+psdVarNo];
            varNoSymm = [ 1:length(v) length(v)+psdVarNoSymm];
            varNoDiag = varNo(find(varNo == varNoSymm));
            
            
            %  I need to isolate the term:  c - A'y in K where K is
            %  the cone of the above dimensions.
            
            vall = [ v ; vpsd ];
            
            [mAT,cm] = spot_decomp_linear(vall,pr.variables);
            [b,~]    = spot_decomp_linear(objective,pr.variables);
            b = -b.';
            
            [i,j,s] = find(mAT);
            
            mAT = sparse(varNo(i),j,s,KvarCnt,size(mAT,2));
            keep = varNo(i) ~= varNoSymm(i);
            mAT = mAT+sparse(varNoSymm(i(keep)),...
                             j(keep),s(keep),KvarCnt,size(mAT,2));
            
            [i,j,s] = find(-cm);
            c = sparse(varNo(i),j,s,KvarCnt,1);
            keep = varNo(i) ~= varNoSymm(i);
            c = c + sparse(varNoSymm(i(keep)),j(keep),s(keep),KvarCnt,1);
            
            A = -mAT.';
        end
        
        
    end
end