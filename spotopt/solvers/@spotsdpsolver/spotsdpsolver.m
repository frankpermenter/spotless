classdef  spotsdpsolver 
    methods (Abstract)
        %inst = isInstalled(solvr);
        can = canSolve(solvr,sdp);
        sol = minimize(solvr,sdp,objective);
    end
    
    methods (Static)
         function [A,b] = decompLinear(lin,vall)
            [veq,peq,Ceq] = decomp(lin);
            constant = ~any(peq~=0,2);
            cnsti = find(constant);

            b = -Ceq(:,cnsti);
            Aeq = Ceq(:,~constant)*peq(~constant,:);
            
            veqIndices = match(vall,veq);
            
            % T*vall = veq;
            T = sparse(1:length(veq),veqIndices,ones(length(veq),1),length(veq),length(vall));
            A = Aeq*T;
        end
        function [As,bs] = linearToSedumi(lin,vall,varNo,KvarCnt)
            [A,bs] = spotsdpsolver.decompLinear(lin,vall);
            [i,j,s] = find(A);
            As = sparse(i,varNo(j),s,size(A,1),KvarCnt);
        end
    end
end