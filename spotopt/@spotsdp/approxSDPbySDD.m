function [approx,oldVariables,newVariables] = approxSDPbySDD(pr)
%
% [approx,oldVar,newExpr] = approxSDPbySDD(pr)
%
% pr      -- spotsdp.
% oldVar  -- Variables replaced.
% newExpr -- Expression of old variables in new variables.
% approx  -- New spotsdp.
% 
% Approx is a new SDP with each sdp constraint and variable
% replaced by a scaled diagonally dominant constraint.
    
    % Step one:  Construct new program with new cone variables
    approx = spotsdp;
    
    fold = pr.freeVariables;
    [approx,fnew] = approx.newFree(pr.numFree);
    
    pold = pr.posVariables;
    [approx,pnew] = approx.newPos(pr.numPos);
    
    lold = pr.lorVariables;
    [approx,lnew] = approx.newLors(pr.lorDim);
    
    rold = pr.rlorVariables;
    [approx,rnew] = approx.newRLors(pr.rlorDim);
    
    psdold = pr.psdVariables;
    psddim = pr.psdDim;
    psdnew = [];
    for i = 1:length(pr.psdDim)
        [approx,psdApprox] = newSDD(approx,pr.psdDim(i));
        psdnew = [ psdnew ; psdApprox ];
    end

    oldVariables = [ fold ; pold ; lold ; rold ; psdold ];
    newVariables = [ fnew ; pnew ; lnew ; rnew ; psdnew ];

    % Step Two: Walk through constraints, replacing new variables.

    [pcstr,dim] = pr.posConstraints();
    approx = approx.withPos(subs(pcstr,oldVariables,newVariables));

    [lcstr,dim] = pr.lorConstraints();
    approx = approx.withLor(subs(lcstr,oldVariables,newVariables),dim);

    [rcstr,dim] = pr.rlorConstraints();
    approx = approx.withRLor(subs(rcstr,oldVariables,newVariables),dim);

    [~,dim] = pr.psdConstraints();
    for i = 1:length(dim)
        scstr = subs(pr.psdConstraint(i),oldVariables,newVariables);
        [approx,slack] = newSDD(approx,pr.psdCstrDim(i));
        approx = approx.withEqs(scstr - slack);
    end
end