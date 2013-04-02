function [dl,dobj] = toDual(pr,pobj)
%
%  [dl,dobj] = toDual(pr,pobj)
%
%  pr   -- spotsdp
%  pobj -- 1-by-1 linear msspoly, function of pr.variables.
%  dual -- spotsdp 
%  dobj -- 1-by-1 msspoly
%
%
%  pr is given by the mixed primal dual description:
%
%  f free,  x in K1,  A1*x + A2*f = b,   F1*x + F2*f - g in K2  
%
%  pobj = <c1,x> + <c2,f>  (to be minimized)
%    
%  dl is the dual conic program derived from the Lagrangian:
%  
%  <c1,x> + <c2,f> - <lambda,x> + <y,A1*x+A2*f-b> - <s,F1*x+F2*f-g>
%
%  i.e.: 
%
%  y free,  c1'+A1'*y-F1'*s in K1*,  s in  K2*, A2'*y-F2'*s + c2' =0
%
%  with dobj = <g,s> - <y,b>.
    ;
    f = pr.freeVariables;
    x = [ pr.posVariables ; pr.lorVariables ; pr.rlorVariables ; pr.psdVariables];
    
    % x in K1, determine inner product for this cone.
    Q1 = inner_product(length(pr.posVariables)+...
                       length(pr.lorVariables)+...
                       length(pr.rlorVariables),pr.psdDim);
    

    dl = spotsdp;
    
    if isempty(pr.equations)
        A1Ty = zeros(length(x),1);
        A2Ty = zeros(length(f),1);
        bTy = 0;
    else
        [A,b] = spot_decomp_linear(pr.equations,[x;f]);
        A1 = A(:,1:length(x)); A2 = A(:,length(x)+1:end);

        % Construct dual with variables.
        [dl,y] = dl.newFree(length(b));
        
        A2Ty = A2'*y;
        A1Ty = Q1\A1'*y;
        bTy  = b'*y;
    end    
    
    [c,~] = spot_decomp_linear(pobj,[x;f]);
    c1 = c(:,1:length(x));  c2 = c(:,length(x)+1:end);
    c1T = c1';
    
    % Next, determine the cone K2*
    [pcstr,pdim] = pr.posConstraints();
    [lcstr,ldim] = pr.lorConstraints();
    [rcstr,rdim] = pr.rlorConstraints();
    [scstr,sdim] = pr.psdConstraints();
    
    [dl,pdl] = dl.newPos(pdim);
    [dl,ldl] = dl.newLors(ldim);
    [dl,rdl] = dl.newRLors(rdim);
    [dl,sdl] = dl.newPSDs(sdim);
    
    s = [ pdl ; ldl ; rdl ; sdl ];
    
    % Determine Equations
    K2 = [ pcstr ; lcstr ; rcstr ; scstr];
    
    % Need to determine the inner-product for this cone.
    % Particularly the PSD cone.
    Q2 = inner_product(length(pcstr)+length(lcstr)+length(rcstr),sdim);
    
    if ~isempty(K2)
        [F,g] = spot_decomp_linear(K2,[x;f]);
        F1 = F(:,1:length(x)); F2 = F(:,length(x)+1:end);
    
        F1Ts = Q1\F1'*Q2*s;
        F2Ts = F2'*Q2*s;
        gTs = g'*Q2*s;
    else
        F1Ts = zeros(size(c1T));
        F2Ts = zeros(size(c2'));
        gTs = 0;
    end
    
    dl = dl.withEqs(A2Ty-F2Ts+c2');
    
    % Finally determine cone K1*
    K1star = c1T+A1Ty-F1Ts;
    
    ipos  = (1:pr.numPos)';
    ilor  = pr.numPos+(1:pr.numLor)';
    irlor = pr.numPos+pr.numLor+(1:pr.numRLor)';
    ipsd  = pr.numPos+pr.numLor+pr.numRLor+(1:pr.numPSD)';
    
    dl = dl.withPos(K1star(ipos));
    dl = dl.withLor(K1star(ilor),pr.lorDim);
    dl = dl.withRLor(K1star(irlor),pr.rlorDim);
    dl = dl.withPSD(K1star(ipsd),pr.psdDim);
    
    dobj = gTs-bTy;
    
    function Q = inner_product(n,dim)
        I = 2*ones(sum(spotsdp.psdDimToNo(dim)),1);
        
        i0 = 0;
        for i = 1:length(dim)
            ii = cumsum(1:dim(i));
            I(i0+ii) = I(i0+ii) - 1;
            i0 = i0 + spotsdp.psdDimToNo(dim(i));
        end
        
        Q = blkdiag(speye(n),spdiags(I,0,length(I),length(I)));
    end
end