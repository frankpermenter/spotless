function [prog,psdApprox] = newDD(prog,n)
%
%   [prout,psd] = pr.newSDD(n)
%
%   Adds a "scaled diagonally dominant matrix variable.
%   O.w. known as a factor-width 2 PSD matrix.
%
    M = nchoosek(n,2);
    [prog,d] = prog.newPos(n);
    [prog,p] = prog.newPos(M);
    [prog,m] = prog.newPos(M);
    
    psdP = [ 1  1 1 ]'*p';
    psdM = [ 1 -1 1 ]'*m';
    

    N2 = mss_s2v(reshape(1:n^2,n,n)); 
    
    % Enumerate unique pairs (i,j), i < j.
    [i,j]=ind2sub([n n],N2);
    ondiag = i == j;
    id = i(ondiag);
    jd = i(ondiag);

    offdiag = i ~= j;
    i = i(offdiag);
    j = j(offdiag);
    
    % Produce row/column that each 2-by-2 should be mapped to.
    I = [ i' ; i' ; j' ];
    J = [ i' ; j' ; j' ];

    % Produce the matrix index that corresponds to these.
    Ind = sub2ind([n n],I(:),J(:));
    VInd = mss_match(N2,Ind);
    S = sparse(VInd,(1:length(VInd))',ones(length(VInd),1),length(N2),length(VInd));
    
    Indd = sub2ind([n n],id,jd);
    VIndd = mss_match(N2,Indd);
    Sd = sparse(VIndd,(1:length(VIndd))',ones(length(VIndd),1),length(N2),length(VIndd));
    
    psdApprox = S*psdP(:) + S*psdM(:) + Sd*d;
end