function [prog,psdApprox] = newSDD(prog,n)
%
%   [prout,psd] = pr.newSDD(n)
%
%   Adds a "scaled diagonally dominant matrix variable.
%   O.w. known as a factor-width 2 PSD matrix.
%
    M = nchoosek(n,2);
    [prog,r] = prog.newRLor([3 M]);
    
    % Map to 2-by-2 SDP matrices in upper triangular format.
    psd2 = [ sqrt(2) 0 0 ; 0 0 1 ; 0 sqrt(2) 0]*r;
    
    N2 = mss_s2v(reshape(1:n^2,n,n)); 
    
    % Enumerate unique pairs (i,j), i < j.
    [i,j]=ind2sub([n n],N2);
    offdiag = i ~= j;
    i = i(offdiag);
    j = j(offdiag);
    
    % Produce row/column that each 2-by-2 should be mapped to.
    I = [ i' ; i' ; j' ];
    J = [ i' ; j' ; j' ];
    
    % Produce the matrix index that corresponds to these.
    Ind = sub2ind([n n],I(:),J(:));
    % Now map matrix indices into the upper triangular vector format.
    VInd = mss_match(N2,Ind);
    S = sparse(VInd,(1:length(VInd))',ones(length(VInd),1),length(N2),length(VInd));
    
    psdApprox = S*psd2(:);
end