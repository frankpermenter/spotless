function [pr,objective,answer] = test_sdpProjection(pr,n)
    A = randn(n,n);
    A = A + A';
    [V,D] = eig(A);
    
    answer = norm(A-V*(D.*(D>=0))*inv(V),'fro');
    
    [pr,p] = pr.newFree(nchoosek(n+1,2));
    P = mss_v2s(p);
    pr = pr.withPSD(P);
    
    [pr,l0] = pr.newFree(1);
    
    pr = pr.withLor([l0 ; A(:)-P(:)]);

    objective = l0;
end