function [pr,objective,answer] = example_sdpProjection(n)
    A = randn(n,n);
    A = A + A';
    [V,D] = eig(A);
    
    answer = norm(A-V*(D.*(D>=0))*inv(V),'fro');
    
    pr = spotprog;
    [pr,p] = pr.new('free',nchoosek(n+1,2));
    P = mss_v2s(p);
    pr = pr.with('psd',P);
    
    [pr,l0] = pr.new('free',1);
    
    pr = pr.with('lor',[l0 ; A(:)-P(:)]);

    objective = l0;
end