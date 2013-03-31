function [pr,objective,answer] = example_maxEig(n)
    A = randn(n,n);
    P = A'*A;
    answer = max(eig(P));
    
    pr = spotprog;
    [pr,tau] = pr.new('free',1);
    pr = pr.with('psd',tau*eye(n)-P);
    objective = tau;
end