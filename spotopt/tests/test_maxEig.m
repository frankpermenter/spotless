function [pr,objective,answer] = test_maxEig(pr,n)
    A = randn(n,n);
    P = A'*A;
    answer = max(eig(P));

    [pr,tau] = pr.newFree(1);
    pr = pr.withPSD(tau*eye(n)-P);
    objective = tau;
end