function [pr,objective,answer] = example_maxEig(n)
    A = randn(n,n);
    P = A'*A;
    answer = max(eig(P));

    pr = spotsdp;
    [pr,tau] = pr.newFree(1);
    pr = pr.withPSD(tau*eye(n)-P);
    objective = tau;
end