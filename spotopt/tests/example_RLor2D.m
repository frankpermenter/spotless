function [pr,objective,answer] = example_RLor2D()
    N = 10;
    answer = 0;
    pr = spotprog;
    [pr,f] = pr.newFree(N);
    for i = 1:N
        A = randn(2,2);
        A = A + A';
        [V,D] = eig(A);
        answer = answer + norm(A-V*(D.*(D>=0))*inv(V),'fro');

        [pr,r] = pr.newRLor(3);
        P = [ sqrt(2)*r(1) r(3) ; r(3) sqrt(2)*r(2) ];
        pr = pr.withLor([f(i);A(:)-P(:)]);
    end
    objective = 1+sum(f);
    answer = answer + 1;
end