function [pr,objective,answer] = example_Chebyshev()
    N=101;
    d = 20;
    xx = linspace(-1,1,N).';
    Phi = repmat(xx,1,d+1).^repmat((0:d),N,1);
    
    pr = spotprog;
    [pr,c] = pr.newFree(d+1);
    
    err = abs(xx) - Phi*c;
    
    [pr,t] = pr.newFree(1);
    
    pr = pr.withPos(t - err);
    pr = pr.withPos(t + err);
    answer = 0.0141;
    objective = t;
end