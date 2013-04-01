function [pr,objective,answer] = example_Chebyshev()
    N=101;
    d = 20;
    xx = linspace(-1,1,N).';
    Phi = repmat(xx,1,d+1).^repmat((0:d),N,1);
    
    pr = spotprog;
    [pr,c] = pr.new('free',d+1);
    
    err = abs(xx) - Phi*c;
    
    [pr,t] = pr.new('free',1);
    
    pr = pr.with('pos',t - err);
    pr = pr.with('pos',t + err);
    answer = 0.0141;
    objective = t;
end