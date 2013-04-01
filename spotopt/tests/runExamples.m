function [] = runExamples(solver)
    if nargin < 1,
        solver = spotsolversedumi();
    end
    
    [err,answer,pinf,dinf] = runFeasibleExample(solver,...
                                             @example_maxEig,20);
    
    if abs(err)/abs(answer) > 1e-6, error('Test failed.'); end
    
    
    [err,answer,pinf,dinf] = runFeasibleExample(solver,...
                                             @example_sdpProjection,20);
    if (abs(err)/abs(answer)) > 1e-7, error('Test failed.'); end

    [err,answer,pinf,dinf] = runFeasibleExample(solver,...
                                             @example_Chebyshev);
    if (abs(err)/abs(answer)) > 2e-2, error('Test failed.'); end

    [err,answer,pinf,dinf] = runFeasibleExample(solver,...
                                             @example_RLor2D);
    if (abs(err)/abs(answer)) > 1e-8, error('Test failed.'); end

end