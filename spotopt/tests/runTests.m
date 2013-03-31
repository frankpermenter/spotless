function [] = runTests(solver)
    if nargin < 1,
        solver = spotsolversedumi();
    end
    
    [err,answer,pinf,dinf] = runFeasibleTest(solver,...
                                             @test_maxEig,20);
    
    if abs(err)/abs(answer) > 1e-6, error('Test failed.'); end
    
    
    [err,answer,pinf,dinf] = runFeasibleTest(solver,...
                                             @test_sdpProjection,20);
    if (abs(err)/abs(answer)) > 1e-7, error('Test failed.'); end
end