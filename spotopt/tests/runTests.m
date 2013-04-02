function [] = runTests(solver)
    if nargin < 1,
        solver = spotsolversedumi();
    end
    
    [err,answer,pinf,dinf] = runFeasibleTest(solver,...
                                             @test_maxEig,spotsdp,20);
    
    if any(abs(err)/abs(answer) > 1e-6), error('Test failed.'); end
    
    
    [err,answer,pinf,dinf] = runFeasibleTest(solver,...
                                             @test_sdpProjection,spotsdp,20);
    if any((abs(err)/abs(answer)) > 1e-7), error('Test failed.'); end
end