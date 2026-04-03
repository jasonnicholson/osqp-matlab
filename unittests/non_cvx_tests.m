classdef non_cvx_tests < matlab.unittest.TestCase
    % NON_CVX_TESTS Test non-convex problem detection for OSQP v1.0.0

    methods(Test)
        function test_non_cvx_small_sigma(testCase)
            P = sparse([2 5; 5 1]);
            q = [3; 4];
            A = sparse([-1 0; 0 -1; -1 -3; 2 5; 3 4]);
            u = [0; 0; -15; 100; 80];
            l = -1e30 * ones(5, 1);

            solver = osqp;
            testCase.verifyError(@() solver.setup(P, q, A, l, u, ...
                'verbose', false, 'sigma', 1e-6), '');

            delete(solver);
        end
    end
end
