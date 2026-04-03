classdef warm_start_tests < matlab.unittest.TestCase
    % WARM_START_TESTS Test warm starting for OSQP v1.0.0

    properties
        P
        q
        A
        u
        l
        m
        n
        opts
        tol
    end

    methods(TestMethodSetup)
        function setup_problem(testCase)
            rng(4);
            testCase.n = 100;
            testCase.m = 200;
            Pt = sprandn(testCase.n, testCase.n, 0.6);
            testCase.P = Pt * Pt' + speye(testCase.n);
            testCase.q = randn(testCase.n, 1);
            testCase.A = sprandn(testCase.m, testCase.n, 0.8);
            testCase.l = -2 * ones(testCase.m, 1);
            testCase.u =  2 * ones(testCase.m, 1);
            testCase.tol = 1e-4;
            testCase.opts = struct('verbose', false, ...
                'eps_abs', 1e-8, 'eps_rel', 1e-8, ...
                'warm_starting', true);
        end
    end

    methods(Test)
        function test_warm_start(testCase)
            % Cold start (warm start from zero)
            solver = osqp;
            solver.setup(testCase.P, testCase.q, ...
                testCase.A, testCase.l, testCase.u, testCase.opts);
            solver.warm_start('x', zeros(testCase.n, 1), ...
                              'y', zeros(testCase.m, 1));
            res_cold = solver.solve();

            % Solve from optimal warm start
            solver.warm_start('x', res_cold.x, 'y', res_cold.y);
            res_warm = solver.solve();

            % Warm start should converge faster
            testCase.verifyLessThan(res_warm.info.iter, res_cold.info.iter);
            testCase.verifyLessThan(res_warm.info.iter, 10);

            delete(solver);
        end
    end
end
