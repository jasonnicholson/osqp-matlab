classdef primal_infeasibility_tests < matlab.unittest.TestCase
    % PRIMAL_INFEASIBILITY_TESTS Test primal infeasibility detection for OSQP v1.0.0

    properties
        opts
        tol
    end

    methods(TestMethodSetup)
        function setup_problem(testCase)
            testCase.opts = struct('verbose', false, ...
                'eps_prim_inf', 1e-5, 'max_iter', 2500);
            testCase.tol = 1e-4;
        end
    end

    methods(Test)
        function test_primal_infeasible_problem(testCase)
            rng(4);
            n = 50;
            m = 500;
            Pt = sprandn(n, n, 0.6);
            P = Pt' * Pt;
            q = randn(n, 1);
            A = sprandn(m, n, 0.8);
            u = 3 + randn(m, 1);
            l = -3 + randn(m, 1);

            % Make infeasible: duplicate row with non-overlapping bounds
            nhalf = floor(n/2);
            A(nhalf, :) = A(nhalf + 1, :);
            l(nhalf) = u(nhalf + 1) + 10 * rand();
            u(nhalf) = l(nhalf) + 0.5;

            solver = osqp;
            solver.setup(P, q, A, l, u, testCase.opts);
            res = solver.solve();

            testCase.verifyEqual(res.info.status_val, ...
                osqp.constant('OSQP_PRIMAL_INFEASIBLE'));

            delete(solver);
        end
    end
end
