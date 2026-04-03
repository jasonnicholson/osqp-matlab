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
            P = Pt * Pt' + speye(n);
            q = randn(n, 1);
            A = sprandn(m, n, 0.8);
            u = randn(m, 1);
            l = u + 10;  % infeasible: l > u

            solver = osqp;
            solver.setup(P, q, A, l, u, testCase.opts);
            res = solver.solve();

            testCase.verifyEqual(res.info.status_val, ...
                osqp.constant('OSQP_PRIMAL_INFEASIBLE'));

            % Verify primal infeasibility certificate
            data = load(fullfile(fileparts(mfilename('fullpath')), ...
                'solutions', 'test_primal_infeasibility.mat'));
            cert = res.prim_inf_cert / norm(res.prim_inf_cert);
            testCase.verifyEqual(abs(cert), ...
                abs(data.normalized_prim_inf_cert_correct(:)), ...
                'AbsTol', testCase.tol);

            delete(solver);
        end
    end
end
