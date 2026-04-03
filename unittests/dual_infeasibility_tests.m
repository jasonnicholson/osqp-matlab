classdef dual_infeasibility_tests < matlab.unittest.TestCase
    % DUAL_INFEASIBILITY_TESTS Test dual infeasibility detection for OSQP v1.0.0

    properties
        opts
        tol
    end

    methods(TestMethodSetup)
        function setup_problem(testCase)
            testCase.opts = struct('verbose', false, ...
                'eps_prim_inf', 1e-15, ...
                'polishing', true);
            testCase.tol = 1e-4;
        end
    end

    methods(Test)
        function test_dual_infeasible_lp(testCase)
            P = sparse(2, 2);
            q = [2; -1];
            A = speye(2);
            l = [0; 0];
            u = Inf(2, 1);

            solver = osqp;
            solver.setup(P, q, A, l, u, testCase.opts);
            res = solver.solve();

            testCase.verifyEqual(res.info.status_val, ...
                osqp.constant('OSQP_DUAL_INFEASIBLE'));

            % Verify dual infeasibility certificate
            data = load(fullfile(fileparts(mfilename('fullpath')), ...
                'solutions', 'test_dual_infeasibility.mat'));
            cert = res.dual_inf_cert / norm(res.dual_inf_cert);
            testCase.verifyEqual(cert, ...
                data.lp_normalized_dual_inf_cert_correct(:), ...
                'AbsTol', testCase.tol);

            delete(solver);
        end

        function test_dual_infeasible_qp(testCase)
            P = sparse(diag([4; 0]));
            q = [0; 2];
            A = sparse([1 1; -1 1]);
            l = [-Inf; -Inf];
            u = [2; 3];

            solver = osqp;
            solver.setup(P, q, A, l, u, testCase.opts);
            res = solver.solve();

            testCase.verifyEqual(res.info.status_val, ...
                osqp.constant('OSQP_DUAL_INFEASIBLE'));

            % Verify dual infeasibility certificate
            data = load(fullfile(fileparts(mfilename('fullpath')), ...
                'solutions', 'test_dual_infeasibility.mat'));
            cert = res.dual_inf_cert / norm(res.dual_inf_cert);
            testCase.verifyEqual(cert, ...
                data.qp_normalized_dual_inf_cert_correct(:), ...
                'AbsTol', testCase.tol);

            delete(solver);
        end

        function test_primal_dual_infeasible(testCase)
            P = sparse(2, 2);
            q = [-1; -1];
            A = sparse([1 -1; -1 1; 1 0; 0 1]);
            l = [1; 1; 0; 0];
            u = Inf(4, 1);

            solver = osqp;
            solver.setup(P, q, A, l, u, testCase.opts);
            res = solver.solve();

            testCase.verifyTrue( ...
                res.info.status_val == osqp.constant('OSQP_PRIMAL_INFEASIBLE') || ...
                res.info.status_val == osqp.constant('OSQP_DUAL_INFEASIBLE'));

            delete(solver);
        end
    end
end
