classdef polishing_tests < matlab.unittest.TestCase
    % POLISHING_TESTS Test polishing for OSQP v1.0.0

    properties
        tol
    end

    methods(TestMethodSetup)
        function setup_problem(testCase)
            testCase.tol = 1e-4;
        end
    end

    methods(Test)
        function test_polish_simple(testCase)
            P = sparse([11 0; 0 0]);
            q = [3; 4];
            A = sparse([-1 0; 0 -1; -1 -3; 2 5; 3 4]);
            u = [0; 0; -15; 100; 80];
            l = -1e30 * ones(5, 1);

            [x_ref, y_ref, obj_ref] = load_high_accuracy('test_polish_simple');

            solver = osqp;
            solver.setup(P, q, A, l, u, ...
                'verbose', false, 'polishing', true, ...
                'polish_refine_iter', 4, ...
                'eps_abs', 1e-8, 'eps_rel', 1e-8);
            res = solver.solve();

            testCase.verifyEqual(res.info.status_polish, ...
                osqp.constant.OSQP_POLISH_SUCCESS);
            testCase.verifyEqual(res.x, x_ref, 'AbsTol', testCase.tol);
            testCase.verifyEqual(res.y, y_ref, 'AbsTol', testCase.tol);
            testCase.verifyEqual(res.info.obj_val, obj_ref, 'AbsTol', testCase.tol);

            delete(solver);
        end

        function test_polish_unconstrained(testCase)
            rng(4);
            n = 30;
            Pt = sprandn(n, n, 0.6);
            P = Pt * Pt' + speye(n);
            q = randn(n, 1);
            A = sparse(0, n);
            l = zeros(0, 1);
            u = zeros(0, 1);

            [x_ref, ~, obj_ref] = load_high_accuracy('test_polish_unconstrained');

            solver = osqp;
            solver.setup(P, q, A, l, u, ...
                'verbose', false, 'polishing', true, ...
                'polish_refine_iter', 4, ...
                'eps_abs', 1e-8, 'eps_rel', 1e-8);
            res = solver.solve();

            testCase.verifyEqual(res.x, x_ref, 'AbsTol', testCase.tol);
            testCase.verifyEqual(res.info.obj_val, obj_ref, 'AbsTol', testCase.tol);

            delete(solver);
        end

        function test_polish_random(testCase)
            rng(6);
            n = 30;
            m = 50;
            Pt = sprandn(n, n, 0.6);
            P = Pt * Pt' + speye(n);
            q = randn(n, 1);
            A = sprandn(m, n, 0.8);
            l = -2 * ones(m, 1);
            u =  2 * ones(m, 1);

            [x_ref, y_ref, obj_ref] = load_high_accuracy('test_polish_random');

            solver = osqp;
            solver.setup(P, q, A, l, u, ...
                'verbose', false, 'polishing', true, ...
                'polish_refine_iter', 4, ...
                'eps_abs', 1e-8, 'eps_rel', 1e-8);
            res = solver.solve();

            testCase.verifyEqual(res.info.status_polish, ...
                osqp.constant.OSQP_POLISH_SUCCESS);
            testCase.verifyEqual(res.x, x_ref, 'AbsTol', testCase.tol);
            testCase.verifyEqual(res.y, y_ref, 'AbsTol', testCase.tol);
            testCase.verifyEqual(res.info.obj_val, obj_ref, 'AbsTol', testCase.tol);

            delete(solver);
        end
    end
end
