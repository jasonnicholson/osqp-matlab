classdef basic_tests < matlab.unittest.TestCase
    % BASIC_TESTS Basic QP tests for OSQP v1.0.0

    properties
        P
        q
        A
        u
        l
        m
        n
        solver
        tol
    end

    methods(TestMethodSetup)
        function setup_problem(testCase)
            testCase.P = sparse([11 0; 0 0]);
            testCase.q = [3; 4];
            testCase.A = sparse([-1 0; 0 -1; -1 -3; 2 5; 3 4]);
            testCase.u = [0; 0; -15; 100; 80];
            testCase.l = -1e30 * ones(5, 1);
            testCase.n = 2;
            testCase.m = 5;
            testCase.tol = 1e-4;

            testCase.solver = osqp;
            testCase.solver.setup(testCase.P, testCase.q, ...
                testCase.A, testCase.l, testCase.u, ...
                'verbose', false, 'polishing', true, ...
                'eps_abs', 1e-8, 'eps_rel', 1e-8);
        end
    end

    methods(TestMethodTeardown)
        function teardown(testCase)
            delete(testCase.solver);
        end
    end

    methods(Test)
        function test_basic_qp(testCase)
            [x_ref, y_ref, obj_ref] = load_high_accuracy('test_basic_QP');
            res = testCase.solver.solve();
            testCase.verifyEqual(res.x, x_ref, 'AbsTol', testCase.tol);
            testCase.verifyEqual(res.y, y_ref, 'AbsTol', testCase.tol);
            testCase.verifyEqual(res.info.obj_val, obj_ref, 'AbsTol', testCase.tol);
        end

        function test_update_q(testCase)
            [x_ref, y_ref, obj_ref] = load_high_accuracy('test_update_q');
            testCase.solver.update('q', [10; 20]);
            res = testCase.solver.solve();
            testCase.verifyEqual(res.x, x_ref, 'AbsTol', testCase.tol);
            testCase.verifyEqual(res.y, y_ref, 'AbsTol', testCase.tol);
            testCase.verifyEqual(res.info.obj_val, obj_ref, 'AbsTol', testCase.tol);
        end

        function test_update_l(testCase)
            [x_ref, y_ref, obj_ref] = load_high_accuracy('test_update_l');
            testCase.solver.update('l', -50 * ones(testCase.m, 1));
            res = testCase.solver.solve();
            testCase.verifyEqual(res.x, x_ref, 'AbsTol', testCase.tol);
            testCase.verifyEqual(res.y, y_ref, 'AbsTol', testCase.tol);
            testCase.verifyEqual(res.info.obj_val, obj_ref, 'AbsTol', testCase.tol);
        end

        function test_update_u(testCase)
            [x_ref, y_ref, obj_ref] = load_high_accuracy('test_update_u');
            testCase.solver.update('u', 1000 * ones(testCase.m, 1));
            res = testCase.solver.solve();
            testCase.verifyEqual(res.x, x_ref, 'AbsTol', testCase.tol);
            testCase.verifyEqual(res.y, y_ref, 'AbsTol', testCase.tol);
            testCase.verifyEqual(res.info.obj_val, obj_ref, 'AbsTol', testCase.tol);
        end

        function test_update_bounds(testCase)
            [x_ref, y_ref, obj_ref] = load_high_accuracy('test_update_bounds');
            testCase.solver.update('l', -50 * ones(testCase.m, 1), ...
                                   'u', 1000 * ones(testCase.m, 1));
            res = testCase.solver.solve();
            testCase.verifyEqual(res.x, x_ref, 'AbsTol', testCase.tol);
            testCase.verifyEqual(res.y, y_ref, 'AbsTol', testCase.tol);
            testCase.verifyEqual(res.info.obj_val, obj_ref, 'AbsTol', testCase.tol);
        end

        function test_update_max_iter(testCase)
            testCase.solver.update_settings('max_iter', 80);
            res = testCase.solver.solve();
            testCase.verifyEqual(res.info.status_val, ...
                osqp.constant.OSQP_MAX_ITER_REACHED);
        end

        function test_update_early_termination(testCase)
            testCase.solver.update_settings('check_termination', 0, ...
                                            'max_iter', 500);
            res = testCase.solver.solve();
            testCase.verifyEqual(res.info.iter, 500);
        end

        function test_update_rho(testCase)
            res0 = testCase.solver.solve();
            testCase.solver.update_settings('rho', 0.7);
            res1 = testCase.solver.solve();
            testCase.verifyEqual(res1.info.obj_val, res0.info.obj_val, ...
                'AbsTol', testCase.tol);
        end

        function test_update_time_limit(testCase)
            testCase.solver.update_settings('time_limit', 1e-6);
            res = testCase.solver.solve();
            testCase.verifyEqual(res.info.status_val, ...
                osqp.constant.OSQP_TIME_LIMIT_REACHED);
        end

        function test_upper_triangular_P(testCase)
            % Solve with full P, should get same result as triu(P)
            solver2 = osqp;
            P_full = testCase.P + tril(testCase.P, -1)';
            solver2.setup(P_full, testCase.q, ...
                testCase.A, testCase.l, testCase.u, ...
                'verbose', false, 'polishing', true, ...
                'eps_abs', 1e-8, 'eps_rel', 1e-8);
            res1 = testCase.solver.solve();
            res2 = solver2.solve();
            testCase.verifyEqual(res1.x, res2.x, 'AbsTol', testCase.tol);
            delete(solver2);
        end
    end
end
