classdef solver_tests < matlab.unittest.TestCase %#ok<*PROP>
    % SOLVER_TESTS  Tests for the osqp.Solver pure-MATLAB ADMM backend.

    properties
        P
        q
        A
        l
        u
        n
        m
        tol
    end

    methods (TestMethodSetup)
        function setup_problem(testCase)
            testCase.P = sparse([11 0; 0 0]);
            testCase.q = [3; 4];
            testCase.A = sparse([-1 0; 0 -1; -1 -3; 2 5; 3 4]);
            testCase.u = [0; 0; -15; 100; 80];
            testCase.l = -1e30 * ones(5, 1);
            testCase.n = 2;
            testCase.m = 5;
            testCase.tol = 1e-4;
        end
    end

    methods (Test)
        function test_basic_qp(testCase)
            % Solve a small QP and verify optimality conditions.
            solver = osqp.Solver();
            solver.setup(testCase.P, testCase.q, testCase.A, ...
                testCase.l, testCase.u, ...
                'verbose', false, 'linear_solver', 'matlab_ldl', ...
                'eps_abs', 1e-8, 'eps_rel', 1e-8);
            res = solver.solve();
            testCase.verifyEqual(res.info.status, 'solved');

            % Check KKT: Px + q + A'y = 0
            Pfull = testCase.P + testCase.P' - diag(diag(testCase.P));
            grad = Pfull * res.x + testCase.q + testCase.A' * res.y;
            testCase.verifyLessThan(norm(grad, inf), testCase.tol);

            % Feasibility
            Ax = testCase.A * res.x;
            testCase.verifyGreaterThanOrEqual(Ax + testCase.tol, testCase.l);
            testCase.verifyLessThanOrEqual(Ax - testCase.tol, testCase.u);
        end

        function test_factory_matlab_backend(testCase)
            % Verify factory function with 'matlab' backend.
            solver = osqp('backend', 'matlab');
            testCase.verifyTrue(isa(solver, 'osqp.Solver'));
            solver.setup(testCase.P, testCase.q, testCase.A, ...
                testCase.l, testCase.u, 'verbose', false, ...
                'linear_solver', 'matlab_ldl', ...
                'eps_abs', 1e-6, 'eps_rel', 1e-6);
            res = solver.solve();
            testCase.verifyEqual(res.info.status, 'solved');
        end

        function test_update_q(testCase)
            solver = osqp.Solver();
            solver.setup(testCase.P, testCase.q, testCase.A, ...
                testCase.l, testCase.u, 'verbose', false, ...
                'linear_solver', 'matlab_ldl', ...
                'eps_abs', 1e-8, 'eps_rel', 1e-8);
            solver.update('q', [10; 20]);
            res = solver.solve();
            testCase.verifyEqual(res.info.status, 'solved');
            testCase.verifyEqual(numel(res.x), testCase.n);
        end

        function test_update_bounds(testCase)
            solver = osqp.Solver();
            solver.setup(testCase.P, testCase.q, testCase.A, ...
                testCase.l, testCase.u, 'verbose', false, ...
                'linear_solver', 'matlab_ldl', ...
                'eps_abs', 1e-8, 'eps_rel', 1e-8);
            solver.update('l', -50 * ones(testCase.m, 1), ...
                'u', 1000 * ones(testCase.m, 1));
            res = solver.solve();
            testCase.verifyEqual(res.info.status, 'solved');
        end

        function test_warm_start(testCase)
            solver = osqp.Solver();
            solver.setup(testCase.P, testCase.q, testCase.A, ...
                testCase.l, testCase.u, 'verbose', false, ...
                'linear_solver', 'matlab_ldl', ...
                'eps_abs', 1e-8, 'eps_rel', 1e-8);
            res1 = solver.solve();

            solver.warm_start('x', res1.x, 'y', res1.y);
            res2 = solver.solve();
            testCase.verifyEqual(res2.info.status, 'solved');
            % Warm started should converge in fewer iterations
            testCase.verifyLessThanOrEqual(res2.info.iter, res1.info.iter);
        end

        function test_cold_start(testCase)
            solver = osqp.Solver();
            solver.setup(testCase.P, testCase.q, testCase.A, ...
                testCase.l, testCase.u, 'verbose', false, ...
                'linear_solver', 'matlab_ldl', ...
                'eps_abs', 1e-8, 'eps_rel', 1e-8);
            solver.solve();
            solver.cold_start();
            res = solver.solve();
            testCase.verifyEqual(res.info.status, 'solved');
        end

        function test_update_settings(testCase)
            solver = osqp.Solver();
            solver.setup(testCase.P, testCase.q, testCase.A, ...
                testCase.l, testCase.u, 'verbose', false, ...
                'linear_solver', 'matlab_ldl');
            solver.update_settings('max_iter', 10);
            res = solver.solve();
            testCase.verifyLessThanOrEqual(res.info.iter, 10);
        end

        function test_get_dimensions(testCase)
            solver = osqp.Solver();
            solver.setup(testCase.P, testCase.q, testCase.A, ...
                testCase.l, testCase.u, 'verbose', false, ...
                'linear_solver', 'matlab_ldl');
            [n, m] = solver.get_dimensions();
            testCase.verifyEqual(n, testCase.n);
            testCase.verifyEqual(m, testCase.m);
        end

        function test_current_settings(testCase)
            solver = osqp.Solver();
            solver.setup(testCase.P, testCase.q, testCase.A, ...
                testCase.l, testCase.u, 'verbose', false, ...
                'linear_solver', 'matlab_ldl', 'max_iter', 500);
            s = solver.current_settings();
            testCase.verifyEqual(s.max_iter, 500);
            testCase.verifyFalse(s.verbose);
        end

        function test_unconstrained(testCase)
            % min 0.5*x'Px + q'x  with no constraints
            P = sparse([4 1; 1 2]);
            q = [1; 1];
            solver = osqp.Solver();
            solver.setup(P, q, sparse(0, 2), [], [], ...
                'verbose', false, 'linear_solver', 'matlab_ldl', ...
                'eps_abs', 1e-10, 'eps_rel', 1e-10);
            res = solver.solve();
            testCase.verifyEqual(res.info.status, 'solved');
            x_exact = -(P \ q);
            testCase.verifyEqual(res.x, x_exact, 'AbsTol', 1e-6);
        end

        function test_polishing(testCase)
            solver = osqp.Solver();
            solver.setup(testCase.P, testCase.q, testCase.A, ...
                testCase.l, testCase.u, ...
                'verbose', false, 'linear_solver', 'matlab_ldl', ...
                'polishing', true, ...
                'eps_abs', 1e-6, 'eps_rel', 1e-6);
            res = solver.solve();
            testCase.verifyEqual(res.info.status, 'solved');
            testCase.verifyTrue(res.info.status_polish ~= 0);
        end

        function test_equality_constraints(testCase)
            % min 0.5*x'*x  s.t. x1 + x2 = 1
            P = speye(2);
            q = zeros(2, 1);
            A = sparse([1 1]);
            l = 1; u = 1;
            solver = osqp.Solver();
            solver.setup(P, q, A, l, u, 'verbose', false, ...
                'linear_solver', 'matlab_ldl', ...
                'eps_abs', 1e-10, 'eps_rel', 1e-10);
            res = solver.solve();
            testCase.verifyEqual(res.info.status, 'solved');
            testCase.verifyEqual(res.x, [0.5; 0.5], 'AbsTol', 1e-6);
        end

        function test_options_object(testCase)
            % Pass an osqp.Options object to setup.
            opts = osqp.Options();
            opts.verbose = false;
            opts.max_iter = 2000;
            opts.eps_abs = 1e-8;
            opts.eps_rel = 1e-8;
            opts.linear_solver = 'matlab_ldl';

            solver = osqp.Solver();
            solver.setup(testCase.P, testCase.q, testCase.A, ...
                testCase.l, testCase.u, opts);
            res = solver.solve();
            testCase.verifyEqual(res.info.status, 'solved');
        end

        function test_results_struct_fields(testCase)
            solver = osqp.Solver();
            solver.setup(testCase.P, testCase.q, testCase.A, ...
                testCase.l, testCase.u, 'verbose', false, ...
                'linear_solver', 'matlab_ldl');
            res = solver.solve();
            testCase.verifyTrue(isfield(res, 'x'));
            testCase.verifyTrue(isfield(res, 'y'));
            testCase.verifyTrue(isfield(res, 'prim_inf_cert'));
            testCase.verifyTrue(isfield(res, 'dual_inf_cert'));
            testCase.verifyTrue(isfield(res, 'info'));
            testCase.verifyTrue(isfield(res.info, 'iter'));
            testCase.verifyTrue(isfield(res.info, 'status'));
            testCase.verifyTrue(isfield(res.info, 'status_val'));
            testCase.verifyTrue(isfield(res.info, 'obj_val'));
            testCase.verifyTrue(isfield(res.info, 'pri_res'));
            testCase.verifyTrue(isfield(res.info, 'dua_res'));
            testCase.verifyTrue(isfield(res.info, 'run_time'));
        end
    end
end
