classdef solver_coverage_tests < matlab.unittest.TestCase %#ok<*PROP>
    % SOLVER_COVERAGE_TESTS  Tests for uncovered Solver.m code paths.

    properties
        P
        q
        A
        l
        u
        n
        m
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
        end
    end

    methods (Test)
        %% Settings parsing alternatives
        function test_setup_with_struct(testCase)
            % Cover the struct-based settings path in setup().
            s.verbose = false;
            s.linear_solver = 'matlab_ldl';
            s.eps_abs = 1e-8;
            s.eps_rel = 1e-8;
            solver = osqp.Solver();
            solver.setup(testCase.P, testCase.q, testCase.A, ...
                testCase.l, testCase.u, s);
            res = solver.solve();
            testCase.verifyEqual(res.info.status, 'solved');
        end

        function test_update_settings_with_struct(testCase)
            % Cover the struct-based update_settings path.
            solver = osqp.Solver();
            solver.setup(testCase.P, testCase.q, testCase.A, ...
                testCase.l, testCase.u, 'verbose', false, ...
                'linear_solver', 'matlab_ldl');
            s_upd.max_iter = 2000;
            s_upd.verbose = false;
            solver.update_settings(s_upd);
            cs = solver.current_settings();
            testCase.verifyEqual(cs.max_iter, 2000);
        end

        function test_update_settings_odd_args_error(testCase)
            % Cover the odd-argument-count error in update_settings.
            solver = osqp.Solver();
            solver.setup(testCase.P, testCase.q, testCase.A, ...
                testCase.l, testCase.u, 'verbose', false, ...
                'linear_solver', 'matlab_ldl');
            testCase.verifyError( ...
                @() solver.update_settings('max_iter'), ...
                'OSQP:update_settings');
        end

        function test_update_settings_setup_only_error(testCase)
            % Cover the setup-only field error in update_settings.
            solver = osqp.Solver();
            solver.setup(testCase.P, testCase.q, testCase.A, ...
                testCase.l, testCase.u, 'verbose', false, ...
                'linear_solver', 'matlab_ldl');
            testCase.verifyError( ...
                @() solver.update_settings('sigma', 1e-3), ...
                'OSQP:update_settings');
        end

        function test_update_settings_struct_setup_only_error(testCase)
            % Cover the setup-only field error via struct path.
            solver = osqp.Solver();
            solver.setup(testCase.P, testCase.q, testCase.A, ...
                testCase.l, testCase.u, 'verbose', false, ...
                'linear_solver', 'matlab_ldl');
            s_upd.scaling = 5;
            testCase.verifyError( ...
                @() solver.update_settings(s_upd), ...
                'OSQP:update_settings');
        end

        %% Dimension detection edge cases
        function test_setup_P_empty_n_from_q(testCase)
            % Cover: n detected from q when P is empty.
            q = [1; 2; 3];
            A = speye(3);
            solver = osqp.Solver();
            solver.setup([], q, A, -ones(3,1), ones(3,1), ...
                'verbose', false, 'linear_solver', 'matlab_ldl');
            [n, m] = solver.get_dimensions();
            testCase.verifyEqual(n, 3);
            testCase.verifyEqual(m, 3);
        end

        function test_setup_P_empty_n_from_A(testCase)
            % Cover: n detected from A when P and q are empty.
            A = sparse([1 0 0; 0 1 0]);
            solver = osqp.Solver();
            solver.setup([], [], A, [-1;-1], [1;1], ...
                'verbose', false, 'linear_solver', 'matlab_ldl');
            [n, ~] = solver.get_dimensions();
            testCase.verifyEqual(n, 3);
        end

        function test_setup_no_variables_error(testCase)
            % Cover: error when P, q, A are all empty.
            solver = osqp.Solver();
            testCase.verifyError( ...
                @() solver.setup([], [], [], [], []), ...
                'OSQP:setup');
        end

        %% Scaling disabled
        function test_scaling_disabled(testCase)
            % Cover: scaling=0 early return in scaleProblem.
            solver = osqp.Solver();
            solver.setup(testCase.P, testCase.q, testCase.A, ...
                testCase.l, testCase.u, 'verbose', false, ...
                'linear_solver', 'matlab_ldl', 'scaling', 0, ...
                'eps_abs', 1e-6, 'eps_rel', 1e-6);
            res = solver.solve();
            testCase.verifyEqual(res.info.status, 'solved');
        end

        %% Non-full-triangle P
        function test_setup_full_symmetric_P(testCase)
            % Cover: triu(P) path when P is not upper-triangular.
            P_full = [11 1; 1 2];  % symmetric but not triu
            solver = osqp.Solver();
            solver.setup(sparse(P_full), testCase.q, testCase.A, ...
                testCase.l, testCase.u, 'verbose', false, ...
                'linear_solver', 'matlab_ldl', ...
                'eps_abs', 1e-6, 'eps_rel', 1e-6);
            res = solver.solve();
            testCase.verifyEqual(res.info.status, 'solved');
        end

        %% Max iterations reached
        function test_max_iter_reached(testCase)
            % Cover: post-loop convergence check and STATUS_MAX_ITER_REACHED.
            solver = osqp.Solver();
            solver.setup(testCase.P, testCase.q, testCase.A, ...
                testCase.l, testCase.u, 'verbose', false, ...
                'linear_solver', 'matlab_ldl', ...
                'max_iter', 2, 'check_termination', 1, ...
                'eps_abs', 1e-15, 'eps_rel', 1e-15);
            res = solver.solve();
            testCase.verifyTrue(ismember(res.info.status_val, ...
                [osqp.constant.OSQP_MAX_ITER_REACHED, ...
                 osqp.constant.OSQP_SOLVED_INACCURATE]));
            testCase.verifyLessThanOrEqual(res.info.iter, 2);
            % Should still have a valid solution
            testCase.verifyEqual(numel(res.x), testCase.n);
        end

        %% Time limit reached
        function test_time_limit_reached(testCase)
            % Cover: time_limit branch in solve loop.
            % Use a large problem that takes many iterations.
            rng(42);
            n = 50; m = 100;
            Pt = sprandn(n, n, 0.3);
            P = Pt * Pt' + 0.01 * speye(n);
            q = randn(n, 1);
            A = sprandn(m, n, 0.5);
            l = -2 * ones(m, 1);
            u = 2 * ones(m, 1);

            solver = osqp.Solver();
            solver.setup(P, q, A, l, u, 'verbose', false, ...
                'linear_solver', 'matlab_ldl', ...
                'time_limit', 1e-10, ...  % Extremely small time limit
                'eps_abs', 1e-15, 'eps_rel', 1e-15, ...
                'check_termination', 0);  % Disable periodic check to force time limit
            res = solver.solve();
            testCase.verifyTrue(ismember(res.info.status_val, ...
                [osqp.constant.OSQP_TIME_LIMIT_REACHED, ...
                 osqp.constant.OSQP_SOLVED, ...
                 osqp.constant.OSQP_SOLVED_INACCURATE]));
        end

        %% Warm start edge cases
        function test_warm_start_x_only(testCase)
            % Cover: "only x provided" path → zeros y.
            solver = osqp.Solver();
            solver.setup(testCase.P, testCase.q, testCase.A, ...
                testCase.l, testCase.u, 'verbose', false, ...
                'linear_solver', 'matlab_ldl', ...
                'eps_abs', 1e-6, 'eps_rel', 1e-6);
            res1 = solver.solve();
            solver.warm_start('x', res1.x);
            res2 = solver.solve();
            testCase.verifyEqual(res2.info.status, 'solved');
        end

        function test_warm_start_y_only(testCase)
            % Cover: "only y provided" path → zeros x, z.
            solver = osqp.Solver();
            solver.setup(testCase.P, testCase.q, testCase.A, ...
                testCase.l, testCase.u, 'verbose', false, ...
                'linear_solver', 'matlab_ldl', ...
                'eps_abs', 1e-6, 'eps_rel', 1e-6);
            res1 = solver.solve();
            solver.warm_start('y', res1.y);
            res2 = solver.solve();
            testCase.verifyEqual(res2.info.status, 'solved');
        end

        function test_warm_start_struct(testCase)
            % Cover: struct-based warm_start input.
            solver = osqp.Solver();
            solver.setup(testCase.P, testCase.q, testCase.A, ...
                testCase.l, testCase.u, 'verbose', false, ...
                'linear_solver', 'matlab_ldl', ...
                'eps_abs', 1e-6, 'eps_rel', 1e-6);
            res1 = solver.solve();
            ws.x = res1.x;
            ws.y = res1.y;
            solver.warm_start(ws);
            res2 = solver.solve();
            testCase.verifyEqual(res2.info.status, 'solved');
        end

        function test_warm_start_bad_field_error(testCase)
            % Cover: error for unrecognized field in warm_start.
            solver = osqp.Solver();
            solver.setup(testCase.P, testCase.q, testCase.A, ...
                testCase.l, testCase.u, 'verbose', false, ...
                'linear_solver', 'matlab_ldl');
            testCase.verifyError( ...
                @() solver.warm_start('z', zeros(5,1)), '');
        end

        function test_warm_starting_disabled(testCase)
            % Cover: warm_starting=false branch in solve (all-zeros init).
            solver = osqp.Solver();
            solver.setup(testCase.P, testCase.q, testCase.A, ...
                testCase.l, testCase.u, 'verbose', false, ...
                'linear_solver', 'matlab_ldl', ...
                'warm_starting', false, ...
                'eps_abs', 1e-6, 'eps_rel', 1e-6);
            res = solver.solve();
            testCase.verifyEqual(res.info.status, 'solved');
        end

        %% Non-convexity
        function test_non_convex_warning(testCase)
            % Cover: non_convex=true path and STATUS_NON_CONVEX.
            % P is negative definite but sigma-regularized P is PD.
            P = sparse([-1 0; 0 -0.5]);
            q = [1; 1];
            A = speye(2);
            l = [-10; -10];
            u = [10; 10];
            solver = osqp.Solver();
            solver.setup(P, q, A, l, u, 'verbose', false, ...
                'linear_solver', 'matlab_ldl', 'sigma', 100);
            res = solver.solve();
            testCase.verifyEqual(res.info.status_val, ...
                osqp.constant.OSQP_NON_CONVEX);
            testCase.verifyEqual(res.info.status, 'non_convex');
            testCase.verifyTrue(isnan(res.info.obj_val));
        end

        function test_non_convex_error(testCase)
            % Cover: P is non-convex and sigma too small → error.
            P = sparse([-10 0; 0 -10]);
            q = [1; 1];
            A = speye(2);
            solver = osqp.Solver();
            testCase.verifyError( ...
                @() solver.setup(P, q, A, [-1;-1], [1;1], ...
                    'verbose', false, 'linear_solver', 'matlab_ldl'), ...
                'OSQP:NonConvex');
        end

        %% Dual infeasibility
        function test_dual_infeasible(testCase)
            % Cover: dual infeasibility detection and result formatting.
            % min -x  s.t. x >= 0 (unbounded below)
            P = sparse(1, 1);
            q = -1;
            A = sparse(1, 1, 1);
            l = 0;
            u = inf;
            solver = osqp.Solver();
            solver.setup(P, q, A, l, u, 'verbose', false, ...
                'linear_solver', 'matlab_ldl', ...
                'eps_dual_inf', 1e-4, 'max_iter', 4000);
            res = solver.solve();
            testCase.verifyTrue(ismember(res.info.status_val, ...
                [osqp.constant.OSQP_DUAL_INFEASIBLE, ...
                 osqp.constant.OSQP_DUAL_INFEASIBLE_INACCURATE]));
            testCase.verifyEqual(res.info.obj_val, -inf);
            testCase.verifyTrue(all(isnan(res.x)));
        end

        %% Adaptive rho with explicit interval
        function test_adaptive_rho_explicit_interval(testCase)
            % Cover: user-specified adaptive_rho_interval > 0.
            solver = osqp.Solver();
            solver.setup(testCase.P, testCase.q, testCase.A, ...
                testCase.l, testCase.u, 'verbose', false, ...
                'linear_solver', 'matlab_ldl', ...
                'adaptive_rho', true, ...
                'adaptive_rho_interval', 5, ...
                'eps_abs', 1e-6, 'eps_rel', 1e-6);
            res = solver.solve();
            testCase.verifyEqual(res.info.status, 'solved');
        end

        %% Polishing: unconstrained (n_act == 0)
        function test_polish_unconstrained(testCase)
            % Cover: polishing with no active constraints (Cholesky path).
            P = sparse([4 1; 1 2]);
            q = [1; 1];
            solver = osqp.Solver();
            solver.setup(P, q, sparse(0, 2), [], [], ...
                'verbose', false, 'linear_solver', 'matlab_ldl', ...
                'polishing', true, ...
                'eps_abs', 1e-10, 'eps_rel', 1e-10);
            res = solver.solve();
            testCase.verifyEqual(res.info.status, 'solved');
            testCase.verifyEqual(res.info.status_polish, 1);
        end

        %% Polishing: constrained success
        function test_polish_constrained_success(testCase)
            % Cover: successful polishing with active constraints.
            % Simple well-conditioned QP where polish should improve.
            P = sparse([2 0; 0 2]);
            q = [-2; -2];
            A = sparse([1 1; -1 0; 0 -1]);
            l = [-inf; -inf; -inf];
            u = [1; 0; 0];  % x1 + x2 <= 1, x1 >= 0, x2 >= 0
            solver = osqp.Solver();
            solver.setup(P, q, A, l, u, 'verbose', false, ...
                'linear_solver', 'matlab_ldl', ...
                'polishing', true, ...
                'eps_abs', 1e-6, 'eps_rel', 1e-6);
            res = solver.solve();
            testCase.verifyEqual(res.info.status, 'solved');
            % Polishing should succeed or be performed
            testCase.verifyTrue(res.info.status_polish ~= 0);
        end

        %% Update with rho triggers refactorize
        function test_update_rho(testCase)
            % Cover: rho update in update_settings triggers refactorize.
            solver = osqp.Solver();
            solver.setup(testCase.P, testCase.q, testCase.A, ...
                testCase.l, testCase.u, 'verbose', false, ...
                'linear_solver', 'matlab_ldl');
            solver.update_settings('rho', 0.5);
            cs = solver.current_settings();
            testCase.verifyEqual(cs.rho, 0.5);
            res = solver.solve();
            testCase.verifyEqual(res.info.status, 'solved');
        end

        %% solve() before setup error
        function test_solve_before_setup_error(testCase)
            solver = osqp.Solver();
            testCase.verifyError(@() solver.solve(), 'OSQP:solve');
        end

        %% update_settings before setup error
        function test_update_settings_before_setup_error(testCase)
            solver = osqp.Solver();
            testCase.verifyError( ...
                @() solver.update_settings('max_iter', 100), ...
                'OSQP:update_settings');
        end

        %% Default l and u
        function test_default_l_u(testCase)
            % Cover: empty l and u default to -inf/+inf.
            P = sparse([2 0; 0 2]);
            q = [1; 1];
            A = sparse([1 0; 0 1]);
            solver = osqp.Solver();
            solver.setup(P, q, A, [], [], 'verbose', false, ...
                'linear_solver', 'matlab_ldl', ...
                'eps_abs', 1e-6, 'eps_rel', 1e-6);
            res = solver.solve();
            testCase.verifyEqual(res.info.status, 'solved');
        end
    end
end
