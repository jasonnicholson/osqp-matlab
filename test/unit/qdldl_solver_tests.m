classdef qdldl_solver_tests < matlab.unittest.TestCase %#ok<*PROP>
    % QDLDL_SOLVER_TESTS  Tests for QDLDLSolver using a mock QDLDLFactorization.

    properties
        P
        q
        A
        l
        u
        n
        m
        tol
        mockPath
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

            % Add mock QDLDLFactorization to path
            testCase.mockPath = fullfile(fileparts(fileparts( ...
                mfilename('fullpath'))), 'mocks');
            addpath(testCase.mockPath);
        end
    end

    methods (TestMethodTeardown)
        function teardown(testCase)
            rmpath(testCase.mockPath);
        end
    end

    methods (Test)
        function test_qdldl_solver_basic(testCase)
            % QDLDLSolver can construct, solve, and produce correct results.
            K = osqp.LinearSolver.buildKKT( ...
                testCase.P, testCase.A, ones(testCase.m, 1), 1e-6, ...
                testCase.n, testCase.m);
            solver = osqp.linsys.QDLDLSolver(K);
            testCase.verifyTrue(isa(solver, 'osqp.linsys.QDLDLSolver'));

            b = ones(testCase.n + testCase.m, 1);
            x = solver \ b;
            testCase.verifyEqual(numel(x), testCase.n + testCase.m);

            % Verify correctness: K*x ≈ b
            Kfull = K + triu(K, 1)';
            testCase.verifyEqual(Kfull * x, b, 'AbsTol', 1e-8);
        end

        function test_qdldl_solver_refactorize(testCase)
            % refactorize rebuilds the factorization with new K.
            K = osqp.LinearSolver.buildKKT( ...
                testCase.P, testCase.A, ones(testCase.m, 1), 1e-6, ...
                testCase.n, testCase.m);
            solver = osqp.linsys.QDLDLSolver(K);

            % Modify rho and refactorize
            K2 = osqp.LinearSolver.buildKKT( ...
                testCase.P, testCase.A, 2*ones(testCase.m, 1), 1e-6, ...
                testCase.n, testCase.m);
            solver.refactorize(K2);

            b = ones(testCase.n + testCase.m, 1);
            x = solver \ b;
            K2full = K2 + triu(K2, 1)';
            testCase.verifyEqual(K2full * x, b, 'AbsTol', 1e-8);
        end

        function test_qdldl_solver_via_solver(testCase)
            % End-to-end: osqp.Solver with linear_solver='qdldl'.
            solver = osqp.Solver();
            solver.setup(testCase.P, testCase.q, testCase.A, ...
                testCase.l, testCase.u, 'verbose', false, ...
                'linear_solver', 'qdldl', ...
                'eps_abs', 1e-8, 'eps_rel', 1e-8);
            res = solver.solve();
            testCase.verifyEqual(res.info.status, 'solved');
            testCase.verifyEqual(numel(res.x), testCase.n);
        end

        function test_qdldl_solver_equality_constraint(testCase)
            % min 0.5*x'*x  s.t. x1 + x2 = 1
            P = speye(2);
            q = zeros(2, 1);
            A = sparse([1 1]);
            l = 1; u = 1;
            solver = osqp.Solver();
            solver.setup(P, q, A, l, u, 'verbose', false, ...
                'linear_solver', 'qdldl', ...
                'eps_abs', 1e-10, 'eps_rel', 1e-10);
            res = solver.solve();
            testCase.verifyEqual(res.info.status, 'solved');
            testCase.verifyEqual(res.x, [0.5; 0.5], 'AbsTol', 1e-6);
        end

        function test_qdldl_solver_update_settings_rho(testCase)
            % Updating rho triggers refactorize with qdldl backend.
            solver = osqp.Solver();
            solver.setup(testCase.P, testCase.q, testCase.A, ...
                testCase.l, testCase.u, 'verbose', false, ...
                'linear_solver', 'qdldl');
            solver.update_settings('rho', 0.5);
            res = solver.solve();
            testCase.verifyEqual(res.info.status, 'solved');
        end
    end
end
