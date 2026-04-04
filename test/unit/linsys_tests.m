classdef linsys_tests < matlab.unittest.TestCase
    % LINSYS_TESTS  Unit tests for osqp.LinearSolver and backends.

    properties (TestParameter)
        % Test with both small and medium problems
        problemSize = struct('small', 2, 'medium', 10)
    end

    methods (Test)
        function testBuildKKT(testCase)
            % Verify KKT assembly: upper triangular, correct size
            [P, A, rho_vec, sigma, n, m] = linsys_tests.makeQP(4);
            K = osqp.LinearSolver.buildKKT(triu(P), A, rho_vec, sigma, n, m);
            testCase.verifyTrue(istriu(K), ...
                'KKT matrix must be upper triangular');
            testCase.verifySize(K, [n+m, n+m]);
            % Verify the full symmetric KKT is correct
            Kfull = K + K' - diag(diag(K));
            testCase.verifyEqual(Kfull, Kfull', 'AbsTol', eps, ...
                'Full KKT must be symmetric');
        end

        function testMatlabLDLSolverSolve(testCase, problemSize)
            n = problemSize;
            [P, A, rho_vec, sigma, ~, m] = linsys_tests.makeQP(n);
            K = osqp.LinearSolver.buildKKT(triu(P), A, rho_vec, sigma, n, m);
            Kfull = K + K' - diag(diag(K));
            b = randn(n + m, 1);

            solver = osqp.linsys.MatlabLDLSolver(K);
            x = solver.solve(b);
            x_ref = Kfull \ b;
            testCase.verifyEqual(x, x_ref, 'AbsTol', 1e-10, ...
                'MatlabLDLSolver solve must match backslash');
        end

        function testMatlabLDLSolverBackslash(testCase, problemSize)
            n = problemSize;
            [P, A, rho_vec, sigma, ~, m] = linsys_tests.makeQP(n);
            K = osqp.LinearSolver.buildKKT(triu(P), A, rho_vec, sigma, n, m);
            Kfull = K + K' - diag(diag(K));
            b = randn(n + m, 1);

            solver = osqp.linsys.MatlabLDLSolver(K);
            x = solver \ b;
            x_ref = Kfull \ b;
            testCase.verifyEqual(x, x_ref, 'AbsTol', 1e-10, ...
                'Backslash (mldivide) must match direct solve');
        end

        function testMatlabLDLSolverRefactorize(testCase)
            n = 4;
            [P, A, rho_vec, sigma, ~, m] = linsys_tests.makeQP(n);
            K1 = osqp.LinearSolver.buildKKT(triu(P), A, rho_vec, sigma, n, m);
            solver = osqp.linsys.MatlabLDLSolver(K1);

            % Change rho and refactorize
            rho_vec2 = rho_vec * 2;
            K2 = osqp.LinearSolver.buildKKT(triu(P), A, rho_vec2, sigma, n, m);
            K2full = K2 + K2' - diag(diag(K2));
            solver.refactorize(K2);

            b = randn(n + m, 1);
            x = solver.solve(b);
            x_ref = K2full \ b;
            testCase.verifyEqual(x, x_ref, 'AbsTol', 1e-10, ...
                'Refactorized solver must match new K');
        end

        function testMatlabLDLSolverIndefiniteKKT(testCase)
            % The KKT matrix is indefinite (negative diagonal block).
            % Verify decomposition handles it correctly.
            n = 5;
            [P, A, rho_vec, sigma, ~, m] = linsys_tests.makeQP(n);
            K = osqp.LinearSolver.buildKKT(triu(P), A, rho_vec, sigma, n, m);
            Kfull = K + K' - diag(diag(K));

            % Verify Kfull is indeed indefinite
            eigvals = eig(full(Kfull));
            testCase.verifyTrue(any(eigvals < 0) && any(eigvals > 0), ...
                'KKT matrix should be indefinite');

            solver = osqp.linsys.MatlabLDLSolver(K);
            b = randn(n + m, 1);
            x = solver.solve(b);
            testCase.verifyEqual(Kfull * x, b, 'AbsTol', 1e-9, ...
                'Solution must satisfy K*x = b for indefinite K');
        end

        function testMatlabLDLNoConstraints(testCase)
            % Test with m=0 (unconstrained QP)
            n = 3;
            P = gallery('lehmer', n);
            A = sparse(0, n);
            rho_vec = zeros(0, 1);
            sigma = 1e-6;
            K = osqp.LinearSolver.buildKKT(triu(P), A, rho_vec, sigma, n, 0);
            testCase.verifySize(K, [n, n]);

            solver = osqp.linsys.MatlabLDLSolver(K);
            b = randn(n, 1);
            x = solver.solve(b);
            Kfull = K + K' - diag(diag(K));
            x_ref = Kfull \ b;
            testCase.verifyEqual(x, x_ref, 'AbsTol', 1e-12);
        end
    end

    methods (Static, Access = private)
        function [P, A, rho_vec, sigma, n, m] = makeQP(n)
            % Generate a random QP with n variables and n constraints
            rng(42 + n);  % deterministic
            M = randn(n);
            P = M' * M + 0.1 * eye(n);  % PSD Hessian
            A = randn(n, n);
            A = sparse(A);
            m = n;
            sigma = 1e-6;
            rho = 0.1;
            rho_vec = rho * ones(m, 1);
            % Make some equality constraints
            rho_vec(1) = rho * 1000;
        end
    end
end
