classdef update_matrices_tests < matlab.unittest.TestCase
    % UPDATE_MATRICES_TESTS Test P and A matrix updates for OSQP v1.0.0

    properties
        P
        P_new
        q
        A
        A_new
        u
        l
        m
        n
        solver
        tol
    end

    methods(TestMethodSetup)
        function setup_problem(testCase)
            rng(1);
            testCase.n = 5;
            testCase.m = 8;
            p = 0.7;

            Pt = sprandn(testCase.n, testCase.n, p);
            testCase.P = Pt * Pt' + speye(testCase.n);
            testCase.q = randn(testCase.n, 1);
            testCase.A = sprandn(testCase.m, testCase.n, p);
            testCase.l = -2 * ones(testCase.m, 1);
            testCase.u =  2 * ones(testCase.m, 1);
            testCase.tol = 1e-4;

            % New matrices with SAME sparsity pattern as originals
            % P_new: replace values in triu(P) pattern, keep PSD
            P_triu = triu(testCase.P);
            [pi, pj, ~] = find(P_triu);
            P_new_vals = randn(nnz(P_triu), 1);
            P_triu_new = sparse(pi, pj, P_new_vals, testCase.n, testCase.n);
            testCase.P_new = P_triu_new + triu(P_triu_new, 1)' + 5*speye(testCase.n);

            % A_new: replace values in A pattern
            [ai, aj, ~] = find(testCase.A);
            A_new_vals = randn(nnz(testCase.A), 1);
            testCase.A_new = sparse(ai, aj, A_new_vals, testCase.m, testCase.n);

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
        function test_solve(testCase)
            [x_ref, y_ref, obj_ref] = load_high_accuracy('test_solve');
            res = testCase.solver.solve();
            testCase.verifyEqual(res.x, x_ref, 'AbsTol', testCase.tol);
            testCase.verifyEqual(res.y, y_ref, 'AbsTol', testCase.tol);
            testCase.verifyEqual(res.info.obj_val, obj_ref, 'AbsTol', testCase.tol);
        end

        function test_update_P(testCase)
            [x_ref, y_ref, obj_ref] = load_high_accuracy('test_update_P');
            P_up = triu(testCase.P_new);
            Px = nonzeros(P_up);
            testCase.solver.update('Px', Px);
            res = testCase.solver.solve();
            testCase.verifyEqual(res.x, x_ref, 'AbsTol', testCase.tol);
            testCase.verifyEqual(res.y, y_ref, 'AbsTol', testCase.tol);
            testCase.verifyEqual(res.info.obj_val, obj_ref, 'AbsTol', testCase.tol);
        end

        function test_update_P_allind(testCase)
            [x_ref, y_ref, obj_ref] = load_high_accuracy('test_update_P_allind');
            P_up = triu(testCase.P_new);
            [~, ~, Px] = find(P_up);
            [ri, ci] = find(P_up);
            Px_idx = zeros(length(ri), 1);
            Pp = [0; cumsum(full(sum(P_up ~= 0, 1)'))];
            for k = 1:length(ri)
                col = ci(k);
                Px_idx(k) = Pp(col) + sum(ci(1:k-1) == col & ri(1:k-1) < ri(k));
            end
            testCase.solver.update('Px', Px, 'Px_idx', Px_idx + 1);
            res = testCase.solver.solve();
            testCase.verifyEqual(res.x, x_ref, 'AbsTol', testCase.tol);
            testCase.verifyEqual(res.y, y_ref, 'AbsTol', testCase.tol);
            testCase.verifyEqual(res.info.obj_val, obj_ref, 'AbsTol', testCase.tol);
        end

        function test_update_A(testCase)
            [x_ref, y_ref, obj_ref] = load_high_accuracy('test_update_A');
            Ax = nonzeros(testCase.A_new);
            testCase.solver.update('Ax', Ax);
            res = testCase.solver.solve();
            testCase.verifyEqual(res.x, x_ref, 'AbsTol', testCase.tol);
            testCase.verifyEqual(res.y, y_ref, 'AbsTol', testCase.tol);
            testCase.verifyEqual(res.info.obj_val, obj_ref, 'AbsTol', testCase.tol);
        end

        function test_update_A_allind(testCase)
            [x_ref, y_ref, obj_ref] = load_high_accuracy('test_update_A_allind');
            [~, ~, Ax] = find(testCase.A_new);
            Ax_idx = (1:length(Ax))';
            testCase.solver.update('Ax', Ax, 'Ax_idx', Ax_idx);
            res = testCase.solver.solve();
            testCase.verifyEqual(res.x, x_ref, 'AbsTol', testCase.tol);
            testCase.verifyEqual(res.y, y_ref, 'AbsTol', testCase.tol);
            testCase.verifyEqual(res.info.obj_val, obj_ref, 'AbsTol', testCase.tol);
        end

        function test_update_P_A_indP_indA(testCase)
            [x_ref, y_ref, obj_ref] = load_high_accuracy('test_update_P_A_indP_indA');
            P_up = triu(testCase.P_new);
            Px = nonzeros(P_up);
            [~, ~, Ax] = find(testCase.A_new);
            Px_idx = (1:length(Px))';
            Ax_idx = (1:length(Ax))';
            testCase.solver.update('Px', Px, 'Px_idx', Px_idx, ...
                                   'Ax', Ax, 'Ax_idx', Ax_idx);
            res = testCase.solver.solve();
            testCase.verifyEqual(res.x, x_ref, 'AbsTol', testCase.tol);
            testCase.verifyEqual(res.y, y_ref, 'AbsTol', testCase.tol);
            testCase.verifyEqual(res.info.obj_val, obj_ref, 'AbsTol', testCase.tol);
        end

        function test_update_P_A_indP(testCase)
            [x_ref, y_ref, obj_ref] = load_high_accuracy('test_update_P_A_indP');
            P_up = triu(testCase.P_new);
            Px = nonzeros(P_up);
            Px_idx = (1:length(Px))';
            Ax = nonzeros(testCase.A_new);
            testCase.solver.update('Px', Px, 'Px_idx', Px_idx, 'Ax', Ax);
            res = testCase.solver.solve();
            testCase.verifyEqual(res.x, x_ref, 'AbsTol', testCase.tol);
            testCase.verifyEqual(res.y, y_ref, 'AbsTol', testCase.tol);
            testCase.verifyEqual(res.info.obj_val, obj_ref, 'AbsTol', testCase.tol);
        end

        function test_update_P_A_indA(testCase)
            [x_ref, y_ref, obj_ref] = load_high_accuracy('test_update_P_A_indA');
            Px = nonzeros(triu(testCase.P_new));
            [~, ~, Ax] = find(testCase.A_new);
            Ax_idx = (1:length(Ax))';
            testCase.solver.update('Px', Px, 'Ax', Ax, 'Ax_idx', Ax_idx);
            res = testCase.solver.solve();
            testCase.verifyEqual(res.x, x_ref, 'AbsTol', testCase.tol);
            testCase.verifyEqual(res.y, y_ref, 'AbsTol', testCase.tol);
            testCase.verifyEqual(res.info.obj_val, obj_ref, 'AbsTol', testCase.tol);
        end

        function test_update_P_A_allind(testCase)
            [x_ref, y_ref, obj_ref] = load_high_accuracy('test_update_P_A_allind');
            Px = nonzeros(triu(testCase.P_new));
            Ax = nonzeros(testCase.A_new);
            testCase.solver.update('Px', Px, 'Ax', Ax);
            res = testCase.solver.solve();
            testCase.verifyEqual(res.x, x_ref, 'AbsTol', testCase.tol);
            testCase.verifyEqual(res.y, y_ref, 'AbsTol', testCase.tol);
            testCase.verifyEqual(res.info.obj_val, obj_ref, 'AbsTol', testCase.tol);
        end
    end
end
