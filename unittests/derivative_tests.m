classdef derivative_tests < matlab.unittest.TestCase
    % DERIVATIVE_TESTS Test adjoint derivatives for OSQP v1.0.0
    %
    % Tests verify adjoint derivative gradients against finite differences.

    properties
        tol
        abs_tol
        eps_fd  % finite difference step size
    end

    methods(TestMethodSetup)
        function setup_problem(testCase)
            testCase.tol = 5e-3;
            testCase.abs_tol = 5e-3;
            testCase.eps_fd = 1e-5;
        end
    end

    methods(Static)
        function [P, q, A, l, u, true_x, true_y] = get_prob(n, m)
            % Generate a random QP with known feasible point
            rng(1);
            Pt = sprandn(n, n, 0.5);
            P = Pt * Pt' + speye(n);
            q = randn(n, 1);
            A = sprandn(m, n, 0.5);
            true_x = randn(n, 1);
            true_y = randn(m, 1);
            l = A * true_x - rand(m, 1);
            u = A * true_x + rand(m, 1);
        end
    end

    methods(Test)
        function test_dl_dq(testCase)
            n = 10; m = 20;
            [P, q, A, l, u, true_x, true_y] = derivative_tests.get_prob(n, m);

            % Compute gradient via adjoint
            [~, dq, ~, ~] = compute_grads(P, q, A, l, u, true_x, true_y);

            % Finite difference
            dq_fd = zeros(n, 1);
            for i = 1:n
                q_p = q; q_p(i) = q_p(i) + testCase.eps_fd;
                q_m = q; q_m(i) = q_m(i) - testCase.eps_fd;
                loss_p = compute_loss(P, q_p, A, l, u, true_x, true_y);
                loss_m = compute_loss(P, q_m, A, l, u, true_x, true_y);
                dq_fd(i) = (loss_p - loss_m) / (2 * testCase.eps_fd);
            end

            testCase.verifyEqual(dq, dq_fd, 'AbsTol', testCase.abs_tol, ...
                'RelTol', testCase.tol);
        end

        function test_dl_dP(testCase)
            n = 5; m = 10;
            [P, q, A, l, u, true_x, true_y] = derivative_tests.get_prob(n, m);

            % Compute gradient via adjoint
            [dP, ~, ~, ~] = compute_grads(P, q, A, l, u, true_x, true_y);

            % Finite difference over upper-triangular entries
            P_up = triu(P);
            [ri, ci, ~] = find(P_up);
            for k = 1:length(ri)
                ii = ri(k); jj = ci(k);
                P_p = P; P_p(ii, jj) = P_p(ii, jj) + testCase.eps_fd;
                if ii ~= jj, P_p(jj, ii) = P_p(jj, ii) + testCase.eps_fd; end
                P_m = P; P_m(ii, jj) = P_m(ii, jj) - testCase.eps_fd;
                if ii ~= jj, P_m(jj, ii) = P_m(jj, ii) - testCase.eps_fd; end
                loss_p = compute_loss(P_p, q, A, l, u, true_x, true_y);
                loss_m = compute_loss(P_m, q, A, l, u, true_x, true_y);
                fd_val = (loss_p - loss_m) / (2 * testCase.eps_fd);
                adj_val = full(dP(ii, jj));
                testCase.verifyEqual(adj_val, fd_val, 'AbsTol', testCase.abs_tol);
            end
        end

        function test_dl_dA(testCase)
            n = 5; m = 10;
            [P, q, A, l, u, true_x, true_y] = derivative_tests.get_prob(n, m);

            % Compute gradient via adjoint
            [~, ~, dA, ~] = compute_grads(P, q, A, l, u, true_x, true_y);

            % Finite difference over nonzero entries of A
            [ri, ci, ~] = find(A);
            for k = 1:length(ri)
                ii = ri(k); jj = ci(k);
                A_p = A; A_p(ii, jj) = A_p(ii, jj) + testCase.eps_fd;
                A_m = A; A_m(ii, jj) = A_m(ii, jj) - testCase.eps_fd;
                loss_p = compute_loss(P, q, A_p, l, u, true_x, true_y);
                loss_m = compute_loss(P, q, A_m, l, u, true_x, true_y);
                fd_val = (loss_p - loss_m) / (2 * testCase.eps_fd);
                adj_val = full(dA(ii, jj));
                testCase.verifyEqual(adj_val, fd_val, 'AbsTol', testCase.abs_tol);
            end
        end

        function test_dl_dl(testCase)
            n = 30; m = 30;
            [P, q, A, l, u, true_x, true_y] = derivative_tests.get_prob(n, m);

            [~, ~, ~, db] = compute_grads(P, q, A, l, u, true_x, true_y);
            dl = db.dl;

            dl_fd = zeros(m, 1);
            for i = 1:m
                l_p = l; l_p(i) = l_p(i) + testCase.eps_fd;
                l_m = l; l_m(i) = l_m(i) - testCase.eps_fd;
                loss_p = compute_loss(P, q, A, l_p, u, true_x, true_y);
                loss_m = compute_loss(P, q, A, l_m, u, true_x, true_y);
                dl_fd(i) = (loss_p - loss_m) / (2 * testCase.eps_fd);
            end

            testCase.verifyEqual(dl, dl_fd, 'AbsTol', testCase.abs_tol, ...
                'RelTol', testCase.tol);
        end

        function test_dl_du(testCase)
            n = 10; m = 20;
            [P, q, A, l, u, true_x, true_y] = derivative_tests.get_prob(n, m);

            [~, ~, ~, db] = compute_grads(P, q, A, l, u, true_x, true_y);
            du = db.du;

            du_fd = zeros(m, 1);
            for i = 1:m
                u_p = u; u_p(i) = u_p(i) + testCase.eps_fd;
                u_m = u; u_m(i) = u_m(i) - testCase.eps_fd;
                loss_p = compute_loss(P, q, A, l, u_p, true_x, true_y);
                loss_m = compute_loss(P, q, A, l, u_m, true_x, true_y);
                du_fd(i) = (loss_p - loss_m) / (2 * testCase.eps_fd);
            end

            testCase.verifyEqual(du, du_fd, 'AbsTol', testCase.abs_tol, ...
                'RelTol', testCase.tol);
        end
    end
end


%% ---- Local helper functions ----

function loss = compute_loss(P, q, A, l, u, true_x, true_y)
    % Solve QP and compute loss = 0.5*||x-true_x||^2 + 0.5*||y-true_y||^2
    solver = osqp;
    solver.setup(P, q, A, l, u, ...
        'verbose', false, ...
        'eps_abs', 1e-9, 'eps_rel', 1e-9, ...
        'max_iter', 500000, 'polishing', true);
    res = solver.solve();
    loss = 0.5 * sum((res.x - true_x).^2) + 0.5 * sum((res.y - true_y).^2);
    delete(solver);
end


function [dP, dq, dA, db] = compute_grads(P, q, A, l, u, true_x, true_y)
    % Solve QP and compute adjoint derivatives
    solver = osqp;
    solver.setup(P, q, A, l, u, ...
        'verbose', false, ...
        'eps_abs', 1e-9, 'eps_rel', 1e-9, ...
        'max_iter', 500000, 'polishing', true);
    res = solver.solve();

    % dloss/dx = x - true_x, dloss/dy = y - true_y
    dx = res.x - true_x;
    dy = res.y - true_y;

    solver.adjoint_derivative_compute(dx, dy);
    [dP, dA] = solver.adjoint_derivative_get_mat();
    [dq, dl, du] = solver.adjoint_derivative_get_vec();
    db = struct('dl', dl, 'du', du);

    delete(solver);
end
