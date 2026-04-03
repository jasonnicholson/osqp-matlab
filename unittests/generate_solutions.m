function generate_solutions()
% GENERATE_SOLUTIONS  Regenerate all reference solution .mat files
%   using the current OSQP v1.0.0 solver.
%
%   Run this from the osqp-matlab root directory.

sol_dir = fullfile(fileparts(mfilename('fullpath')), 'solutions');
if ~exist(sol_dir, 'dir')
    mkdir(sol_dir);
end

%% ---- basic_tests data ----
P = sparse([11 0; 0 0]);
q = [3; 4];
A = sparse([-1 0; 0 -1; -1 -3; 2 5; 3 4]);
u = [0; 0; -15; 100; 80];
l = -1e30 * ones(5, 1);
n = 2;
m = 5;

opts = {'verbose', false, 'polishing', true, 'eps_abs', 1e-8, 'eps_rel', 1e-8};

% test_basic_QP
solver = osqp;
solver.setup(P, q, A, l, u, opts{:});
res = solver.solve();
save_sol(sol_dir, 'test_basic_QP', res);
delete(solver);

% test_update_q
solver = osqp;
solver.setup(P, q, A, l, u, opts{:});
solver.update('q', [10; 20]);
res = solver.solve();
save_sol(sol_dir, 'test_update_q', res);
delete(solver);

% test_update_l
solver = osqp;
solver.setup(P, q, A, l, u, opts{:});
solver.update('l', -50 * ones(m, 1));
res = solver.solve();
save_sol(sol_dir, 'test_update_l', res);
delete(solver);

% test_update_u
solver = osqp;
solver.setup(P, q, A, l, u, opts{:});
solver.update('u', 1000 * ones(m, 1));
res = solver.solve();
save_sol(sol_dir, 'test_update_u', res);
delete(solver);

% test_update_bounds
solver = osqp;
solver.setup(P, q, A, l, u, opts{:});
solver.update('l', -50 * ones(m, 1), 'u', 1000 * ones(m, 1));
res = solver.solve();
save_sol(sol_dir, 'test_update_bounds', res);
delete(solver);

%% ---- update_matrices_tests data ----
rng(1);
n = 5; m = 8; p = 0.7;
Pt = sprandn(n, n, p);
P2 = Pt * Pt' + speye(n);
q2 = randn(n, 1);
A2 = sprandn(m, n, p);
l2 = -2 * ones(m, 1);
u2 =  2 * ones(m, 1);

% Generate P_new and A_new with SAME sparsity patterns (mirrors test setup)
P2triu = triu(P2);
[pi, pj, ~] = find(P2triu);
P_new_vals = randn(nnz(P2triu), 1);
P2triu_new = sparse(pi, pj, P_new_vals, n, n);
P2_new = P2triu_new + triu(P2triu_new, 1)' + 5*speye(n);

[ai, aj, ~] = find(A2);
A_new_vals = randn(nnz(A2), 1);
A2_new = sparse(ai, aj, A_new_vals, m, n);

% test_solve
solver = osqp;
solver.setup(P2, q2, A2, l2, u2, opts{:});
res = solver.solve();
save_sol(sol_dir, 'test_solve', res);
delete(solver);

% test_update_P
solver = osqp;
solver.setup(P2, q2, A2, l2, u2, opts{:});
Px = nonzeros(triu(P2_new));
solver.update('Px', Px);
res = solver.solve();
save_sol(sol_dir, 'test_update_P', res);
delete(solver);

% test_update_P_allind
solver = osqp;
solver.setup(P2, q2, A2, l2, u2, opts{:});
P_up = triu(P2_new);
[~, ~, Px] = find(P_up);
Px_idx = (1:length(Px))';
solver.update('Px', Px, 'Px_idx', Px_idx);
res = solver.solve();
save_sol(sol_dir, 'test_update_P_allind', res);
delete(solver);

% test_update_A
solver = osqp;
solver.setup(P2, q2, A2, l2, u2, opts{:});
Ax = nonzeros(A2_new);
solver.update('Ax', Ax);
res = solver.solve();
save_sol(sol_dir, 'test_update_A', res);
delete(solver);

% test_update_A_allind
solver = osqp;
solver.setup(P2, q2, A2, l2, u2, opts{:});
[~, ~, Ax] = find(A2_new);
Ax_idx = (1:length(Ax))';
solver.update('Ax', Ax, 'Ax_idx', Ax_idx);
res = solver.solve();
save_sol(sol_dir, 'test_update_A_allind', res);
delete(solver);

% test_update_P_A_indP_indA
solver = osqp;
solver.setup(P2, q2, A2, l2, u2, opts{:});
P_up = triu(P2_new);
Px = nonzeros(P_up);
[~, ~, Ax] = find(A2_new);
Px_idx = (1:length(Px))';
Ax_idx = (1:length(Ax))';
solver.update('Px', Px, 'Px_idx', Px_idx, 'Ax', Ax, 'Ax_idx', Ax_idx);
res = solver.solve();
save_sol(sol_dir, 'test_update_P_A_indP_indA', res);
delete(solver);

% test_update_P_A_indP
solver = osqp;
solver.setup(P2, q2, A2, l2, u2, opts{:});
solver.update('Px', nonzeros(triu(P2_new)), 'Px_idx', (1:nnz(triu(P2_new)))', ...
              'Ax', nonzeros(A2_new));
res = solver.solve();
save_sol(sol_dir, 'test_update_P_A_indP', res);
delete(solver);

% test_update_P_A_indA
solver = osqp;
solver.setup(P2, q2, A2, l2, u2, opts{:});
Ax = nonzeros(A2_new);
Ax_idx = (1:length(Ax))';
solver.update('Px', nonzeros(triu(P2_new)), 'Ax', Ax, 'Ax_idx', Ax_idx);
res = solver.solve();
save_sol(sol_dir, 'test_update_P_A_indA', res);
delete(solver);

% test_update_P_A_allind
solver = osqp;
solver.setup(P2, q2, A2, l2, u2, opts{:});
solver.update('Px', nonzeros(triu(P2_new)), 'Ax', nonzeros(A2_new));
res = solver.solve();
save_sol(sol_dir, 'test_update_P_A_allind', res);
delete(solver);

%% ---- primal_infeasibility data ----
rng(4);
n = 50; m = 500;
Pt = sprandn(n, n, 0.6);
P_pi = Pt' * Pt;
q_pi = randn(n, 1);
A_pi = sprandn(m, n, 0.8);
u_pi = 3 + randn(m, 1);
l_pi = -3 + randn(m, 1);
nhalf = floor(n/2);
A_pi(nhalf, :) = A_pi(nhalf + 1, :);
l_pi(nhalf) = u_pi(nhalf + 1) + 10 * rand();
u_pi(nhalf) = l_pi(nhalf) + 0.5;

solver = osqp;
solver.setup(P_pi, q_pi, A_pi, l_pi, u_pi, ...
    'verbose', false, 'eps_prim_inf', 1e-6, ...
    'polishing', true, 'eps_abs', 1e-8, 'eps_rel', 1e-8);
res = solver.solve();
save_sol(sol_dir, 'test_primal_infeasibility', res);
delete(solver);

%% ---- dual_infeasibility data ----
P_di = sparse(diag([4; 0]));
q_di = [0; 2];
A_di = sparse([1 1; -1 1]);
l_di = [-Inf; -Inf];
u_di = [2; 3];

solver = osqp;
solver.setup(P_di, q_di, A_di, l_di, u_di, ...
    'verbose', false, 'polishing', true, ...
    'eps_abs', 1e-8, 'eps_rel', 1e-8);
res = solver.solve();
save_sol(sol_dir, 'test_dual_infeasibility', res);
delete(solver);

%% ---- feasibility_tests data ----
rng(4);
n = 30; m = 30;
P_f = sparse(n, n);
q_f = zeros(n, 1);
A_f = sprandn(m, n, 0.8);
b_f = randn(m, 1);
l_f = b_f;
u_f = b_f;

solver = osqp;
solver.setup(P_f, q_f, A_f, l_f, u_f, ...
    'verbose', false, 'polishing', true, ...
    'eps_abs', 1e-8, 'eps_rel', 1e-8);
res = solver.solve();
save_sol(sol_dir, 'test_feasibility_problem', res);
delete(solver);

%% ---- unconstrained data ----
rng(4);
n = 30;
Pt = sprandn(n, n, 0.6);
P_uc = Pt * Pt' + speye(n);
q_uc = randn(n, 1);
A_uc = sparse(0, n);
l_uc = zeros(0, 1);
u_uc = zeros(0, 1);

solver = osqp;
solver.setup(P_uc, q_uc, A_uc, l_uc, u_uc, ...
    'verbose', false, 'polishing', true, ...
    'eps_abs', 1e-8, 'eps_rel', 1e-8);
res = solver.solve();
save_sol(sol_dir, 'test_unconstrained_problem', res);
delete(solver);

%% ---- polishing data ----
% test_polish_simple
P_ps = sparse([11 0; 0 0]);
q_ps = [3; 4];
A_ps = sparse([-1 0; 0 -1; -1 -3; 2 5; 3 4]);
u_ps = [0; 0; -15; 100; 80];
l_ps = -1e30 * ones(5, 1);

solver = osqp;
solver.setup(P_ps, q_ps, A_ps, l_ps, u_ps, ...
    'verbose', false, 'polishing', true, ...
    'eps_abs', 1e-8, 'eps_rel', 1e-8);
res = solver.solve();
save_sol(sol_dir, 'test_polish_simple', res);
delete(solver);

% test_polish_unconstrained
rng(4);
n = 30;
Pt = sprandn(n, n, 0.6);
P_pu = Pt * Pt' + speye(n);
q_pu = randn(n, 1);
A_pu = sparse(0, n);
l_pu = zeros(0, 1);
u_pu = zeros(0, 1);

solver = osqp;
solver.setup(P_pu, q_pu, A_pu, l_pu, u_pu, ...
    'verbose', false, 'polishing', true, ...
    'polish_refine_iter', 4, ...
    'eps_abs', 1e-8, 'eps_rel', 1e-8);
res = solver.solve();
save_sol(sol_dir, 'test_polish_unconstrained', res);
delete(solver);

% test_polish_random
rng(6);
n = 30; m = 50;
Pt = sprandn(n, n, 0.6);
P_pr = Pt * Pt' + speye(n);
q_pr = randn(n, 1);
A_pr = sprandn(m, n, 0.8);
l_pr = -2 * ones(m, 1);
u_pr =  2 * ones(m, 1);

solver = osqp;
solver.setup(P_pr, q_pr, A_pr, l_pr, u_pr, ...
    'verbose', false, 'polishing', true, ...
    'polish_refine_iter', 4, ...
    'eps_abs', 1e-8, 'eps_rel', 1e-8);
res = solver.solve();
save_sol(sol_dir, 'test_polish_random', res);
delete(solver);

fprintf('All reference solutions generated in %s\n', sol_dir);
end

function save_sol(sol_dir, name, res)
    x_val = res.x; %#ok<NASGU>
    y_val = res.y; %#ok<NASGU>
    obj   = res.info.obj_val; %#ok<NASGU>
    save(fullfile(sol_dir, [name '.mat']), 'x_val', 'y_val', 'obj');
    fprintf('  Saved %s.mat\n', name);
end
