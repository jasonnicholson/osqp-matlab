% run_ml_unit_tests.m — Faithfully replicate the unit tests using matlab_ldl backend

fprintf('=============================================================\n');
fprintf('  OSQP MATLAB Backend (matlab_ldl) Unit Tests\n');
fprintf('=============================================================\n\n');

tol = 1e-4;
npass = 0;
nfail = 0;

ML = {'linear_solver', 'matlab_ldl'};

% Force MATLAB backend via factory function
newSolver = @() osqp('backend', 'matlab');

%% ===== BASIC TESTS (same problem as basic_tests.m) =====
fprintf('--- basic_tests ---\n');

P = sparse([11 0; 0 0]);
q = [3; 4];
A = sparse([-1 0; 0 -1; -1 -3; 2 5; 3 4]);
u = [0; 0; -15; 100; 80];
l = -1e30 * ones(5, 1);
n = 2; m = 5;

basic_opts = {'verbose', false, 'polishing', true, 'eps_abs', 1e-8, 'eps_rel', 1e-8, ML{:}};

tests = {'test_basic_QP', 'test_update_q', 'test_update_l', 'test_update_u', 'test_update_bounds'};
updates = {
    {}
    {'q', [10; 20]}
    {'l', -50*ones(m,1)}
    {'u', 1000*ones(m,1)}
    {'l', -50*ones(m,1), 'u', 1000*ones(m,1)}
};

for t = 1:numel(tests)
    try
        solver = newSolver(); solver.setup(P, q, A, l, u, basic_opts{:});
        if ~isempty(updates{t})
            solver.update(updates{t}{:});
        end
        [x_ref, y_ref, obj_ref] = load_high_accuracy(tests{t});
        res = solver.solve();
        assert(max(abs(res.x - x_ref)) < tol);
        assert(max(abs(res.y - y_ref)) < tol);
        assert(abs(res.info.obj_val - obj_ref) < tol);
        fprintf('  %s: PASS\n', tests{t}); npass = npass + 1;
    catch e
        fprintf('  %s: FAIL (%s)\n', tests{t}, e.message); nfail = nfail + 1;
    end
    delete(solver);
end

% test_update_max_iter
try
    solver = newSolver(); solver.setup(P, q, A, l, u, basic_opts{:});
    solver.update_settings('max_iter', 80);
    res = solver.solve();
    assert(res.info.status_val == osqp.constant.OSQP_MAX_ITER_REACHED);
    fprintf('  test_max_iter: PASS\n'); npass = npass + 1;
catch e
    fprintf('  test_max_iter: FAIL (%s)\n', e.message); nfail = nfail + 1;
end
delete(solver);

% test_early_term
try
    solver = newSolver(); solver.setup(P, q, A, l, u, basic_opts{:});
    solver.update_settings('check_termination', 0, 'max_iter', 500);
    res = solver.solve();
    assert(res.info.iter == 500);
    fprintf('  test_early_term: PASS\n'); npass = npass + 1;
catch e
    fprintf('  test_early_term: FAIL (%s)\n', e.message); nfail = nfail + 1;
end
delete(solver);

%% ===== WARM START TESTS (rng(4), n=100, m=200) =====
fprintf('--- warm_start_tests ---\n');

try
    rng(4);
    n_ws = 100; m_ws = 200;
    Pt = sprandn(n_ws, n_ws, 0.6);
    P_ws = Pt * Pt' + speye(n_ws);
    q_ws = randn(n_ws, 1);
    A_ws = sprandn(m_ws, n_ws, 0.8);
    l_ws = -2 * ones(m_ws, 1);
    u_ws =  2 * ones(m_ws, 1);

    solver = newSolver();
    solver.setup(P_ws, q_ws, A_ws, l_ws, u_ws, 'verbose', false, ...
        'eps_abs', 1e-8, 'eps_rel', 1e-8, 'warm_starting', true, ML{:});

    solver.warm_start('x', zeros(n_ws, 1), 'y', zeros(m_ws, 1));
    res_cold = solver.solve();

    solver.warm_start('x', res_cold.x, 'y', res_cold.y);
    res_warm = solver.solve();

    assert(res_warm.info.iter < res_cold.info.iter);
    fprintf('  test_warm_start (cold=%d warm=%d): PASS\n', res_cold.info.iter, res_warm.info.iter);
    npass = npass + 1;
catch e
    fprintf('  test_warm_start: FAIL (%s)\n', e.message); nfail = nfail + 1;
end
delete(solver);

%% ===== POLISHING TESTS =====
fprintf('--- polishing_tests ---\n');

% test_polish_simple (same basic problem)
try
    solver = newSolver();
    solver.setup(sparse([11 0; 0 0]), [3; 4], sparse([-1 0; 0 -1; -1 -3; 2 5; 3 4]), ...
        -1e30*ones(5,1), [0; 0; -15; 100; 80], 'verbose', false, 'polishing', true, ...
        'polish_refine_iter', 4, 'eps_abs', 1e-8, 'eps_rel', 1e-8, ML{:});
    [x_ref, y_ref, obj_ref] = load_high_accuracy('test_polish_simple');
    res = solver.solve();
    assert(max(abs(res.x - x_ref)) < tol);
    assert(max(abs(res.y - y_ref)) < tol);
    assert(abs(res.info.obj_val - obj_ref) < tol);
    fprintf('  test_polish_simple: PASS\n'); npass = npass + 1;
catch e
    fprintf('  test_polish_simple: FAIL (%s)\n', e.message); nfail = nfail + 1;
end
delete(solver);

% test_polish_unconstrained
try
    rng(4);
    n_u = 30;
    Pt = sprandn(n_u, n_u, 0.6);
    P_u = Pt * Pt' + speye(n_u);
    q_u = randn(n_u, 1);
    A_u = sparse(0, n_u);
    l_u = zeros(0, 1);
    u_u = zeros(0, 1);

    [x_ref, ~, obj_ref] = load_high_accuracy('test_polish_unconstrained');
    solver = newSolver();
    solver.setup(P_u, q_u, A_u, l_u, u_u, 'verbose', false, 'polishing', true, ...
        'polish_refine_iter', 4, 'eps_abs', 1e-8, 'eps_rel', 1e-8, ML{:});
    res = solver.solve();
    assert(max(abs(res.x - x_ref)) < tol);
    assert(abs(res.info.obj_val - obj_ref) < tol);
    fprintf('  test_polish_unconstrained: PASS\n'); npass = npass + 1;
catch e
    fprintf('  test_polish_unconstrained: FAIL (%s)\n', e.message); nfail = nfail + 1;
end
delete(solver);

% test_polish_random
try
    rng(6);
    n_r = 30; m_r = 50;
    Pt = sprandn(n_r, n_r, 0.6);
    P_r = Pt * Pt' + speye(n_r);
    q_r = randn(n_r, 1);
    A_r = sprandn(m_r, n_r, 0.8);
    l_r = -2 * ones(m_r, 1);
    u_r =  2 * ones(m_r, 1);

    [x_ref, y_ref, obj_ref] = load_high_accuracy('test_polish_random');
    solver = newSolver();
    solver.setup(P_r, q_r, A_r, l_r, u_r, 'verbose', false, 'polishing', true, ...
        'polish_refine_iter', 4, 'eps_abs', 1e-8, 'eps_rel', 1e-8, ML{:});
    res = solver.solve();
    assert(max(abs(res.x - x_ref)) < tol);
    assert(max(abs(res.y - y_ref)) < tol);
    assert(abs(res.info.obj_val - obj_ref) < tol);
    fprintf('  test_polish_random: PASS\n'); npass = npass + 1;
catch e
    fprintf('  test_polish_random: FAIL (%s)\n', e.message); nfail = nfail + 1;
end
delete(solver);

%% ===== FEASIBILITY TESTS =====
fprintf('--- feasibility_tests ---\n');

try
    rng(4);
    n_f = 30; m_f = 30;
    Pf = sparse(n_f, n_f);
    qf = zeros(n_f, 1);
    Af = sprandn(m_f, n_f, 0.8);
    bf = randn(m_f, 1);
    lf = bf;
    uf = bf;

    [x_ref, y_ref, ~] = load_high_accuracy('test_feasibility_problem');
    solver = newSolver();
    solver.setup(Pf, qf, Af, lf, uf, 'verbose', false, 'polishing', true, ...
        'eps_abs', 1e-8, 'eps_rel', 1e-8, ML{:});
    res = solver.solve();
    assert(max(abs(res.x - x_ref)) < tol);
    assert(max(abs(res.y - y_ref)) < tol);
    fprintf('  test_feasibility_problem: PASS\n'); npass = npass + 1;
catch e
    fprintf('  test_feasibility_problem: FAIL (%s)\n', e.message); nfail = nfail + 1;
end
delete(solver);

%% ===== PRIMAL INFEASIBILITY TESTS =====
fprintf('--- primal_infeasibility_tests ---\n');

try
    rng(3);
    n_pi = 50; m_pi = 500;
    Pt = sprandn(n_pi, n_pi, 0.6);
    P_pi = Pt * Pt' + speye(n_pi);
    q_pi = randn(n_pi, 1);
    A_pi = sprandn(m_pi, n_pi, 0.8);
    u_pi = randn(m_pi, 1);
    l_pi = u_pi + 100;  % makes it infeasible (l > u)

    solver = newSolver();
    solver.setup(P_pi, q_pi, A_pi, l_pi, u_pi, 'verbose', false, ...
        'eps_prim_inf', 1e-4, ML{:});
    res = solver.solve();
    assert(res.info.status_val == osqp.constant.OSQP_PRIMAL_INFEASIBLE);
    fprintf('  test_primal_infeasibility: PASS\n'); npass = npass + 1;
catch e
    fprintf('  test_primal_infeasibility: FAIL (%s)\n', e.message); nfail = nfail + 1;
end
delete(solver);

%% ===== DUAL INFEASIBILITY TESTS =====
fprintf('--- dual_infeasibility_tests ---\n');

try
    n_di = 2;
    P_di = sparse(n_di, n_di);
    q_di = [-1; -1];
    A_di = speye(n_di);
    l_di = zeros(n_di, 1);
    u_di = inf * ones(n_di, 1);

    solver = newSolver();
    solver.setup(P_di, q_di, A_di, l_di, u_di, 'verbose', false, ...
        'eps_dual_inf', 1e-4, ML{:});
    res = solver.solve();
    assert(res.info.status_val == osqp.constant.OSQP_DUAL_INFEASIBLE);
    fprintf('  test_dual_infeasibility: PASS\n'); npass = npass + 1;
catch e
    fprintf('  test_dual_infeasibility: FAIL (%s)\n', e.message); nfail = nfail + 1;
end
delete(solver);

%% ===== UNCONSTRAINED TESTS =====
fprintf('--- unconstrained_tests ---\n');

try
    rng(4);
    n_uc = 30;
    Pt = sprandn(n_uc, n_uc, 0.6);
    P_uc = Pt * Pt' + speye(n_uc);
    q_uc = randn(n_uc, 1);
    A_uc = sparse(0, n_uc);
    l_uc = zeros(0, 1);
    u_uc = zeros(0, 1);

    solver = newSolver();
    solver.setup(P_uc, q_uc, A_uc, l_uc, u_uc, 'verbose', false, ...
        'eps_abs', 1e-8, 'eps_rel', 1e-8, ML{:});
    res = solver.solve();

    % Unconstrained QP: x = -P\q
    x_expected = -(P_uc \ q_uc);
    assert(max(abs(res.x - x_expected)) < tol);
    fprintf('  test_unconstrained: PASS\n'); npass = npass + 1;
catch e
    fprintf('  test_unconstrained: FAIL (%s)\n', e.message); nfail = nfail + 1;
end
delete(solver);

%% ===== UPDATE MATRICES TESTS (rng(1), n=5, m=8) =====
fprintf('--- update_matrices_tests ---\n');

rng(1);
n_um = 5; m_um = 8; p_um = 0.7;
Pt = sprandn(n_um, n_um, p_um);
P_um = Pt * Pt' + speye(n_um);
q_um = randn(n_um, 1);
A_um = sprandn(m_um, n_um, p_um);
l_um = -2 * ones(m_um, 1);
u_um =  2 * ones(m_um, 1);

P_triu_um = triu(P_um);
[pi2, pj2, ~] = find(P_triu_um);
P_new_vals = randn(nnz(P_triu_um), 1);
P_triu_new = sparse(pi2, pj2, P_new_vals, n_um, n_um);
P_new_um = P_triu_new + triu(P_triu_new, 1)' + 5*speye(n_um);

[ai2, aj2, ~] = find(A_um);
A_new_vals = randn(nnz(A_um), 1);
A_new_um = sparse(ai2, aj2, A_new_vals, m_um, n_um);

um_opts = {'verbose', false, 'polishing', true, 'eps_abs', 1e-8, 'eps_rel', 1e-8, ML{:}};

% test_solve (original problem)
try
    solver = newSolver(); solver.setup(P_um, q_um, A_um, l_um, u_um, um_opts{:});
    [x_ref, y_ref, obj_ref] = load_high_accuracy('test_solve');
    res = solver.solve();
    assert(max(abs(res.x - x_ref)) < tol);
    assert(max(abs(res.y - y_ref)) < tol);
    assert(abs(res.info.obj_val - obj_ref) < tol);
    fprintf('  test_solve: PASS\n'); npass = npass + 1;
catch e
    fprintf('  test_solve: FAIL (%s)\n', e.message); nfail = nfail + 1;
end
delete(solver);

% test_update_P
try
    solver = newSolver(); solver.setup(P_um, q_um, A_um, l_um, u_um, um_opts{:});
    P_up = triu(P_new_um);
    Px = nonzeros(P_up);
    solver.update('Px', Px);
    [x_ref, y_ref, obj_ref] = load_high_accuracy('test_update_P');
    res = solver.solve();
    assert(max(abs(res.x - x_ref)) < tol);
    assert(max(abs(res.y - y_ref)) < tol);
    assert(abs(res.info.obj_val - obj_ref) < tol);
    fprintf('  test_update_P: PASS\n'); npass = npass + 1;
catch e
    fprintf('  test_update_P: FAIL (%s)\n', e.message); nfail = nfail + 1;
end
delete(solver);

% test_update_P_allind
try
    solver = newSolver(); solver.setup(P_um, q_um, A_um, l_um, u_um, um_opts{:});
    P_up = triu(P_new_um);
    [~, ~, Px] = find(P_up);
    Px_idx = (1:length(Px))';
    solver.update('Px', Px, 'Px_idx', Px_idx);
    [x_ref, y_ref, obj_ref] = load_high_accuracy('test_update_P_allind');
    res = solver.solve();
    assert(max(abs(res.x - x_ref)) < tol);
    assert(max(abs(res.y - y_ref)) < tol);
    assert(abs(res.info.obj_val - obj_ref) < tol);
    fprintf('  test_update_P_allind: PASS\n'); npass = npass + 1;
catch e
    fprintf('  test_update_P_allind: FAIL (%s)\n', e.message); nfail = nfail + 1;
end
delete(solver);

% test_update_A
try
    solver = newSolver(); solver.setup(P_um, q_um, A_um, l_um, u_um, um_opts{:});
    Ax = nonzeros(A_new_um);
    solver.update('Ax', Ax);
    [x_ref, y_ref, obj_ref] = load_high_accuracy('test_update_A');
    res = solver.solve();
    assert(max(abs(res.x - x_ref)) < tol);
    assert(max(abs(res.y - y_ref)) < tol);
    assert(abs(res.info.obj_val - obj_ref) < tol);
    fprintf('  test_update_A: PASS\n'); npass = npass + 1;
catch e
    fprintf('  test_update_A: FAIL (%s)\n', e.message); nfail = nfail + 1;
end
delete(solver);

% test_update_A_allind
try
    solver = newSolver(); solver.setup(P_um, q_um, A_um, l_um, u_um, um_opts{:});
    [~, ~, Ax] = find(A_new_um);
    Ax_idx = (1:length(Ax))';
    solver.update('Ax', Ax, 'Ax_idx', Ax_idx);
    [x_ref, y_ref, obj_ref] = load_high_accuracy('test_update_A_allind');
    res = solver.solve();
    assert(max(abs(res.x - x_ref)) < tol);
    assert(max(abs(res.y - y_ref)) < tol);
    assert(abs(res.info.obj_val - obj_ref) < tol);
    fprintf('  test_update_A_allind: PASS\n'); npass = npass + 1;
catch e
    fprintf('  test_update_A_allind: FAIL (%s)\n', e.message); nfail = nfail + 1;
end
delete(solver);

% test_update_P_A_allind
try
    solver = newSolver(); solver.setup(P_um, q_um, A_um, l_um, u_um, um_opts{:});
    Px = nonzeros(triu(P_new_um));
    Ax = nonzeros(A_new_um);
    solver.update('Px', Px, 'Ax', Ax);
    [x_ref, y_ref, obj_ref] = load_high_accuracy('test_update_P_A_allind');
    res = solver.solve();
    assert(max(abs(res.x - x_ref)) < tol);
    assert(max(abs(res.y - y_ref)) < tol);
    assert(abs(res.info.obj_val - obj_ref) < tol);
    fprintf('  test_update_P_A_allind: PASS\n'); npass = npass + 1;
catch e
    fprintf('  test_update_P_A_allind: FAIL (%s)\n', e.message); nfail = nfail + 1;
end
delete(solver);

% test_update_P_A_indP
try
    solver = newSolver(); solver.setup(P_um, q_um, A_um, l_um, u_um, um_opts{:});
    P_up = triu(P_new_um);
    Px = nonzeros(P_up);
    Px_idx = (1:length(Px))';
    Ax = nonzeros(A_new_um);
    solver.update('Px', Px, 'Px_idx', Px_idx, 'Ax', Ax);
    [x_ref, y_ref, obj_ref] = load_high_accuracy('test_update_P_A_indP');
    res = solver.solve();
    assert(max(abs(res.x - x_ref)) < tol);
    assert(max(abs(res.y - y_ref)) < tol);
    assert(abs(res.info.obj_val - obj_ref) < tol);
    fprintf('  test_update_P_A_indP: PASS\n'); npass = npass + 1;
catch e
    fprintf('  test_update_P_A_indP: FAIL (%s)\n', e.message); nfail = nfail + 1;
end
delete(solver);

% test_update_P_A_indA
try
    solver = newSolver(); solver.setup(P_um, q_um, A_um, l_um, u_um, um_opts{:});
    Px = nonzeros(triu(P_new_um));
    [~, ~, Ax] = find(A_new_um);
    Ax_idx = (1:length(Ax))';
    solver.update('Px', Px, 'Ax', Ax, 'Ax_idx', Ax_idx);
    [x_ref, y_ref, obj_ref] = load_high_accuracy('test_update_P_A_indA');
    res = solver.solve();
    assert(max(abs(res.x - x_ref)) < tol);
    assert(max(abs(res.y - y_ref)) < tol);
    assert(abs(res.info.obj_val - obj_ref) < tol);
    fprintf('  test_update_P_A_indA: PASS\n'); npass = npass + 1;
catch e
    fprintf('  test_update_P_A_indA: FAIL (%s)\n', e.message); nfail = nfail + 1;
end
delete(solver);

% test_update_P_A_indP_indA
try
    solver = newSolver(); solver.setup(P_um, q_um, A_um, l_um, u_um, um_opts{:});
    P_up = triu(P_new_um);
    Px = nonzeros(P_up);
    [~, ~, Ax] = find(A_new_um);
    Px_idx = (1:length(Px))';
    Ax_idx = (1:length(Ax))';
    solver.update('Px', Px, 'Px_idx', Px_idx, 'Ax', Ax, 'Ax_idx', Ax_idx);
    [x_ref, y_ref, obj_ref] = load_high_accuracy('test_update_P_A_indP_indA');
    res = solver.solve();
    assert(max(abs(res.x - x_ref)) < tol);
    assert(max(abs(res.y - y_ref)) < tol);
    assert(abs(res.info.obj_val - obj_ref) < tol);
    fprintf('  test_update_P_A_indP_indA: PASS\n'); npass = npass + 1;
catch e
    fprintf('  test_update_P_A_indP_indA: FAIL (%s)\n', e.message); nfail = nfail + 1;
end
delete(solver);

%% ===== NON-CONVEX TESTS =====
fprintf('--- non_cvx_tests ---\n');

try
    P_lp = sparse(2, 2);
    q_lp = [1; -1];
    A_lp = sparse([1 0; 0 1; 1 1]);
    l_lp = [0; 0; -inf];
    u_lp = [inf; inf; 1];
    solver = newSolver();
    solver.setup(P_lp, q_lp, A_lp, l_lp, u_lp, 'verbose', false, ML{:});
    res = solver.solve();
    assert(res.info.status_val == osqp.constant.OSQP_SOLVED || ...
           res.info.status_val == osqp.constant.OSQP_SOLVED_INACCURATE);
    fprintf('  test_lp: PASS\n'); npass = npass + 1;
catch e
    fprintf('  test_lp: FAIL (%s)\n', e.message); nfail = nfail + 1;
end
delete(solver);

%% ===== COLD START TEST =====
fprintf('--- cold_start_tests ---\n');

try
    solver = newSolver();
    solver.setup(sparse([11 0; 0 0]), [3; 4], sparse([-1 0; 0 -1; -1 -3; 2 5; 3 4]), ...
        -1e30*ones(5,1), [0; 0; -15; 100; 80], 'verbose', false, ...
        'eps_abs', 1e-8, 'eps_rel', 1e-8, ML{:});
    res1 = solver.solve();
    solver.cold_start();
    res2 = solver.solve();
    assert(max(abs(res1.x - res2.x)) < tol);
    fprintf('  test_cold_start: PASS\n'); npass = npass + 1;
catch e
    fprintf('  test_cold_start: FAIL (%s)\n', e.message); nfail = nfail + 1;
end
delete(solver);

%% ===== SUMMARY =====
fprintf('\n=============================================================\n');
fprintf('  TOTAL: %d  |  PASS: %d  |  FAIL: %d\n', npass + nfail, npass, nfail);
fprintf('=============================================================\n');
