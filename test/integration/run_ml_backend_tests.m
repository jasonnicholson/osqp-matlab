% run_ml_backend_tests.m — Quick smoke tests for the matlab_ldl backend

P = sparse([11 0; 0 0]);
q = [3; 4];
A = sparse([-1 0; 0 -1; -1 -3; 2 5; 3 4]);
u = [0; 0; -15; 100; 80];
l = -1e30 * ones(5, 1);
tol = 1e-4;
nfail = 0;

%% test_basic_qp
solver = osqp(backend='matlab');
solver.setup(P, q, A, l, u, 'verbose', false, 'polishing', true, ...
    'eps_abs', 1e-8, 'eps_rel', 1e-8, 'linear_solver', 'matlab_ldl');
[x_ref, y_ref, obj_ref] = load_high_accuracy('test_basic_QP');
res = solver.solve();
xerr = max(abs(res.x - x_ref));
yerr = max(abs(res.y - y_ref));
oerr = abs(res.info.obj_val - obj_ref);
fprintf('test_basic_qp:  x_err=%.2e  y_err=%.2e  obj_err=%.2e  %s\n', ...
    xerr, yerr, oerr, pass_fail(xerr < tol && yerr < tol && oerr < tol));
if xerr >= tol || yerr >= tol || oerr >= tol, nfail = nfail + 1; end

%% test_update_q
solver.update('q', [10; 20]);
[x_ref, y_ref, obj_ref] = load_high_accuracy('test_update_q');
res = solver.solve();
xerr = max(abs(res.x - x_ref));
yerr = max(abs(res.y - y_ref));
oerr = abs(res.info.obj_val - obj_ref);
fprintf('test_update_q:  x_err=%.2e  y_err=%.2e  obj_err=%.2e  %s\n', ...
    xerr, yerr, oerr, pass_fail(xerr < tol && yerr < tol && oerr < tol));
if xerr >= tol || yerr >= tol || oerr >= tol, nfail = nfail + 1; end

%% test_update_l
solver2 = osqp(backend='matlab');
solver2.setup(P, q, A, l, u, 'verbose', false, 'polishing', true, ...
    'eps_abs', 1e-8, 'eps_rel', 1e-8, 'linear_solver', 'matlab_ldl');
solver2.update('l', -50 * ones(5, 1));
[x_ref, y_ref, obj_ref] = load_high_accuracy('test_update_l');
res = solver2.solve();
xerr = max(abs(res.x - x_ref));
yerr = max(abs(res.y - y_ref));
oerr = abs(res.info.obj_val - obj_ref);
fprintf('test_update_l:  x_err=%.2e  y_err=%.2e  obj_err=%.2e  %s\n', ...
    xerr, yerr, oerr, pass_fail(xerr < tol && yerr < tol && oerr < tol));
if xerr >= tol || yerr >= tol || oerr >= tol, nfail = nfail + 1; end

%% test_update_max_iter
solver3 = osqp(backend='matlab');
solver3.setup(P, q, A, l, u, 'verbose', false, 'polishing', true, ...
    'eps_abs', 1e-8, 'eps_rel', 1e-8, 'linear_solver', 'matlab_ldl');
solver3.update_settings('max_iter', 80);
res = solver3.solve();
ok = (res.info.status_val == osqp.constant.OSQP_MAX_ITER_REACHED);
fprintf('test_max_iter:  status=%d  %s\n', res.info.status_val, pass_fail(ok));
if ~ok, nfail = nfail + 1; end

%% test_early_termination
solver4 = osqp(backend='matlab');
solver4.setup(P, q, A, l, u, 'verbose', false, 'polishing', true, ...
    'eps_abs', 1e-8, 'eps_rel', 1e-8, 'linear_solver', 'matlab_ldl');
solver4.update_settings('check_termination', 0, 'max_iter', 500);
res = solver4.solve();
ok = (res.info.iter == 500);
fprintf('test_early_term: iter=%d  %s\n', res.info.iter, pass_fail(ok));
if ~ok, nfail = nfail + 1; end

%% test_warm_start
solver5 = osqp(backend='matlab');
solver5.setup(P, q, A, l, u, 'verbose', false, 'eps_abs', 1e-8, ...
    'eps_rel', 1e-8, 'linear_solver', 'matlab_ldl');
solver5.warm_start('x', zeros(2, 1), 'y', zeros(5, 1));
res_cold = solver5.solve();
iter_cold = res_cold.info.iter;
solver5.warm_start('x', res_cold.x, 'y', res_cold.y);
res_warm = solver5.solve();
iter_warm = res_warm.info.iter;
ok = (iter_warm <= iter_cold);
fprintf('test_warm_start: cold=%d warm=%d  %s\n', iter_cold, iter_warm, pass_fail(ok));
if ~ok, nfail = nfail + 1; end

%% Summary
fprintf('\n=== %d failures ===\n', nfail);

function s = pass_fail(ok)
    if ok, s = 'PASS'; else, s = 'FAIL'; end
end
