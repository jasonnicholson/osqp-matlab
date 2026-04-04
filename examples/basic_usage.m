% BASIC_USAGE Basic usage of OSQP v1.0.0
%
% Demonstrates setting up, solving, and updating a QP problem.

% Problem data
P = sparse([4 1; 1 2]);
q = [1; 1];
A = sparse([1 1; 1 0; 0 1]);
l = [1; 0; 0];
u = [1; 0.7; 0.7];

% Create solver and setup problem
solver = osqp;
solver.setup(P, q, A, l, u, ...
    'verbose', false, ...
    'polishing', true, ...
    'eps_abs', 1e-8, ...
    'eps_rel', 1e-8);

% Solve
results = solver.solve();

% Check solution
fprintf('Status: %s\n', results.info.status);
fprintf('Objective value: %.6f\n', results.info.obj_val);
fprintf('Primal solution: [%.6f, %.6f]\n', results.x(1), results.x(2));
fprintf('Dual solution:   [%.6f, %.6f, %.6f]\n', results.y(1), results.y(2), results.y(3));
fprintf('Iterations: %d\n', results.info.iter);
fprintf('Solve time: %.4f ms\n', results.info.solve_time * 1000);

% Update settings and re-solve
solver.update_settings('max_iter', 200);
results2 = solver.solve();
fprintf('Re-solve iterations: %d\n', results2.info.iter);

% Print version and capabilities
fprintf('\nOSQP version: %s\n', osqp.version());
fprintf('Capabilities: %d\n', osqp.capabilities());

% Clean up
delete(solver);
