% UPDATE_VECTORS Demonstrate vector updates in OSQP v1.0.0
%
% Shows how to update q, l, u vectors without re-setting up the solver.

% Problem data
P = sparse([4 1; 1 2]);
q = [1; 1];
A = sparse([1 1; 1 0; 0 1]);
l = [1; 0; 0];
u = [1; 0.7; 0.7];

% Create solver
solver = osqp;
solver.setup(P, q, A, l, u, 'verbose', false, 'polishing', true);

% Initial solve
res1 = solver.solve();
fprintf('Initial solution: [%.4f, %.4f], obj = %.4f\n', ...
    res1.x(1), res1.x(2), res1.info.obj_val);

% Update linear cost
solver.update('q', [10; 20]);
res2 = solver.solve();
fprintf('After q update:   [%.4f, %.4f], obj = %.4f\n', ...
    res2.x(1), res2.x(2), res2.info.obj_val);

% Update bounds
solver.update('l', [0.5; 0; 0], 'u', [1.5; 0.9; 0.9]);
res3 = solver.solve();
fprintf('After bound update: [%.4f, %.4f], obj = %.4f\n', ...
    res3.x(1), res3.x(2), res3.info.obj_val);

% Clean up
delete(solver);
