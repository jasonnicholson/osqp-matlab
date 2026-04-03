% UPDATE_MATRICES Demonstrate matrix updates in OSQP v1.0.0
%
% Shows how to update P and A matrix values (same sparsity pattern).

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

% Update P values (upper triangular, same sparsity pattern)
% P_new = [6 2; 2 3] → triu = [6 2; 0 3] → nonzeros = [6, 2, 3]
Px_new = [6; 2; 3];
solver.update('Px', Px_new);
res2 = solver.solve();
fprintf('After P update:   [%.4f, %.4f], obj = %.4f\n', ...
    res2.x(1), res2.x(2), res2.info.obj_val);

% Update A values (same sparsity pattern)
% A_new = [2 1; 1 0; 0 2] → nonzeros = [2, 1, 1, 2]
Ax_new = [2; 1; 1; 2];
solver.update('Ax', Ax_new);
res3 = solver.solve();
fprintf('After A update:   [%.4f, %.4f], obj = %.4f\n', ...
    res3.x(1), res3.x(2), res3.info.obj_val);

% Clean up
delete(solver);
