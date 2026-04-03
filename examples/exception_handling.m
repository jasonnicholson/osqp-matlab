% EXCEPTION_HANDLING Error handling in OSQP v1.0.0
%
% Demonstrates how to handle errors from OSQP using try/catch.

%% Non-convex problem (setup error)
P_non_cvx = sparse([2 5; 5 1]);  % Not positive semidefinite
q = [3; 4];
A = sparse([-1 0; 0 -1; -1 -3; 2 5; 3 4]);
l = -1e30 * ones(5, 1);
u = [0; 0; -15; 100; 80];

solver = osqp;
try
    solver.setup(P_non_cvx, q, A, l, u, 'verbose', false, 'sigma', 1e-6);
    fprintf('Setup succeeded (unexpected)\n');
catch ME
    fprintf('Setup failed as expected: %s\n', ME.message);
end
delete(solver);

%% Infeasible problem (solve returns status)
P = sparse([4 1; 1 2]);
q = [1; 1];
A = sparse([1 0; 0 1; 1 0]);
l = [1; 1; -0.5];  % Conflicting: x1 >= 1 and x1 <= -0.5
u = [-0.5; 2; 1];

solver = osqp;
solver.setup(P, q, A, l, u, 'verbose', false);
res = solver.solve();

if res.info.status_val == osqp.constant('OSQP_SOLVED')
    fprintf('Problem solved: obj = %.4f\n', res.info.obj_val);
elseif res.info.status_val == osqp.constant('OSQP_PRIMAL_INFEASIBLE')
    fprintf('Problem is primal infeasible\n');
elseif res.info.status_val == osqp.constant('OSQP_DUAL_INFEASIBLE')
    fprintf('Problem is dual infeasible\n');
else
    fprintf('Solver status: %s (%d)\n', res.info.status, res.info.status_val);
end

delete(solver);
