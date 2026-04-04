% ADJOINT_DERIVATIVES Demonstrate adjoint derivatives in OSQP v1.0.0
%
% Shows how to compute gradients of a loss function with respect to
% problem data using adjoint differentiation.

% Check if derivatives are supported
solver = osqp;
if ~solver.has_capability('OSQP_CAPABILITY_DERIVATIVES')
    fprintf('Derivatives are not available in this build.\n');
    delete(solver);
    return;
end

% Problem data: min 0.5 x'Px + q'x s.t. l <= Ax <= u
P = sparse([4 1; 1 2]);
q = [1; 1];
A = sparse([1 1; 1 0; 0 1]);
l = [1; 0; 0];
u = [1; 0.7; 0.7];

% Set up and solve with tight tolerances
solver = osqp;
solver.setup(P, q, A, l, u, ...
    'verbose', false, ...
    'polishing', true, ...
    'eps_abs', 1e-9, ...
    'eps_rel', 1e-9);
res = solver.solve();

fprintf('Solution: x = [%.6f, %.6f]\n', res.x(1), res.x(2));
fprintf('Dual:     y = [%.6f, %.6f, %.6f]\n', res.y(1), res.y(2), res.y(3));

%% Compute adjoint derivatives
% Loss function: L(x, y) = 0.5 * ||x - x_target||^2
x_target = [0.5; 0.5];
dx = res.x - x_target;  % dL/dx
dy = zeros(3, 1);       % dL/dy = 0

solver.adjoint_derivative_compute(dx, dy);

% Get derivative matrices (dL/dP, dL/dA)
[dP, dA] = solver.adjoint_derivative_get_mat();
fprintf('\ndL/dP:\n');
disp(full(dP));
fprintf('dL/dA:\n');
disp(full(dA));

% Get derivative vectors (dL/dq, dL/dl, dL/du)
[dq, dl, du] = solver.adjoint_derivative_get_vec();
fprintf('dL/dq: [%.6f, %.6f]\n', dq(1), dq(2));
fprintf('dL/dl: [%.6f, %.6f, %.6f]\n', dl(1), dl(2), dl(3));
fprintf('dL/du: [%.6f, %.6f, %.6f]\n', du(1), du(2), du(3));

% Clean up
delete(solver);
