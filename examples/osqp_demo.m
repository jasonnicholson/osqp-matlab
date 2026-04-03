% Demo showing the usage of OSQP from Matlab and the code generation features.
% This problem is the same one that is presented in the osqp_demo.c file.
% Updated for OSQP v1.0.0.

% The problem data
P = sparse([4., 1.; 1., 2.]);
q = [1; 1];
A = sparse([1., 1; 1, 0; 0, 1]);
l = [1.0; 0.0; 0.0];
u = [1.0; 0.7; 0.7];

% Create the solver
solver = osqp;
solver.setup(P, q, A, l, u, 'verbose', true)

% Solve the problem
results = solver.solve()

% Print solution
fprintf('Status: %s\n', results.info.status);
fprintf('Objective: %.4f\n', results.info.obj_val);
fprintf('x = [%.4f, %.4f]\n', results.x(1), results.x(2));

% Generate embedded C code (if codegen capability is available)
if solver.has_capability('OSQP_CAPABILITY_CODEGEN')
    solver.codegen('osqp_demo', 'parameters', 'vectors', 'force_rewrite', true);
    fprintf('Code generated in osqp_demo/\n');
end

% Clean up
delete(solver);
