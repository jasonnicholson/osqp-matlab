% CODE_GENERATION Demonstrate native code generation in OSQP v1.0.0
%
% Generates embedded C code using the native osqp_codegen() API.

% Problem data
P = sparse([4 1; 1 2]);
q = [1; 1];
A = sparse([1 1; 1 0; 0 1]);
l = [1; 0; 0];
u = [1; 0.7; 0.7];

% Create and solve
solver = osqp;
solver.setup(P, q, A, l, u, 'verbose', false);
res = solver.solve();
fprintf('Solution: [%.4f, %.4f]\n', res.x(1), res.x(2));

% Check codegen capability
if ~solver.has_capability('OSQP_CAPABILITY_CODEGEN')
    fprintf('Code generation is not available in this build.\n');
    delete(solver);
    return;
end

% Generate code — vectors only (embedded mode 1)
solver.codegen('codegen_vectors', ...
    'prefix', '', ...
    'parameters', 'vectors', ...
    'force_rewrite', true, ...
    'float_type', false, ...
    'printing_enable', false, ...
    'profiling_enable', false, ...
    'interrupt_enable', false);
fprintf('Generated code (vectors only) in codegen_vectors/\n');

% Generate code — vectors + matrices (embedded mode 2)
solver.codegen('codegen_matrices', ...
    'prefix', 'qp_', ...
    'parameters', 'matrices', ...
    'force_rewrite', true);
fprintf('Generated code (vectors + matrices) in codegen_matrices/\n');

% Clean up
delete(solver);
