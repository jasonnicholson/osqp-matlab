function [m, n, P, q, A, l, u] = validateData(P, q, A, l, u, infty)
% validateData Infer and validate problem data for OSQP solver ::
%
%   [m, n, P, q, A, l, u] = validateData(P, q, A, l, u, infty)
%
% Inputs:
%   P: Quadratic cost matrix (sparse, n x n). If empty, set to zero.
%   q: Linear cost vector (n x 1). If empty, set to zeros.
%   A: Linear constraint matrix (sparse, m x n). If empty, no constraints.
%   l: Lower bound vector (m x 1). If empty and A specified, set to -inf.
%   u: Upper bound vector (m x 1). If empty and A specified, set to inf.
%   infty: Infinity value for clamping bounds.
%
% Output arguments:
%   m: Number of constraints (rows in A).
%   n: Number of variables (columns in A or length of q).
%   P: Validated quadratic cost matrix (sparse, upper triangular).
%   q: Validated linear cost vector.
%   A: Validated constraint matrix (sparse).
%   l: Validated lower bound vector, clamped to -infty.
%   u: Validated upper bound vector, clamped to infty.
%
% Notes:
%   - This function mirrors _infer_mnpqalu from the Python interface.
%   - P and A must be sparse matrices if provided as numeric matrices.
%   - P is made upper triangular if it has lower triangular elements.
%   - Bounds L and U are clamped.


  if isempty(P)
    if ~isempty(q)
      n = length(q);
    elseif ~isempty(A)
      n = size(A, 2);
    else
      error('The problem does not have any variables');
    end
  else
    n = size(P, 1);
  end

  m = 0;
  if ~isempty(A)
    m = size(A, 1);
  end

  if isempty(A)
    assert(isempty(l) && isempty(u), 'If A is unspecified, leave l/u unspecified too.');
  else
    assert(~isempty(l) || ~isempty(u), 'If A is specified, specify at least one of l/u.');
    if isempty(l)
      l = -inf * ones(size(A, 1), 1);
    end
    if isempty(u)
      u = inf * ones(size(A, 1), 1);
    end
  end

  if isempty(P)
    P = sparse(n, n);
  end

  if isempty(q)
    q = zeros(n, 1);
  end

  if isempty(A)
    A = sparse(m, n);
    l = zeros(m, 1);
    u = zeros(m, 1);
  end

  assert(length(q) == n, 'Incorrect dimension of q');
  assert(length(l) == m, 'Incorrect dimension of l');
  assert(length(u) == m, 'Incorrect dimension of u');

  % Type checks: P and A should be sparse matrices
  if ~issparse(P) && isnumeric(P) && ismatrix(P)
    error('P is required to be a sparse matrix');
  end
  if ~issparse(A) && isnumeric(A) && ismatrix(A)
    error('A is required to be a sparse matrix');
  end

  % If P has lower triangular elements, take upper triangular part
  if ~istriu(P)
    P = triu(P);
  end

  % Assume P and A are already in CSC form (sparse in MATLAB is CSC-like)
  % No need to convert, as MATLAB sparse is CSC

  % Clamp l and u to bounds
  u = min(u, infty);
  l = max(l, -infty);

end
