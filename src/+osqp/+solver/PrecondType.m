classdef PrecondType
  % OSQP.PRECONDTYPE  Preconditioner types for the CG method.
  %   Mirrors osqp_precond_type in osqp_api_constants.h.

  enumeration
    OSQP_NO_PRECONDITIONER       (0)
    OSQP_DIAGONAL_PRECONDITIONER (1)  % diagonal (Jacobi) preconditioner
  end
end
