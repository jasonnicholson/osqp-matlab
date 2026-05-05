classdef LinSysSolverType < double
  % OSQP.LINSSSOLVERTYPE  Linear system solver types.
  %   Mirrors osqp_linsys_solver_type in osqp_api_constants.h.

  enumeration
    OSQP_UNKNOWN_SOLVER  (0)
    OSQP_DIRECT_SOLVER   (1)
    OSQP_INDIRECT_SOLVER (2)
  end
end
