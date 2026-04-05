classdef StatusType
  % OSQP.STATUSTYPE  Solver status codes.
  %   Mirrors osqp_status_type in osqp_api_constants.h.
  %
  %   Example:
  %     if res.info.status_val == int32(osqp.StatusType.OSQP_SOLVED)

  enumeration
    OSQP_SOLVED                       (1)
    OSQP_SOLVED_INACCURATE            (2)
    OSQP_PRIMAL_INFEASIBLE            (3)
    OSQP_PRIMAL_INFEASIBLE_INACCURATE (4)
    OSQP_DUAL_INFEASIBLE              (5)
    OSQP_DUAL_INFEASIBLE_INACCURATE   (6)
    OSQP_MAX_ITER_REACHED             (7)
    OSQP_TIME_LIMIT_REACHED           (8)
    OSQP_NON_CVX                      (9)   % problem non-convex
    OSQP_SIGINT                       (10)  % interrupted by user
    OSQP_UNSOLVED                     (11)  % only setup has been called
  end
end
