classdef PolishStatusType < double
  % OSQP.POLISHSTATUSTYPE  Polish status codes.
  %   Mirrors osqp_polish_status_type in osqp_api_constants.h.
  %

  enumeration
    OSQP_POLISH_LINSYS_ERROR        (-2)
    OSQP_POLISH_FAILED              (-1)
    OSQP_POLISH_NOT_PERFORMED       (0)
    OSQP_POLISH_SUCCESS             (1)
    OSQP_POLISH_NO_ACTIVE_SET_FOUND (2)  % no active set detected, polishing skipped
  end
end
