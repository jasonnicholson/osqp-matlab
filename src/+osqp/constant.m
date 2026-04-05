classdef constant
  % OSQP.CONSTANT  OSQP solver #define constants, mirroring osqp_api_constants.h.
  %
  %   Enumeration types (solver status, polish status, etc.) are defined as separate classdef enumerations:
  %
  %     - osqp.StatusType
  %     - osqp.PolishStatusType
  %     - osqp.CapabilitiesType
  %     - osqp.LinSysSolverType
  %     - osqp.PrecondType
  %     - osqp.ErrorType
  %
  %   Access as class properties:
  %     osqp.constant.OSQP_INFTY
  %     osqp.constant.OSQP_RHO_MIN
  %

  properties
    % Default solver settings
    OSQP_VERBOSE            = 1
    OSQP_WARM_STARTING      = 1
    OSQP_SCALING            = 10
    OSQP_POLISHING          = 0

    % ADMM parameters
    OSQP_RHO                = 0.1
    OSQP_SIGMA              = 1e-6
    OSQP_ALPHA              = 1.6
    OSQP_RHO_MIN            = 1e-6
    OSQP_RHO_MAX            = 1e6
    OSQP_RHO_TOL            = 1e-4    % tolerance for detecting equality constraint
    OSQP_RHO_EQ_OVER_RHO_INEQ = 1e3
    OSQP_RHO_IS_VEC         = 1       % default (non-CUDA build)

    % CG method parameters
    OSQP_CG_MAX_ITER        = 20
    OSQP_CG_TOL_REDUCTION   = 10
    OSQP_CG_TOL_FRACTION    = 0.15
    OSQP_CG_TOL_MIN         = 1e-7
    OSQP_CG_POLISH_TOL      = 1e-5

    % Adaptive rho update methods
    OSQP_ADAPTIVE_RHO_UPDATE_DISABLED    = 0
    OSQP_ADAPTIVE_RHO_UPDATE_ITERATIONS  = 1
    OSQP_ADAPTIVE_RHO_UPDATE_TIME        = 2
    OSQP_ADAPTIVE_RHO_UPDATE_KKT_ERROR   = 3
    OSQP_ADAPTIVE_RHO_UPDATE_DEFAULT     = 1   % = OSQP_ADAPTIVE_RHO_UPDATE_ITERATIONS
    OSQP_ADAPTIVE_RHO_INTERVAL           = 50
    OSQP_ADAPTIVE_RHO_TOLERANCE          = 5.0 % default (non-CUDA build)
    OSQP_ADAPTIVE_RHO_FRACTION           = 0.4
    OSQP_ADAPTIVE_RHO_MULTIPLE_TERMINATION = 4
    OSQP_ADAPTIVE_RHO_FIXED              = 100

    % Termination parameters
    OSQP_MAX_ITER           = 4000
    OSQP_EPS_ABS            = 1e-3
    OSQP_EPS_REL            = 1e-3
    OSQP_EPS_PRIM_INF       = 1e-4
    OSQP_EPS_DUAL_INF       = 1e-4
    OSQP_SCALED_TERMINATION = 0
    OSQP_TIME_LIMIT         = 1e10
    OSQP_CHECK_DUALGAP      = 1        % default (non-float build)
    OSQP_CHECK_TERMINATION  = 25       % default (non-CUDA build)
    OSQP_DELTA              = 1e-6
    OSQP_POLISH_REFINE_ITER = 3

    % Hard-coded values
    OSQP_NAN                = NaN
    OSQP_INFTY              = 1e30     % default (non-CUDA, non-float build)
    OSQP_DIVISION_TOL       = 1e-30    % = 1.0 / OSQP_INFTY
    OSQP_PRINT_INTERVAL     = 200
    OSQP_MIN_SCALING        = 1e-4
    OSQP_MAX_SCALING        = 1e4
    OSQP_ZERO_DEADZONE      = 1e-15    % default (non-float build)
  end

  methods
    function obj = constant()

      arguments
      end
    end
  end
end
