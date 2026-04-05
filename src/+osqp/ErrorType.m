classdef ErrorType
    % OSQP.ERRORTYPE  Solver error codes.
    %   Mirrors osqp_error_type in osqp_api_constants.h.

    enumeration
        OSQP_NO_ERROR                  (0)
        OSQP_DATA_VALIDATION_ERROR     (1)
        OSQP_SETTINGS_VALIDATION_ERROR (2)
        OSQP_LINSYS_SOLVER_INIT_ERROR  (3)
        OSQP_NONCVX_ERROR              (4)
        OSQP_MEM_ALLOC_ERROR           (5)
        OSQP_WORKSPACE_NOT_INIT_ERROR  (6)
        OSQP_ALGEBRA_LOAD_ERROR        (7)
        OSQP_FOPEN_ERROR               (8)
        OSQP_CODEGEN_DEFINES_ERROR     (9)
        OSQP_DATA_NOT_INITIALIZED      (10)
        OSQP_FUNC_NOT_IMPLEMENTED      (11)  % function not implemented in this library
        OSQP_LAST_ERROR_PLACE          (12)  % sentinel — must remain last
    end
end
