classdef constant
    % OSQP.CONSTANT  OSQP solver constants.
    %
    %   Access as class properties:
    %     osqp.constant.OSQP_INFTY
    %     osqp.constant.OSQP_SOLVED
    %
    %   Legacy function-call syntax (backward compatible):
    %     val = osqp.constant('OSQP_SOLVED')

    properties (Constant)
        OSQP_INFTY                        = 1e30

        % Solver status (matches osqp_api_constants.h enum, v1.0.0)
        OSQP_SOLVED                       = 1
        OSQP_SOLVED_INACCURATE            = 2
        OSQP_PRIMAL_INFEASIBLE            = 3
        OSQP_PRIMAL_INFEASIBLE_INACCURATE = 4
        OSQP_DUAL_INFEASIBLE              = 5
        OSQP_DUAL_INFEASIBLE_INACCURATE   = 6
        OSQP_MAX_ITER_REACHED             = 7
        OSQP_TIME_LIMIT_REACHED           = 8
        OSQP_NON_CONVEX                   = 9
        OSQP_SIGINT                       = 10
        OSQP_UNSOLVED                     = 11

        % Polish status
        OSQP_POLISH_LINSYS_ERROR          = -2
        OSQP_POLISH_FAILED                = -1
        OSQP_POLISH_NOT_PERFORMED         = 0
        OSQP_POLISH_SUCCESS               = 1
        OSQP_POLISH_NO_ACTIVE_SET_FOUND   = 2

        % Algorithm constants (matching C code constants.h)
        OSQP_INFTY_THRESH                 = 1e30 * 1e-4   % OSQP_INFTY * OSQP_MIN_SCALING
        OSQP_RHO_MIN                      = 1e-6
        OSQP_RHO_MAX                      = 1e6
        OSQP_RHO_EQ_OVER_RHO_INEQ        = 1e3
        OSQP_RHO_TOL                      = 1e-4
        OSQP_MIN_SCALING                  = 1e-4
        OSQP_MAX_SCALING                  = 1e4
        OSQP_DIVISION_TOL                 = 1e-10
        OSQP_ADAPTIVE_RHO_MULTIPLE_TERMINATION = 4
    end

    methods (Static)
        function out = lookup(name)
            % LOOKUP  Query a constant by name string.
            %
            %   val = osqp.constant.lookup('OSQP_SOLVED')
            try
                out = osqp.constant.(upper(name));
            catch
                if osqp.CInterface.isMexAvailable()
                    out = osqp.CInterface.constant(name);
                else
                    error('OSQP:constant', 'Unknown constant: %s', name);
                end
            end
        end
    end

    methods
        function obj = constant(name)
            % For backward compatibility: osqp.constant('OSQP_SOLVED')
            % returns the value via subsref overriding if called with
            % an argument.
            if nargin > 0
                error('OSQP:constant', ...
                    ['Use osqp.constant.%s or osqp.constant.lookup(''%s'') ' ...
                     'instead of osqp.constant(''%s'').'], name, name, name);
            end
        end
    end
end
