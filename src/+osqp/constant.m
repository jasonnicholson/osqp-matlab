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
        OSQP_INFTY                       = 1e30
        OSQP_SOLVED                      = 1
        OSQP_SOLVED_INACCURATE           = 2
        OSQP_PRIMAL_INFEASIBLE           = -3
        OSQP_PRIMAL_INFEASIBLE_INACCURATE = 3
        OSQP_DUAL_INFEASIBLE             = -4
        OSQP_DUAL_INFEASIBLE_INACCURATE  = 4
        OSQP_MAX_ITER_REACHED            = -2
        OSQP_TIME_LIMIT_REACHED          = -6
        OSQP_NON_CONVEX                  = -7
        OSQP_UNSOLVED                    = -10
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
