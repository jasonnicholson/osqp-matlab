classdef Options
    % OSQP.OPTIONS  Validated settings for the OSQP solver.
    %
    %   opts = osqp.Options          % all defaults
    %   opts.max_iter = 1000;        % modify after construction
    %   s = opts.toStruct();         % for MEX passthrough
    %
    %   Default values match OSQP v1.0.0.

    % -----------------------------------------------------------------
    %  ADMM parameters
    % -----------------------------------------------------------------
    properties
        rho     (1,1) double {mustBePositive}              = 0.1
        sigma   (1,1) double {mustBePositive}              = 1e-6
        alpha   (1,1) double {mustBeInRange(alpha,0,2)}    = 1.6
    end

    % -----------------------------------------------------------------
    %  Termination criteria
    % -----------------------------------------------------------------
    properties
        max_iter          (1,1) double {mustBePositive, mustBeInteger} = 4000
        eps_abs           (1,1) double {mustBeNonnegative}             = 1e-3
        eps_rel           (1,1) double {mustBeNonnegative}             = 1e-3
        eps_prim_inf      (1,1) double {mustBeNonnegative}             = 1e-4
        eps_dual_inf      (1,1) double {mustBeNonnegative}             = 1e-4
        check_termination (1,1) double {mustBeNonnegative, mustBeInteger} = 25
        time_limit        (1,1) double {mustBeNonnegative}             = 0
        scaled_termination(1,1) logical                                = false
        check_dualgap     (1,1) logical                                = false
    end

    % -----------------------------------------------------------------
    %  Scaling
    % -----------------------------------------------------------------
    properties
        scaling (1,1) double {mustBeNonnegative, mustBeInteger} = 10
    end

    % -----------------------------------------------------------------
    %  Adaptive rho
    % -----------------------------------------------------------------
    properties
        adaptive_rho           (1,1) logical = true
        adaptive_rho_interval  (1,1) double {mustBeNonnegative, mustBeInteger} = 0
        adaptive_rho_tolerance (1,1) double {mustBeGreaterThanOrEqual(adaptive_rho_tolerance,1)} = 5
        adaptive_rho_fraction  (1,1) double {mustBePositive}  = 0.4
    end

    % -----------------------------------------------------------------
    %  Polishing
    % -----------------------------------------------------------------
    properties
        polishing          (1,1) logical                                 = false
        delta              (1,1) double {mustBePositive}                 = 1e-6
        polish_refine_iter (1,1) double {mustBeNonnegative, mustBeInteger} = 3
    end

    % -----------------------------------------------------------------
    %  Linear solver
    % -----------------------------------------------------------------
    properties
        linear_solver (1,:) char {mustBeMember(linear_solver, ...
            {'matlab_ldl','qdldl','qdldl_c'})} = 'qdldl_c'
    end

    % -----------------------------------------------------------------
    %  Output
    % -----------------------------------------------------------------
    properties
        verbose (1,1) logical = true
    end

    % -----------------------------------------------------------------
    %  Warm starting
    % -----------------------------------------------------------------
    properties
        warm_starting (1,1) logical = true
    end

    % -----------------------------------------------------------------
    %  C-interface-only settings (passed through, not used by Solver.m)
    % -----------------------------------------------------------------
    properties
        device            (1,1) double {mustBeNonnegative, mustBeInteger} = 0
        linsys_solver     (1,1) double {mustBeNonnegative, mustBeInteger} = 0
        allocate_solution (1,1) logical                                  = true
        profiler_level    (1,1) double {mustBeNonnegative, mustBeInteger} = 0
        rho_is_vec        (1,1) logical                                  = false
        cg_max_iter       (1,1) double {mustBePositive, mustBeInteger}    = 20
        cg_tol_reduction  (1,1) double {mustBePositive, mustBeInteger}    = 200
        cg_tol_fraction   (1,1) double {mustBePositive}                  = 2.5e-5
        cg_precond        (1,1) double {mustBeNonnegative, mustBeInteger} = 0
    end

    % -----------------------------------------------------------------
    %  Constants
    % -----------------------------------------------------------------
    properties (Constant)
        OSQP_INFTY = 1e30

        % Settings that can only be set at setup time (not updatable)
        SETUP_ONLY_SETTINGS = { ...
            'sigma', 'scaling', ...
            'adaptive_rho', 'adaptive_rho_interval', ...
            'adaptive_rho_fraction', 'adaptive_rho_tolerance', ...
            'linear_solver', ...
            'device', 'linsys_solver', 'allocate_solution', ...
            'profiler_level', 'rho_is_vec', ...
            'cg_max_iter', 'cg_tol_reduction', 'cg_tol_fraction', 'cg_precond'}

        % Settings that can be updated after setup
        UPDATABLE_SETTINGS = { ...
            'max_iter', 'eps_abs', 'eps_rel', ...
            'eps_prim_inf', 'eps_dual_inf', ...
            'time_limit', 'rho', 'alpha', ...
            'delta', 'polishing', 'polish_refine_iter', ...
            'verbose', 'check_termination', 'check_dualgap', ...
            'warm_starting', 'scaled_termination'}
    end

    methods
        function s = toStruct(obj)
            % TOSTRUCT  Convert to plain struct for MEX passthrough.
            mc = metaclass(obj);
            props = mc.PropertyList;
            s = struct();
            for k = 1:numel(props)
                p = props(k);
                if p.Constant || p.Dependent
                    continue
                end
                s.(p.Name) = obj.(p.Name);
            end
        end

        function tf = isUpdatable(~, name)
            % ISUPDATABLE  Return true if the named setting can be changed after setup.
            arguments
                ~
                name (1,:) char
            end
            tf = ismember(name, osqp.Options.UPDATABLE_SETTINGS);
        end

        function tf = isSetupOnly(~, name)
            % ISSETUPONLY  Return true if the named setting can only be set at setup.
            arguments
                ~
                name (1,:) char
            end
            tf = ismember(name, osqp.Options.SETUP_ONLY_SETTINGS);
        end
    end

    methods (Static)
        function obj = fromStruct(s)
            % FROMSTRUCT  Create Options from a struct (e.g. returned by MEX).
            arguments
                s (1,1) struct
            end
            obj = osqp.Options();
            fnames = fieldnames(s);
            mc = metaclass(obj);
            propNames = {mc.PropertyList.Name};
            for k = 1:numel(fnames)
                fn = fnames{k};
                if ismember(fn, propNames)
                    try
                        obj.(fn) = s.(fn);
                    catch
                        % Skip fields that fail validation (e.g. unknown enums)
                    end
                end
            end
        end
    end
end
