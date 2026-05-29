classdef SolverOptions
    % OSQP.SOLVEROPTIONS  Validated settings for the OSQP solver.
    %
    %   opts = osqp.SolverOptions          % all defaults
    %   opts.max_iter = 1000;              % modify after construction
    %   s = opts.toStruct();               % for MEX passthrough
    %
    %   Default values match OSQP.

    % -----------------------------------------------------------------
    %  ADMM parameters
    % -----------------------------------------------------------------
    properties
        % ADMM penalty parameter rho.
        rho     (1,1) double {mustBePositive} = 0.1
        % Dual variable regularization parameter sigma.
        sigma   (1,1) double {mustBePositive} = 1e-6
        % Relaxation parameter alpha in the interval (0, 2).
        alpha   (1,1) double {mustBeBetween(alpha, 0, 2)} = 1.6
    end

    % -----------------------------------------------------------------
    %  Termination criteria
    % -----------------------------------------------------------------
    properties
        % Maximum number of ADMM iterations.
        max_iter (1,1) double {mustBePositive, mustBeInteger} = 4000
        % Absolute feasibility tolerance.
        eps_abs  (1,1) double {mustBeNonnegative} = 1e-3
        % Relative feasibility tolerance.
        eps_rel (1,1) double {mustBeNonnegative}  = 1e-3
        % Primal infeasibility tolerance.
        eps_prim_inf (1,1) double {mustBeNonnegative} = 1e-4
        % Dual infeasibility tolerance.
        eps_dual_inf (1,1) double {mustBeNonnegative} = 1e-4
        % Iteration interval for termination checks.
        check_termination (1,1) double {mustBeNonnegative, mustBeInteger} = 25
        % Solver time limit in seconds; 0 disables the limit.
        time_limit (1,1) double {mustBeNonnegative} = 0
        % Use scaled stopping criteria when true.
        scaled_termination (1,1) logical = false
        % Include duality gap in termination checks when true.
        check_dualgap (1,1) logical = false
    end

    % -----------------------------------------------------------------
    %  Scaling
    % -----------------------------------------------------------------
    properties
        % Number of Ruiz scaling iterations.
        scaling (1,1) double {mustBeNonnegative, mustBeInteger} = 10
    end

    % -----------------------------------------------------------------
    %  Adaptive rho
    % -----------------------------------------------------------------
    properties
        % Enable adaptive rho updates when true.
        adaptive_rho (1,1) logical = true
        % Iteration interval for adaptive rho updates; 0 enables auto mode.
        adaptive_rho_interval  (1,1) double {mustBeNonnegative, mustBeInteger} = 0
        % Minimum ratio threshold used to trigger adaptive rho updates.
        adaptive_rho_tolerance (1,1) double {mustBeGreaterThanOrEqual(adaptive_rho_tolerance,1)} = 5
        % Fraction of setup time used to choose automatic rho update interval.
        adaptive_rho_fraction  (1,1) double {mustBePositive}  = 0.4
    end

    % -----------------------------------------------------------------
    %  Polishing
    % -----------------------------------------------------------------
    properties
        % Enable solution polishing when true.
        polishing (1,1) logical = false
        % Regularization parameter used by polishing.
        delta (1,1) double {mustBePositive} = 1e-6
        % Number of iterative refinement steps during polishing.
        polish_refine_iter (1,1) double {mustBeNonnegative, mustBeInteger} = 3
    end

    % -----------------------------------------------------------------
    %  Linear solver
    % -----------------------------------------------------------------
    properties
        % TODO Better handle how the Cinterface and native solver linear algebra backends are different.
        % Linear system backend used by the solver.
        linear_solver (1,:) char {mustBeMember(linear_solver, {'matlab_ldl', 'qdldl', 'qdldl_c', 'mkl pardiso', 'cuda pcg'})} = 'qdldl'
    end

    % -----------------------------------------------------------------
    %  Output
    % -----------------------------------------------------------------
    properties
        % Print solver progress to the command window when true.
        verbose (1,1) logical = true
    end

    % -----------------------------------------------------------------
    %  Warm starting
    % -----------------------------------------------------------------
    properties
        % Reuse previous primal/dual iterates to warm start when true.
        warm_starting (1,1) logical = true
    end

    % -----------------------------------------------------------------
    %  C-interface-only settings (passed through, not used by Solver.m)
    % -----------------------------------------------------------------
    properties
        % Device identifier used by selected backend implementations.
        device (1,1) double {mustBeNonnegative, mustBeInteger} = 0
        % Numeric code for the low-level linear system solver.
        linsys_solver (1,1) double {mustBeNonnegative, mustBeInteger} = 0
        % Allocate internal solution buffers at setup when true.
        allocate_solution (1,1) logical = true
        % Profiler detail level for C-interface instrumentation.
        profiler_level (1,1) double {mustBeNonnegative, mustBeInteger} = 0
        % Interpret rho as a vector when true.
        rho_is_vec (1,1) logical = false
        % Maximum conjugate-gradient iterations per linear solve.
        cg_max_iter (1,1) double {mustBePositive, mustBeInteger}    = 20
        % Integer factor controlling CG tolerance reduction.
        cg_tol_reduction  (1,1) double {mustBePositive, mustBeInteger} = 200
        % Fraction used to set the CG convergence tolerance.
        cg_tol_fraction (1,1) double {mustBePositive} = 2.5e-5
        % Preconditioner selection code for the CG solver.
        cg_precond (1,1) double {mustBeNonnegative, mustBeInteger} = 0
    end

    % -----------------------------------------------------------------
    %  Constants
    % -----------------------------------------------------------------
    properties (Constant)
        % Settings that can only be set at setup time (not updatable)
        % Names of settings that are only valid during initial setup.
        SETUP_ONLY_SETTINGS = { ...
            'sigma', 'scaling', ...
            'adaptive_rho', 'adaptive_rho_interval', ...
            'adaptive_rho_fraction', 'adaptive_rho_tolerance', ...
            'linear_solver', ...
            'device', 'linsys_solver', 'allocate_solution', ...
            'profiler_level', 'rho_is_vec', ...
            'cg_max_iter', 'cg_tol_reduction', 'cg_tol_fraction', 'cg_precond'}

        % Settings that can be updated after setup
        % Names of settings that can be updated after initial setup.
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
            tf = ismember(name, osqp.SolverOptions.UPDATABLE_SETTINGS);
        end

        function tf = isSetupOnly(~, name)
            % ISSETUPONLY  Return true if the named setting can only be set at setup.
            arguments
                ~
                name (1,:) char
            end
            tf = ismember(name, osqp.SolverOptions.SETUP_ONLY_SETTINGS);
        end
    end

    methods (Static)
        function obj = fromStruct(s)
            % FROMSTRUCT  Create Options from a struct (e.g. returned by MEX).
            arguments
                s (1,1) struct
            end
            obj = osqp.SolverOptions();
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
