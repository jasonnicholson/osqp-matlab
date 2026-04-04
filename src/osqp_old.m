classdef osqp < handle
    % OSQP interface class for OSQP solver v1.0.0
    %
    % This class provides a complete interface to the C implementation
    % of the OSQP solver, with an optional pure-MATLAB backend using
    % MATLAB's built-in ldl() for the KKT linear system solves.
    %
    % To use the pure-MATLAB backend, pass 'linear_solver','matlab_ldl'
    % as a setup option:
    %   solver = osqp;
    %   solver.setup(P, q, A, l, u, 'linear_solver', 'matlab_ldl');
    %
    % osqp Properties:
    %   objectHandle - pointer to the C structure of OSQP solver
    %
    % osqp Methods:
    %
    %   setup             - configure solver with problem data
    %   solve             - solve the QP
    %   update            - modify problem vectors and/or matrices
    %   warm_start        - set warm starting variables x and y
    %   cold_start        - reset solver iterates to zero
    %
    %   default_settings  - create default settings structure
    %   current_settings  - get the current solver settings structure
    %   update_settings   - update the current solver settings structure
    %
    %   get_dimensions    - get the number of variables and constraints
    %   version           - return OSQP version
    %   constant          - return an OSQP internal constant
    %   capabilities      - return solver capability bitmask
    %   has_capability    - check if a specific capability is available
    %
    %   codegen           - generate embeddable C code for the problem
    %
    %   adjoint_derivative_compute - compute adjoint derivatives
    %   adjoint_derivative_get_mat - retrieve derivative matrices dP, dA
    %   adjoint_derivative_get_vec - retrieve derivative vectors dq, dl, du

    properties (SetAccess = private, Hidden = true)
        objectHandle % Handle to underlying C instance
    end

    % Properties for the pure-MATLAB backend
    properties (Access = private, Hidden = true)
        use_matlab_backend = false

        % Problem data (unscaled)
        ml_P_triu
        ml_q
        ml_A
        ml_l
        ml_u
        ml_n = 0
        ml_m = 0

        % Settings
        ml_settings

        % Scaling
        ml_scl

        % Scaled problem
        ml_Ps
        ml_qs
        ml_As
        ml_ls
        ml_us

        % ADMM iterates
        ml_x
        ml_z
        ml_y

        % Rho vectors
        ml_rho_vec
        ml_rho_inv_vec
        ml_constr_type

        % KKT factorization (MATLABLDLSolver)
        ml_kkt_factor

        % Flags / timing
        ml_isSetup = false
        ml_non_convex = false
        ml_setup_time = 0
        ml_solve_time = 0
        ml_update_time = 0
        ml_polish_time = 0
    end

    % Constants for status values
    properties (Constant, Hidden = true)
        OSQP_INFTY = 1e30
        ML_STATUS_DUAL_INFEASIBLE_INACCURATE  = 4
        ML_STATUS_PRIMAL_INFEASIBLE_INACCURATE = 3
        ML_STATUS_SOLVED_INACCURATE            = 2
        ML_STATUS_SOLVED                       = 1
        ML_STATUS_MAX_ITER_REACHED             = -2
        ML_STATUS_PRIMAL_INFEASIBLE            = -3
        ML_STATUS_DUAL_INFEASIBLE              = -4
        ML_STATUS_TIME_LIMIT_REACHED           = -6
        ML_STATUS_NON_CONVEX                   = -7
        ML_STATUS_UNSOLVED                     = -10
    end

    methods(Static)

        function out = default_settings()
            % DEFAULT_SETTINGS get the default solver settings structure
            try
                out = osqp_mex('default_settings', 'static');
                out.linsys_solver = linsys_solver_to_string(out.linsys_solver);
            catch
                out = osqp.matlab_ldl_default_settings();
            end
        end

        function out = constant(constant_name)
            % CONSTANT Return solver constant
            %   C = CONSTANT(CONSTANT_NAME) return constant called CONSTANT_NAME
            switch upper(constant_name)
                case 'OSQP_INFTY'
                    out = 1e30;
                case 'OSQP_SOLVED'
                    out = 1;
                case 'OSQP_SOLVED_INACCURATE'
                    out = 2;
                case 'OSQP_PRIMAL_INFEASIBLE'
                    out = -3;
                case 'OSQP_PRIMAL_INFEASIBLE_INACCURATE'
                    out = 3;
                case 'OSQP_DUAL_INFEASIBLE'
                    out = -4;
                case 'OSQP_DUAL_INFEASIBLE_INACCURATE'
                    out = 4;
                case 'OSQP_MAX_ITER_REACHED'
                    out = -2;
                case 'OSQP_TIME_LIMIT_REACHED'
                    out = -6;
                case 'OSQP_NON_CONVEX'
                    out = -7;
                case 'OSQP_UNSOLVED'
                    out = -10;
                otherwise
                    try
                        out = osqp_mex('constant', 'static', constant_name);
                    catch
                        error('OSQP:constant', 'Unknown constant: %s', constant_name);
                    end
            end
        end

        function out = version()
            % VERSION Return OSQP version string
            try
                out = osqp_mex('version', 'static');
            catch
                out = '1.0.0-matlab';
            end
        end

        function out = capabilities()
            % CAPABILITIES Return capability bitmask
            try
                out = osqp_mex('capabilities', 'static');
            catch
                out = 0;
            end
        end

    end

    methods

        %% Constructor
        function this = osqp(varargin)
            % OSQP Construct OSQP solver class
            %   solver = osqp()
            %   solver = osqp('linear_solver', 'matlab_ldl')
            this.use_matlab_backend = false;
            if nargin >= 2
                for k = 1:2:nargin
                    if ischar(varargin{k}) && strcmpi(varargin{k}, 'linear_solver') ...
                            && ischar(varargin{k+1}) && strcmpi(varargin{k+1}, 'matlab_ldl')
                        this.use_matlab_backend = true;
                    end
                end
            end
            if ~this.use_matlab_backend
                this.objectHandle = osqp_mex('new', varargin{:});
            end
        end

        %% Destructor
        function delete(this)
            % DELETE Destroy OSQP solver class
            if ~this.use_matlab_backend && ~isempty(this.objectHandle)
                osqp_mex('delete', this.objectHandle);
            end
        end

        %% current_settings
        function out = current_settings(this)
            % CURRENT_SETTINGS get the current solver settings structure
            if this.use_matlab_backend
                out = this.ml_settings;
                return;
            end
            out = osqp_mex('current_settings', this.objectHandle);
            out.linsys_solver = linsys_solver_to_string(out.linsys_solver);
        end

        %% update_settings
        function update_settings(this, varargin)
            % UPDATE_SETTINGS update the current solver settings structure
            %
            %   update_settings('setting1', val1, 'setting2', val2, ...)
            %   update_settings(settings_struct)

            if this.use_matlab_backend
                ml_update_settings_impl(this, varargin{:});
                return;
            end
            newSettings = validate_settings(this, false, varargin{:});
            osqp_mex('update_settings', this.objectHandle, newSettings);
        end

        %% get_dimensions
        function [n, m] = get_dimensions(this)
            % GET_DIMENSIONS get the number of variables and constraints
            if this.use_matlab_backend
                n = this.ml_n;
                m = this.ml_m;
                return;
            end
            [n, m] = osqp_mex('get_dimensions', this.objectHandle);
        end

        %% has_capability
        function out = has_capability(this, cap_name)
            % HAS_CAPABILITY check whether a specific capability is available
            %   tf = HAS_CAPABILITY('OSQP_CAPABILITY_CODEGEN')
            if this.use_matlab_backend
                out = false;
                return;
            end
            cap_val = osqp.constant(cap_name);
            cap_flags = osqp.capabilities();
            out = bitand(uint32(cap_flags), uint32(cap_val)) ~= 0;
        end

        %% setup
        function varargout = setup(this, varargin)
            % SETUP configure solver with problem data
            %
            %   setup(P, q, A, l, u, options)

            nargs = length(varargin);
            assert(nargs >= 5, 'Incorrect number of inputs');
            [P, q, A, l, u] = deal(varargin{1:5});

            % Check if matlab_ldl is requested in setup options
            if ~this.use_matlab_backend && nargs > 5
                opts = varargin(6:end);
                if isstruct(opts{1})
                    s = opts{1};
                    if isfield(s, 'linear_solver') && strcmpi(s.linear_solver, 'matlab_ldl')
                        this.use_matlab_backend = true;
                    end
                else
                    for k = 1:2:numel(opts)
                        if ischar(opts{k}) && strcmpi(opts{k}, 'linear_solver') ...
                                && k+1 <= numel(opts) && ischar(opts{k+1}) ...
                                && strcmpi(opts{k+1}, 'matlab_ldl')
                            this.use_matlab_backend = true;
                        end
                    end
                end
            end

            if this.use_matlab_backend
                ml_setup_impl(this, P, q, A, l, u, varargin{6:end});
                if nargout > 0
                    varargout{1} = [];
                end
                return;
            end

            % Get number of variables n
            if isempty(P)
                if ~isempty(q)
                    n = length(q);
                elseif ~isempty(A)
                    n = size(A, 2);
                else
                    error('The problem does not have any variables');
                end
            else
                n = size(P, 1);
            end

            % Get number of constraints m
            if isempty(A)
                m = 0;
            else
                m = size(A, 1);
                assert(size(A, 2) == n, 'Incorrect dimension of A');
            end

            % Create sparse matrices and full vectors if empty
            if isempty(P)
                P = sparse(n, n);
            else
                P = sparse(P);
            end
            if ~istriu(P)
                P = triu(P);
            end
            if isempty(q)
                q = zeros(n, 1);
            else
                q = full(q(:));
            end

            if (isempty(A) && (~isempty(l) || ~isempty(u))) || ...
               (~isempty(A) && (isempty(l) && isempty(u)))
                error('A must be supplied together with at least one bound l or u');
            end

            if ~isempty(A) && isempty(l)
                l = -Inf(m, 1);
            end
            if ~isempty(A) && isempty(u)
                u = Inf(m, 1);
            end

            if isempty(A)
                A = sparse(m, n);
                l = -Inf(m, 1);
                u = Inf(m, 1);
            else
                l = full(l(:));
                u = full(u(:));
                A = sparse(A);
            end

            % Check dimensions
            assert(length(q) == n, 'Incorrect dimension of q');
            assert(length(l) == m, 'Incorrect dimension of l');
            assert(length(u) == m, 'Incorrect dimension of u');

            % Clamp infinity values
            u = min(u, osqp.constant('OSQP_INFTY'));
            l = max(l, -osqp.constant('OSQP_INFTY'));

            % Build settings
            theSettings = validate_settings(this, true, varargin{6:end});

            [varargout{1:nargout}] = osqp_mex('setup', this.objectHandle, ...
                n, m, P, q, A, l, u, theSettings);
        end

        %% solve
        function varargout = solve(this, varargin)
            % SOLVE solve the QP
            %
            %   results = solve()
            %   Returns a struct with fields: x, y, prim_inf_cert,
            %   dual_inf_cert, info

            nargoutchk(0, 1);

            if this.use_matlab_backend
                out = ml_solve_impl(this);
            else
                [out.x, out.y, out.prim_inf_cert, out.dual_inf_cert, out.info] = ...
                    osqp_mex('solve', this.objectHandle);
            end
            if nargout
                varargout{1} = out;
            end
        end

        %% update
        function update(this, varargin)
            % UPDATE modify problem vectors and/or matrices
            %
            %   update('q', q_new, 'l', l_new, 'u', u_new)
            %   update('Px', Px_new, 'Px_idx', Px_idx_new)
            %   update('Ax', Ax_new, 'Ax_idx', Ax_idx_new)
            %   update(struct_with_fields)

            if this.use_matlab_backend
                ml_update_impl(this, varargin{:});
                return;
            end

            allowedFields = {'q', 'l', 'u', 'Px', 'Px_idx', 'Ax', 'Ax_idx'};

            if isempty(varargin)
                return;
            elseif isscalar(varargin)
                if ~isstruct(varargin{1})
                    error('Single input should be a structure with new problem data');
                end
                newData = varargin{1};
            else
                newData = struct(varargin{:});
            end

            % Check for unknown fields
            newFields = fieldnames(newData);
            badIdx = find(~ismember(newFields, allowedFields));
            if ~isempty(badIdx)
                error('Unrecognized input field ''%s''', newFields{badIdx(1)});
            end

            % Extract fields
            try q = double(full(newData.q(:))); catch, q = []; end
            try l = double(full(newData.l(:))); catch, l = []; end
            try u = double(full(newData.u(:))); catch, u = []; end
            try Px = double(full(newData.Px(:))); catch, Px = []; end
            try Px_idx = double(full(newData.Px_idx(:))); catch, Px_idx = []; end
            try Ax = double(full(newData.Ax(:))); catch, Ax = []; end
            try Ax_idx = double(full(newData.Ax_idx(:))); catch, Ax_idx = []; end

            [n, m] = get_dimensions(this);

            assert(isempty(q) || length(q) == n, 'input ''q'' is the wrong size');
            assert(isempty(l) || length(l) == m, 'input ''l'' is the wrong size');
            assert(isempty(u) || length(u) == m, 'input ''u'' is the wrong size');
            assert(isempty(Px) || isempty(Px_idx) || length(Px) == length(Px_idx), ...
                'inputs ''Px'' and ''Px_idx'' must be the same size');
            assert(isempty(Ax) || isempty(Ax_idx) || length(Ax) == length(Ax_idx), ...
                'inputs ''Ax'' and ''Ax_idx'' must be the same size');

            % Convert to 0-based indexing
            if ~isempty(Px_idx)
                Px_idx = Px_idx - 1;
            end
            if ~isempty(Ax_idx)
                Ax_idx = Ax_idx - 1;
            end

            % Clamp infinity values
            if ~isempty(u)
                u = min(u, osqp.constant('OSQP_INFTY'));
            end
            if ~isempty(l)
                l = max(l, -osqp.constant('OSQP_INFTY'));
            end

            osqp_mex('update', this.objectHandle, ...
                q, l, u, Px, Px_idx, length(Px), Ax, Ax_idx, length(Ax));
        end

        %% warm_start
        function warm_start(this, varargin)
            % WARM_START warm start primal and/or dual variables
            %
            %   warm_start('x', x, 'y', y)
            %   warm_start('x', x)
            %   warm_start('y', y)

            [n, m] = get_dimensions(this);

            allowedFields = {'x', 'y'};

            if isempty(varargin)
                return;
            elseif isscalar(varargin)
                if ~isstruct(varargin{1})
                    error('Single input should be a structure');
                end
                newData = varargin{1};
            else
                newData = struct(varargin{:});
            end

            newFields = fieldnames(newData);
            badIdx = find(~ismember(newFields, allowedFields));
            if ~isempty(badIdx)
                error('Unrecognized input field ''%s''', newFields{badIdx(1)});
            end

            try x = double(full(newData.x(:))); catch, x = []; end
            try y = double(full(newData.y(:))); catch, y = []; end

            assert(isempty(x) || length(x) == n, 'input ''x'' is the wrong size');
            assert(isempty(y) || length(y) == m, 'input ''y'' is the wrong size');

            if this.use_matlab_backend
                x_updated = false;
                y_updated = false;
                if ~isempty(x)
                    this.ml_x = x;
                    this.ml_z = this.ml_A * x;
                    x_updated = true;
                end
                if ~isempty(y)
                    this.ml_y = y;
                    y_updated = true;
                end
                if x_updated && ~y_updated
                    this.ml_y = zeros(m, 1);
                elseif ~x_updated && y_updated
                    this.ml_x = zeros(n, 1);
                    this.ml_z = zeros(m, 1);
                end
                return;
            end

            osqp_mex('warm_start', this.objectHandle, x, y);
        end

        %% cold_start
        function cold_start(this)
            % COLD_START reset solver iterates to zero
            if this.use_matlab_backend
                this.ml_x = zeros(this.ml_n, 1);
                this.ml_z = zeros(this.ml_m, 1);
                this.ml_y = zeros(this.ml_m, 1);
                return;
            end
            osqp_mex('cold_start', this.objectHandle);
        end

        %% codegen
        function codegen(this, target_dir, varargin)
            % CODEGEN generate C code for the parametric problem
            if this.use_matlab_backend
                error('OSQP:codegen', 'Code generation is not supported with the matlab_ldl backend.');
            end
            %
            %   codegen(target_dir, options)
            %
            % Options (name-value pairs):
            %   'prefix'             - prefix for generated code (default '')
            %   'parameters'         - 'vectors' or 'matrices' (default 'vectors')
            %   'force_rewrite'      - overwrite target_dir (default false)
            %   'float_type'         - use single precision (default false)
            %   'printing_enable'    - enable printing (default false)
            %   'profiling_enable'   - enable profiling (default false)
            %   'interrupt_enable'   - enable interrupt checking (default false)
            %   'derivatives_enable' - enable derivatives (default false)

            p = inputParser;
            addRequired(p, 'target_dir', @ischar);
            addParameter(p, 'prefix', '', @ischar);
            addParameter(p, 'parameters', 'vectors', ...
                @(x) ischar(validatestring(x, {'vectors', 'matrices'})));
            addParameter(p, 'force_rewrite', false, @islogical);
            addParameter(p, 'float_type', false, @islogical);
            addParameter(p, 'printing_enable', false, @islogical);
            addParameter(p, 'profiling_enable', false, @islogical);
            addParameter(p, 'interrupt_enable', false, @islogical);
            addParameter(p, 'derivatives_enable', false, @islogical);

            parse(p, target_dir, varargin{:});

            if strcmp(p.Results.parameters, 'vectors')
                embedded_mode = 1;
            else
                embedded_mode = 2;
            end

            % Handle existing directory
            if exist(target_dir, 'dir')
                if p.Results.force_rewrite
                    rmdir(target_dir, 's');
                else
                    while true
                        prompt = sprintf('Directory "%s" already exists. Replace? y/n [y]: ', target_dir);
                        str = input(prompt, 's');
                        if any(strcmpi(str, {'', 'y'}))
                            rmdir(target_dir, 's');
                            break;
                        elseif strcmpi(str, 'n')
                            return;
                        end
                    end
                end
            end

            mkdir(target_dir);

            % Build codegen defines struct
            defines = struct();
            defines.embedded_mode     = embedded_mode;
            defines.float_type        = p.Results.float_type;
            defines.printing_enable   = p.Results.printing_enable;
            defines.profiling_enable  = p.Results.profiling_enable;
            defines.interrupt_enable  = p.Results.interrupt_enable;
            defines.derivatives_enable = p.Results.derivatives_enable;

            % Call native osqp_codegen to generate workspace files
            osqp_mex('codegen', this.objectHandle, ...
                target_dir, p.Results.prefix, defines);

            % Copy OSQP codegen source files into target directory
            [osqp_path, ~, ~] = fileparts(which('osqp.m'));
            osqp_root = fileparts(osqp_path);  % go up from src/ to project root
            codegen_src_dir = fullfile(osqp_root, 'build', 'codegen_src');
            if ~exist(codegen_src_dir, 'dir')
                % Try alternate location under _deps
                codegen_src_dir = fullfile(osqp_root, 'build', '_deps', ...
                    'osqp-build', 'codegen_src');
            end
            if ~exist(codegen_src_dir, 'dir')
                % Try installed location
                codegen_src_dir = fullfile(osqp_path, 'codegen', 'sources');
            end
            if exist(codegen_src_dir, 'dir')
                % Copy src/ and inc/ directories
                if exist(fullfile(codegen_src_dir, 'src'), 'dir')
                    copyfile(fullfile(codegen_src_dir, 'src'), ...
                        fullfile(target_dir, 'src'));
                end
                if exist(fullfile(codegen_src_dir, 'inc'), 'dir')
                    copyfile(fullfile(codegen_src_dir, 'inc'), ...
                        fullfile(target_dir, 'inc'));
                end
                if exist(fullfile(codegen_src_dir, 'Makefile'), 'file')
                    copyfile(fullfile(codegen_src_dir, 'Makefile'), ...
                        fullfile(target_dir, 'Makefile'));
                end
            end
        end

        %% adjoint_derivative_compute
        function adjoint_derivative_compute(this, dx, dy)
            % ADJOINT_DERIVATIVE_COMPUTE compute adjoint derivatives
            if this.use_matlab_backend
                error('OSQP:adjoint', 'Adjoint derivatives are not supported with the matlab_ldl backend.');
            end
            %
            %   adjoint_derivative_compute(dx, dy)
            %
            %   dx: n-vector of cost perturbations in primal
            %   dy: m-vector of cost perturbations in dual

            [n, m] = get_dimensions(this);
            assert(length(dx) == n, 'dx has wrong size');
            assert(length(dy) == m, 'dy has wrong size');

            osqp_mex('adjoint_derivative_compute', this.objectHandle, ...
                double(full(dx(:))), double(full(dy(:))));
        end

        %% adjoint_derivative_get_mat
        function [dP, dA] = adjoint_derivative_get_mat(this)
            % ADJOINT_DERIVATIVE_GET_MAT retrieve derivative matrices
            if this.use_matlab_backend
                error('OSQP:adjoint', 'Adjoint derivatives are not supported with the matlab_ldl backend.');
            end
            %
            %   [dP, dA] = adjoint_derivative_get_mat()
            %
            %   dP: n-by-n sparse matrix (upper triangular)
            %   dA: m-by-n sparse matrix

            [dP, dA] = osqp_mex('adjoint_derivative_get_mat', this.objectHandle);
        end

        %% adjoint_derivative_get_vec
        function [dq, dl, du] = adjoint_derivative_get_vec(this)
            % ADJOINT_DERIVATIVE_GET_VEC retrieve derivative vectors
            if this.use_matlab_backend
                error('OSQP:adjoint', 'Adjoint derivatives are not supported with the matlab_ldl backend.');
            end
            %
            %   [dq, dl, du] = adjoint_derivative_get_vec()

            [dq, dl, du] = osqp_mex('adjoint_derivative_get_vec', this.objectHandle);
        end

    end

    % =====================================================================
    % Static helper for MATLAB backend default settings
    % =====================================================================
    methods (Static, Access = private)
        function s = matlab_ldl_default_settings()
            s.rho                      = 0.1;
            s.sigma                    = 1e-6;
            s.scaling                  = 10;
            s.adaptive_rho             = true;
            s.adaptive_rho_interval    = 0;
            s.adaptive_rho_tolerance   = 5;
            s.adaptive_rho_fraction    = 0.4;
            s.max_iter                 = 4000;
            s.eps_abs                  = 1e-3;
            s.eps_rel                  = 1e-3;
            s.eps_prim_inf             = 1e-4;
            s.eps_dual_inf             = 1e-4;
            s.alpha                    = 1.6;
            s.delta                    = 1e-6;
            s.polishing                = false;
            s.polish_refine_iter       = 3;
            s.verbose                  = true;
            s.scaled_termination       = false;
            s.check_termination        = 25;
            s.warm_starting            = true;
            s.time_limit               = 0;
            s.linsys_solver            = 'matlab_ldl';
            s.rho_is_vec               = 0;
            s.linear_solver            = 'matlab_ldl';
        end
    end

    % =====================================================================
    % Private methods implementing the pure-MATLAB backend
    % =====================================================================
    methods (Access = private)

        function ml_setup_impl(this, P, q, A, l, u, varargin)
            t_start = tic;

            % Parse settings
            this.ml_settings = osqp.matlab_ldl_default_settings();
            if ~isempty(varargin)
                if isstruct(varargin{1})
                    s = varargin{1};
                    fnames = fieldnames(s);
                    for k = 1:numel(fnames)
                        if isfield(this.ml_settings, fnames{k})
                            this.ml_settings.(fnames{k}) = s.(fnames{k});
                        end
                    end
                else
                    if mod(numel(varargin), 2) == 0
                        for k = 1:2:numel(varargin)
                            if isfield(this.ml_settings, varargin{k})
                                this.ml_settings.(varargin{k}) = varargin{k+1};
                            end
                        end
                    end
                end
            end

            % Determine dimensions
            if isempty(P)
                if ~isempty(q)
                    this.ml_n = numel(q);
                elseif ~isempty(A)
                    this.ml_n = size(A, 2);
                else
                    error('OSQP:setup', 'Problem has no variables.');
                end
            else
                this.ml_n = size(P, 1);
            end
            if isempty(A)
                this.ml_m = 0;
            else
                this.ml_m = size(A, 1);
            end
            n = this.ml_n;
            m = this.ml_m;

            % Default missing data
            if isempty(P), P = sparse(n, n); end
            if isempty(q), q = zeros(n, 1); end
            if isempty(A), A = sparse(0, n); l = zeros(0, 1); u = zeros(0, 1); end
            if isempty(l), l = -inf(m, 1); end
            if isempty(u), u = inf(m, 1); end

            % Validate
            q = double(full(q(:)));
            l = double(full(l(:)));
            u = double(full(u(:)));
            P = sparse(P);
            if ~istriu(P), P = triu(P); end

            this.ml_P_triu = P;
            this.ml_q = q;
            this.ml_A = sparse(A);
            this.ml_l = max(l, -osqp.OSQP_INFTY);
            this.ml_u = min(u,  osqp.OSQP_INFTY);

            % Convexity check
            Pfull = this.ml_P_triu + this.ml_P_triu' - diag(diag(this.ml_P_triu));
            this.ml_non_convex = false;
            if n > 0 && nnz(Pfull) > 0
                thresh = 1e-7;
                [~, flag] = chol(Pfull + thresh * speye(n), 'lower');
                if flag ~= 0
                    Psig = Pfull + this.ml_settings.sigma * speye(n);
                    [~, flag2] = chol(Psig + thresh * speye(n), 'lower');
                    if flag2 ~= 0
                        error('OSQP:NonConvex', ...
                            'P is non-convex and sigma is too small.');
                    end
                    this.ml_non_convex = true;
                end
            end

            % Initialize iterates
            this.ml_x = zeros(n, 1);
            this.ml_z = zeros(m, 1);
            this.ml_y = zeros(m, 1);

            % Classify constraints
            this.ml_constr_type = osqp.ml_classify_constraints(this.ml_l, this.ml_u);

            % Scale problem
            [this.ml_scl, this.ml_Ps, this.ml_qs, this.ml_As, this.ml_ls, this.ml_us] = ...
                osqp.ml_scale_problem(this.ml_P_triu, this.ml_q, this.ml_A, ...
                this.ml_l, this.ml_u, this.ml_settings);

            % Build rho vector
            [this.ml_rho_vec, this.ml_rho_inv_vec] = osqp.ml_make_rho_vec( ...
                this.ml_constr_type, this.ml_settings.rho, m);

            % Factorize KKT
            this.ml_kkt_factor = osqp.ml_factorize_kkt( ...
                this.ml_Ps, this.ml_As, this.ml_rho_vec, ...
                this.ml_settings.sigma, n, m);

            this.ml_isSetup = true;
            this.ml_setup_time = toc(t_start);
            this.ml_update_time = 0;
            this.ml_polish_time = 0;
        end

        function results = ml_solve_impl(this)
            if ~this.ml_isSetup
                error('OSQP:solve', 'Solver not set up.');
            end

            n = this.ml_n;
            m = this.ml_m;
            s = this.ml_settings;

            results = osqp.ml_empty_results(n, m);

            % Non-convex check
            if this.ml_non_convex
                results.info.status_val = osqp.ML_STATUS_NON_CONVEX;
                results.info.status     = 'non_convex';
                results.info.obj_val    = nan;
                return;
            end

            t_solve = tic;

            % Apply scaling to initial iterates
            if s.warm_starting
                xs = this.ml_scl.Dinv .* this.ml_x;
                zs = this.ml_scl.E    .* this.ml_z;
                ys = this.ml_scl.Einv .* this.ml_y * this.ml_scl.c;
            else
                xs = zeros(n, 1);
                zs = zeros(m, 1);
                ys = zeros(m, 1);
            end

            status_val = osqp.ML_STATUS_UNSOLVED;
            rho_updates = 0;
            adaptive_rho_interval = 0;

            if s.adaptive_rho && s.adaptive_rho_interval > 0
                adaptive_rho_interval = s.adaptive_rho_interval;
            end

            prim_res = inf;
            dual_res = inf;

            for iter = 1:s.max_iter
                xs_prev = xs;
                zs_prev = zs;
                ys_prev = ys;

                % Step 1: Solve KKT for x_tilde, z_tilde
                rhs_x = s.sigma * xs_prev - this.ml_qs;
                if m > 0
                    rhs_z = zs_prev - this.ml_rho_inv_vec .* ys;
                    rhs   = [rhs_x; rhs_z];
                else
                    rhs = rhs_x;
                end

                sol = this.ml_kkt_factor \ rhs;
                xtilde = sol(1:n);

                if m > 0
                    ztilde = zs_prev + this.ml_rho_inv_vec .* (sol(n+1:end) - ys);
                else
                    ztilde = zeros(0, 1);
                end

                % Step 2: Relaxation + projection
                xs = s.alpha * xtilde + (1 - s.alpha) * xs_prev;
                if m > 0
                    zs_relaxed = s.alpha * ztilde + (1 - s.alpha) * zs_prev;
                    zs = min(max(zs_relaxed + this.ml_rho_inv_vec .* ys, ...
                        this.ml_ls), this.ml_us);
                end

                % Step 3: Dual update
                if m > 0
                    ys = ys + this.ml_rho_vec .* ( ...
                        s.alpha * ztilde + (1 - s.alpha) * zs_prev - zs);
                end

                % Convergence check
                if s.check_termination > 0 && mod(iter, s.check_termination) == 0
                    [prim_res, dual_res] = osqp.ml_compute_residuals( ...
                        xs, zs, ys, this.ml_Ps, this.ml_qs, this.ml_As, n, m);
                    [converged, status_val] = osqp.ml_check_convergence( ...
                        prim_res, dual_res, xs, zs, ys, xs_prev, ys_prev, ...
                        this.ml_Ps, this.ml_qs, this.ml_As, this.ml_ls, this.ml_us, ...
                        n, m, s);
                    if converged
                        break;
                    end
                end

                % Time limit
                if s.time_limit > 0 && toc(t_solve) >= s.time_limit
                    status_val = osqp.ML_STATUS_TIME_LIMIT_REACHED;
                    break;
                end

                % Adaptive rho
                if s.adaptive_rho
                    if adaptive_rho_interval == 0
                        if iter == 1
                            t_iter_start = tic;
                        elseif iter == 2
                            t_one_iter = toc(t_iter_start);
                            t_setup = max(this.ml_setup_time, 1e-10);
                            if t_one_iter > 0
                                adaptive_rho_interval = max(1, ...
                                    round(s.adaptive_rho_fraction * t_setup / t_one_iter));
                            else
                                adaptive_rho_interval = 25;
                            end
                        end
                    end
                    if adaptive_rho_interval > 0 && mod(iter, adaptive_rho_interval) == 0
                        [pr2, dr2] = osqp.ml_compute_residuals( ...
                            xs, zs, ys, this.ml_Ps, this.ml_qs, this.ml_As, n, m);
                        new_rho = osqp.ml_compute_new_rho( ...
                            pr2, dr2, this.ml_settings.rho, s.adaptive_rho_tolerance);
                        if new_rho ~= this.ml_settings.rho
                            this.ml_settings.rho = new_rho;
                            [this.ml_rho_vec, this.ml_rho_inv_vec] = osqp.ml_make_rho_vec( ...
                                this.ml_constr_type, new_rho, m);
                            this.ml_kkt_factor = osqp.ml_factorize_kkt( ...
                                this.ml_Ps, this.ml_As, this.ml_rho_vec, ...
                                s.sigma, n, m);
                            rho_updates = rho_updates + 1;
                        end
                    end
                end
            end % for iter

            this.ml_solve_time = toc(t_solve);

            % Final convergence test if loop ended without check
            if status_val == osqp.ML_STATUS_UNSOLVED
                [prim_res, dual_res] = osqp.ml_compute_residuals( ...
                    xs, zs, ys, this.ml_Ps, this.ml_qs, this.ml_As, n, m);
                [~, status_val] = osqp.ml_check_convergence( ...
                    prim_res, dual_res, xs, zs, ys, xs_prev, ys_prev, ...
                    this.ml_Ps, this.ml_qs, this.ml_As, this.ml_ls, this.ml_us, ...
                    n, m, s);
                if status_val == osqp.ML_STATUS_UNSOLVED
                    status_val = osqp.ML_STATUS_MAX_ITER_REACHED;
                end
            end

            % Unscale solution
            x_out = this.ml_scl.D    .* xs;
            z_out = this.ml_scl.Einv .* zs;
            y_out = this.ml_scl.E    .* ys * this.ml_scl.cinv;

            % Store for warm starting
            this.ml_x = x_out;
            this.ml_z = z_out;
            this.ml_y = y_out;

            % Fill results
            results.info.iter        = iter;
            results.info.status_val  = status_val;
            results.info.status      = osqp.ml_status_str(status_val);
            results.info.rho_updates = rho_updates;
            results.info.rho_estimate = this.ml_settings.rho;
            results.info.setup_time  = this.ml_setup_time;
            results.info.solve_time  = this.ml_solve_time;
            results.info.update_time = this.ml_update_time;

            solution_present = ismember(status_val, [ ...
                osqp.ML_STATUS_SOLVED, ...
                osqp.ML_STATUS_SOLVED_INACCURATE, ...
                osqp.ML_STATUS_MAX_ITER_REACHED]);

            if solution_present
                results.x = x_out;
                results.y = y_out;
                results.prim_inf_cert = nan(m, 1);
                results.dual_inf_cert = nan(n, 1);
                Pfull = this.ml_P_triu + this.ml_P_triu' - diag(diag(this.ml_P_triu));
                results.info.obj_val = 0.5 * (x_out' * (Pfull * x_out)) + ...
                    this.ml_q' * x_out;
                results.info.pri_res = norm(this.ml_A * x_out - z_out, inf);
                results.info.dua_res = norm(Pfull * x_out + this.ml_q + this.ml_A' * y_out, inf);
            elseif status_val == osqp.ML_STATUS_PRIMAL_INFEASIBLE || ...
                    status_val == osqp.ML_STATUS_PRIMAL_INFEASIBLE_INACCURATE
                results.x = nan(n, 1);
                results.y = nan(m, 1);
                delta_y = this.ml_scl.E .* (ys - ys_prev);
                if norm(delta_y, inf) > 0
                    delta_y = delta_y / norm(delta_y, inf);
                end
                results.prim_inf_cert = delta_y;
                results.dual_inf_cert = nan(n, 1);
                results.info.obj_val = inf;
                results.info.pri_res = inf;
                results.info.dua_res = inf;
            elseif status_val == osqp.ML_STATUS_DUAL_INFEASIBLE || ...
                    status_val == osqp.ML_STATUS_DUAL_INFEASIBLE_INACCURATE
                results.x = nan(n, 1);
                results.y = nan(m, 1);
                results.prim_inf_cert = nan(m, 1);
                delta_x = this.ml_scl.D .* (xs - xs_prev);
                if norm(delta_x, inf) > 0
                    delta_x = delta_x / norm(delta_x, inf);
                end
                results.dual_inf_cert = delta_x;
                results.info.obj_val = -inf;
                results.info.pri_res = inf;
                results.info.dua_res = inf;
            else
                results.x = nan(n, 1);
                results.y = nan(m, 1);
                results.prim_inf_cert = nan(m, 1);
                results.dual_inf_cert = nan(n, 1);
                results.info.obj_val = nan;
                results.info.pri_res = inf;
                results.info.dua_res = inf;
            end

            % Polishing
            results.info.status_polish = 0;
            if s.polishing && solution_present
                t_polish = tic;
                results = osqp.ml_polish_solution(results, ...
                    this.ml_P_triu, this.ml_q, this.ml_A, this.ml_l, this.ml_u, s);
                this.ml_polish_time = toc(t_polish);
            end
            results.info.polish_time = this.ml_polish_time;
            results.info.run_time = this.ml_setup_time + this.ml_update_time + ...
                this.ml_solve_time + this.ml_polish_time;

            if s.verbose
                fprintf('OSQP(matlab_ldl): %s | iter=%d | obj=%.4e | pri=%.2e | dua=%.2e\n', ...
                    results.info.status, results.info.iter, ...
                    results.info.obj_val, results.info.pri_res, results.info.dua_res);
            end
        end

        function ml_update_impl(this, varargin)
            if ~this.ml_isSetup
                error('OSQP:update', 'Call setup() before update().');
            end
            t_update = tic;

            allowedFields = {'q', 'l', 'u', 'Px', 'Px_idx', 'Ax', 'Ax_idx'};

            if isempty(varargin)
                return;
            elseif isscalar(varargin)
                if isstruct(varargin{1})
                    newData = varargin{1};
                else
                    error('Single input should be a structure with new problem data');
                end
            else
                newData = struct(varargin{:});
            end

            newFields = fieldnames(newData);
            badIdx = find(~ismember(newFields, allowedFields));
            if ~isempty(badIdx)
                error('Unrecognized input field ''%s''', newFields{badIdx(1)});
            end

            refactor = false;

            % q
            if isfield(newData, 'q')
                q_new = double(full(newData.q(:)));
                if numel(q_new) ~= this.ml_n
                    error('OSQP:update', 'q must have length n=%d', this.ml_n);
                end
                this.ml_q = q_new;
            end

            % l and u
            if isfield(newData, 'l') || isfield(newData, 'u')
                if isfield(newData, 'l')
                    l_new = double(full(newData.l(:)));
                    if numel(l_new) ~= this.ml_m
                        error('OSQP:update', 'l must have length m=%d', this.ml_m);
                    end
                    this.ml_l = max(l_new, -osqp.OSQP_INFTY);
                end
                if isfield(newData, 'u')
                    u_new = double(full(newData.u(:)));
                    if numel(u_new) ~= this.ml_m
                        error('OSQP:update', 'u must have length m=%d', this.ml_m);
                    end
                    this.ml_u = min(u_new, osqp.OSQP_INFTY);
                end
                this.ml_constr_type = osqp.ml_classify_constraints(this.ml_l, this.ml_u);
                [this.ml_rho_vec, this.ml_rho_inv_vec] = osqp.ml_make_rho_vec( ...
                    this.ml_constr_type, this.ml_settings.rho, this.ml_m);
                refactor = true;
            end

            % Px
            if isfield(newData, 'Px')
                Px_new = double(full(newData.Px(:)));
                [ri, ci, ~] = find(this.ml_P_triu);
                nz_orig = nnz(this.ml_P_triu);
                if isfield(newData, 'Px_idx')
                    idx = double(full(newData.Px_idx(:)));
                    if numel(Px_new) ~= numel(idx)
                        error('OSQP:update', 'Px and Px_idx must have same length');
                    end
                    old_vals = nonzeros(this.ml_P_triu);
                    old_vals(idx) = Px_new;
                    this.ml_P_triu = sparse(ri, ci, old_vals, this.ml_n, this.ml_n);
                else
                    if numel(Px_new) ~= nz_orig
                        error('OSQP:update', 'Px must have %d elements', nz_orig);
                    end
                    this.ml_P_triu = sparse(ri, ci, Px_new, this.ml_n, this.ml_n);
                end
                refactor = true;
            end

            % Ax
            if isfield(newData, 'Ax')
                Ax_new = double(full(newData.Ax(:)));
                if isfield(newData, 'Ax_idx')
                    idx = double(full(newData.Ax_idx(:)));
                    if numel(Ax_new) ~= numel(idx)
                        error('OSQP:update', 'Ax and Ax_idx must have same length');
                    end
                    [ri, ci, vals] = find(this.ml_A);
                    vals(idx) = Ax_new;
                    this.ml_A = sparse(ri, ci, vals, this.ml_m, this.ml_n);
                else
                    nz_orig = nnz(this.ml_A);
                    if numel(Ax_new) ~= nz_orig
                        error('OSQP:update', 'Ax must have %d elements', nz_orig);
                    end
                    [ri, ci] = find(this.ml_A);
                    this.ml_A = sparse(ri, ci, Ax_new, this.ml_m, this.ml_n);
                end
                refactor = true;
            end

            % Rescale and refactorize if needed
            if refactor
                [this.ml_scl, this.ml_Ps, this.ml_qs, this.ml_As, this.ml_ls, this.ml_us] = ...
                    osqp.ml_scale_problem(this.ml_P_triu, this.ml_q, this.ml_A, ...
                    this.ml_l, this.ml_u, this.ml_settings);
                this.ml_kkt_factor = osqp.ml_factorize_kkt( ...
                    this.ml_Ps, this.ml_As, this.ml_rho_vec, ...
                    this.ml_settings.sigma, this.ml_n, this.ml_m);
            else
                % Only q changed; recompute scaled q
                this.ml_qs = this.ml_scl.c * (this.ml_scl.D .* this.ml_q);
            end

            this.ml_update_time = this.ml_update_time + toc(t_update);
        end

        function ml_update_settings_impl(this, varargin)
            if ~this.ml_isSetup
                error('OSQP:update_settings', 'Call setup() before update_settings().');
            end

            updatable = {'max_iter','eps_abs','eps_rel','eps_prim_inf','eps_dual_inf', ...
                'time_limit','rho','alpha','delta','polishing','polish_refine_iter', ...
                'verbose','check_termination','warm_starting','scaled_termination'};

            rho_updated = false;
            if isscalar(varargin) && isstruct(varargin{1})
                s = varargin{1};
                fnames = fieldnames(s);
                for k = 1:numel(fnames)
                    if ~ismember(fnames{k}, updatable)
                        error('OSQP:update_settings', ...
                            'Setting ''%s'' cannot be updated.', fnames{k});
                    end
                    this.ml_settings.(fnames{k}) = s.(fnames{k});
                end
                rho_updated = isfield(s, 'rho');
            else
                if mod(numel(varargin), 2) ~= 0
                    error('OSQP:update_settings', 'Arguments must be name/value pairs.');
                end
                for k = 1:2:numel(varargin)
                    name = varargin{k};
                    if ~ismember(name, updatable)
                        error('OSQP:update_settings', ...
                            'Setting ''%s'' cannot be updated.', name);
                    end
                    this.ml_settings.(name) = varargin{k+1};
                    if strcmp(name, 'rho'), rho_updated = true; end
                end
            end

            if rho_updated
                [this.ml_rho_vec, this.ml_rho_inv_vec] = osqp.ml_make_rho_vec( ...
                    this.ml_constr_type, this.ml_settings.rho, this.ml_m);
                this.ml_kkt_factor = osqp.ml_factorize_kkt( ...
                    this.ml_Ps, this.ml_As, this.ml_rho_vec, ...
                    this.ml_settings.sigma, this.ml_n, this.ml_m);
            end
        end

    end % methods (private)

    % =====================================================================
    % Static helpers for the MATLAB backend
    % =====================================================================
    methods (Static, Access = private)

        function str = ml_status_str(val)
            switch val
                case  4,  str = 'dual_infeasible_inaccurate';
                case  3,  str = 'primal_infeasible_inaccurate';
                case  2,  str = 'solved_inaccurate';
                case  1,  str = 'solved';
                case -2,  str = 'maximum_iterations_reached';
                case -3,  str = 'primal_infeasible';
                case -4,  str = 'dual_infeasible';
                case -6,  str = 'time_limit_reached';
                case -7,  str = 'non_convex';
                otherwise, str = 'unsolved';
            end
        end

        function r = ml_empty_results(n, m)
            r.x = nan(n, 1);
            r.y = nan(m, 1);
            r.prim_inf_cert = nan(m, 1);
            r.dual_inf_cert = nan(n, 1);
            r.info.iter          = 0;
            r.info.status        = 'unsolved';
            r.info.status_val    = osqp.ML_STATUS_UNSOLVED;
            r.info.status_polish = 0;
            r.info.obj_val       = nan;
            r.info.pri_res       = nan;
            r.info.dua_res       = nan;
            r.info.setup_time    = 0;
            r.info.solve_time    = 0;
            r.info.update_time   = 0;
            r.info.polish_time   = 0;
            r.info.run_time      = 0;
            r.info.rho_updates   = 0;
            r.info.rho_estimate  = 0;
        end

        function ctype = ml_classify_constraints(l, u)
            m = numel(l);
            ctype = zeros(m, 1);
            ctype(l == u) = 1;
        end

        function [rho_vec, rho_inv_vec] = ml_make_rho_vec(ctype, rho, m)
            if m == 0
                rho_vec = zeros(0, 1); rho_inv_vec = zeros(0, 1); return;
            end
            rho_vec = rho * ones(m, 1);
            rho_vec(ctype == 1) = rho * 1e3;
            rho_inv_vec = 1 ./ rho_vec;
        end

        function [scl, Ps, qs, As, ls, us] = ml_scale_problem(P_triu, q, A, l, u, settings)
            n = size(P_triu, 1);
            m = size(A, 1);
            num_iter = settings.scaling;

            D = ones(n, 1);
            E = ones(m, 1);
            c = 1.0;

            if num_iter == 0
                Ps = P_triu; qs = q; As = A; ls = l; us = u;
                scl.D = D; scl.E = E; scl.Dinv = D; scl.Einv = E;
                scl.c = 1; scl.cinv = 1;
                return;
            end

            Pfull = P_triu + P_triu' - diag(diag(P_triu));

            for iter_s = 1:num_iter
                if n > 0
                    Dsp = spdiags(D, 0, n, n);
                    PDE = Dsp * Pfull * Dsp;
                    AD  = A * Dsp;
                    combined = max(full(max(abs(PDE)))', full(max(abs(AD)))');
                    pos = combined > 0;
                    D(pos) = D(pos) ./ sqrt(combined(pos));
                end
                if m > 0
                    Dsp = spdiags(D, 0, n, n);
                    Esp = spdiags(E, 0, m, m);
                    EAD = Esp * A * Dsp;
                    EAD_norms = full(max(abs(EAD), [], 2));
                    pos = EAD_norms > 0;
                    E(pos) = E(pos) ./ sqrt(EAD_norms(pos));
                end
                DPD = D' .* diag(sparse(Pfull))' .* D';
                mean_norm = mean(abs(DPD));
                q_norms   = abs(q .* D);
                cost_scale = max(mean_norm, mean(q_norms));
                if cost_scale > 0
                    c = c / sqrt(cost_scale);
                end
            end

            D = max(D, 1e-6);
            E = max(E, 1e-6);
            c = max(c, 1e-6);
            Dinv = 1 ./ D;
            Einv = 1 ./ E;
            cinv  = 1 / c;

            Dsp = diag(sparse(D));
            Esp = diag(sparse(E));
            Pfull_scaled = c * (Dsp * Pfull * Dsp);
            Ps = triu(Pfull_scaled);
            qs = c * (D .* q);
            As = Esp * A * Dsp;
            ls = max(E .* l, -osqp.OSQP_INFTY);
            us = min(E .* u,  osqp.OSQP_INFTY);

            scl.D = D; scl.E = E; scl.Dinv = Dinv; scl.Einv = Einv;
            scl.c = c; scl.cinv = cinv;
        end

        function F = ml_factorize_kkt(P_triu, A, rho_vec, sigma, n, m)
            Pfull = P_triu + P_triu' - diag(diag(P_triu));
            Ktl   = Pfull + sigma * speye(n);
            if m > 0
                Kbr = -diag(sparse(1 ./ rho_vec));
                K   = [Ktl, A'; A, Kbr];
            else
                K = Ktl;
            end
            K = (K + K') / 2;
            F = MATLABLDLSolver(K);
        end

        function [prim_res, dual_res] = ml_compute_residuals(xs, zs, ys, Ps, qs, As, ~, m)
            if m > 0
                prim_res = norm(As * xs - zs, inf);
            else
                prim_res = 0;
            end
            Pfull = Ps + Ps' - diag(diag(Ps));
            if m > 0
                dual_res = norm(Pfull * xs + qs + As' * ys, inf);
            else
                dual_res = norm(Pfull * xs + qs, inf);
            end
        end

        function [converged, status_val] = ml_check_convergence( ...
                prim_res, dual_res, xs, zs, ys, xs_prev, ys_prev, ...
                Ps, qs, As, ls, us, ~, m, s)
            status_val = osqp.ML_STATUS_UNSOLVED;
            converged  = false;

            Pfull = Ps + Ps' - diag(diag(Ps));

            % Adaptive tolerances
            if m > 0
                Axs_norm = norm(As * xs, inf);
                zs_norm  = norm(zs, inf);
                eps_prim = s.eps_abs + s.eps_rel * max(Axs_norm, zs_norm);
            else
                eps_prim = s.eps_abs;
            end
            Pxs_norm  = norm(Pfull * xs, inf);
            Atys_norm = 0;
            if m > 0
                Atys_norm = norm(As' * ys, inf);
            end
            qs_norm  = norm(qs, inf);
            eps_dual = s.eps_abs + s.eps_rel * max([Pxs_norm; Atys_norm; qs_norm]);

            % Optimality
            if prim_res <= eps_prim && dual_res <= eps_dual
                status_val = osqp.ML_STATUS_SOLVED;
                converged  = true;
                return;
            end

            % Primal infeasibility (certificate: delta_y)
            if m > 0
                dy      = ys - ys_prev;
                norm_dy = norm(dy, inf);
                if norm_dy > s.eps_prim_inf
                    dy_n  = dy / norm_dy;
                    ATdy  = norm(As' * dy_n, inf);
                    lu_ok = osqp.ml_primal_inf_check(dy_n, ls, us, s.eps_prim_inf);
                    if ATdy <= s.eps_prim_inf && lu_ok
                        status_val = osqp.ML_STATUS_PRIMAL_INFEASIBLE;
                        converged  = true;
                        return;
                    end
                end
            end

            % Dual infeasibility (certificate: delta_x)
            dx      = xs - xs_prev;
            norm_dx = norm(dx, inf);
            if norm_dx > s.eps_dual_inf
                dx_n = dx / norm_dx;
                Pdx  = norm(Pfull * dx_n, inf);
                qdx  = qs' * dx_n;
                if Pdx <= s.eps_dual_inf && qdx < -s.eps_dual_inf
                    if m > 0
                        Adx = As * dx_n;
                        ok  = osqp.ml_dual_inf_check(Adx, ls, us, s.eps_dual_inf);
                    else
                        ok = true;
                    end
                    if ok
                        status_val = osqp.ML_STATUS_DUAL_INFEASIBLE;
                        converged  = true;
                        return;
                    end
                end
            end
        end

        function ok = ml_primal_inf_check(v, l, u, eps)
            INF = osqp.OSQP_INFTY;
            vpos = max(v, 0);
            vneg = min(v, 0);
            val  = 0;
            for i = 1:numel(l)
                if abs(u(i)) < INF, val = val + u(i) * vpos(i); end
                if abs(l(i)) < INF, val = val + l(i) * vneg(i); end
            end
            ok = val < -eps;
        end

        function ok = ml_dual_inf_check(Adx, l, u, eps)
            INF = osqp.OSQP_INFTY;
            ok = true;
            for i = 1:numel(l)
                li = l(i); ui = u(i);
                li_fin = abs(li) < INF;
                ui_fin = abs(ui) < INF;
                ai = Adx(i);
                if li_fin && ui_fin
                    if abs(ai) > eps, ok = false; return; end
                elseif li_fin && ~ui_fin
                    if ai < -eps, ok = false; return; end
                elseif ~li_fin && ui_fin
                    if ai > eps, ok = false; return; end
                end
            end
        end

        function new_rho = ml_compute_new_rho(prim_res, dual_res, rho, tol)
            if prim_res == 0 || dual_res == 0
                new_rho = rho; return;
            end
            ratio = sqrt(prim_res / dual_res);
            new_rho_candidate = rho * ratio;
            if new_rho_candidate > tol * rho || new_rho_candidate < rho / tol
                new_rho = max(min(new_rho_candidate, 1e6), 1e-6);
            else
                new_rho = rho;
            end
        end

        function results = ml_polish_solution(results, P_triu, q, A, l, u, s)
            n = numel(q);
            m = size(A, 1);
            x = results.x;
            y = results.y;
            Pfull = P_triu + P_triu' - diag(diag(P_triu));
            INF = osqp.OSQP_INFTY;

            Ax = A * x;
            tol_act = max(s.delta, s.eps_abs);
            act_l = (y < -s.delta & abs(l) < INF) | ...
                (abs(l) < INF & abs(Ax - l) <= tol_act);
            act_u = (y >  s.delta & abs(u) < INF) | ...
                (abs(u) < INF & abs(Ax - u) <= tol_act);
            active_idx = find(act_l | act_u);
            n_act = numel(active_idx);

            if n_act == 0
                try
                    Preg = Pfull + s.delta * speye(n);
                    x_pol = Preg \ (-q);
                    Ax_pol = A * x_pol;
                    feas_l = all(Ax_pol >= l - tol_act | abs(l) >= INF);
                    feas_u = all(Ax_pol <= u + tol_act | abs(u) >= INF);
                    if feas_l && feas_u
                        obj_pol = 0.5 * x_pol' * Pfull * x_pol + q' * x_pol;
                        obj_cur = 0.5 * x' * Pfull * x + q' * x;
                        if obj_pol < obj_cur + abs(obj_cur) * 1e-6
                            results.x = x_pol;
                            results.y = zeros(m, 1);
                            results.info.obj_val = obj_pol;
                            results.info.status_polish = 1;
                        else
                            results.info.status_polish = -1;
                        end
                    else
                        results.info.status_polish = -1;
                    end
                catch
                    results.info.status_polish = -1;
                end
                return;
            end

            A_act = A(active_idx, :);
            b_act = zeros(n_act, 1);
            for k = 1:n_act
                i = active_idx(k);
                if act_l(i) && act_u(i)
                    b_act(k) = (l(i) + u(i)) / 2;
                elseif act_l(i)
                    b_act(k) = l(i);
                else
                    b_act(k) = u(i);
                end
            end

            try
                Kreg = [Pfull + s.delta*speye(n), A_act'; ...
                    A_act, -s.delta*speye(n_act)];
                rhs0 = [-q; b_act];
                sol0 = Kreg \ rhs0;
                x_pol = sol0(1:n);
                y_act = sol0(n+1:end);

                for refine = 1:s.polish_refine_iter
                    res = rhs0 - Kreg * [x_pol; y_act];
                    delta_sol = Kreg \ res;
                    x_pol = x_pol + delta_sol(1:n);
                    y_act = y_act + delta_sol(n+1:end);
                end

                y_pol = zeros(m, 1);
                y_pol(active_idx) = y_act;

                Ax_pol = A * x_pol;
                feas_l = all(Ax_pol >= l - tol_act | abs(l) >= INF);
                feas_u = all(Ax_pol <= u + tol_act | abs(u) >= INF);
                if feas_l && feas_u
                    obj_pol = 0.5 * x_pol' * Pfull * x_pol + q' * x_pol;
                    % Compute polished residuals and compare with ADMM
                    pri_res_pol = norm(Ax_pol - min(max(Ax_pol, l), u), inf);
                    dua_res_pol = norm(Pfull * x_pol + q + A' * y_pol, inf);
                    pri_res_cur = results.info.pri_res;
                    dua_res_cur = results.info.dua_res;
                    % Accept only if polished residuals are not worse
                    if max(pri_res_pol, dua_res_pol) <= ...
                            max(pri_res_cur, dua_res_cur) * 10 + 1e-10
                        results.x = x_pol;
                        results.y = y_pol;
                        results.info.obj_val = obj_pol;
                        results.info.pri_res = pri_res_pol;
                        results.info.dua_res = dua_res_pol;
                        results.info.status_polish = 1;
                    else
                        results.info.status_polish = -1;
                    end
                else
                    results.info.status_polish = -1;
                end
            catch
                results.info.status_polish = -1;
            end
        end

    end % methods (static, private)

end


%% ---- Local helper functions ----

function currentSettings = validate_settings(this, isInit, varargin)
    % Build a settings struct by merging defaults/current with overrides

    % Setup-only fields (cannot change after setup)
    setupOnlyFields = {'scaling', 'rho_is_vec', 'sigma', ...
        'adaptive_rho', 'adaptive_rho_interval', ...
        'adaptive_rho_fraction', 'adaptive_rho_tolerance'};

    if isInit
        currentSettings = osqp_mex('default_settings', this.objectHandle);
    else
        currentSettings = osqp_mex('current_settings', this.objectHandle);
    end

    if isempty(varargin)
        return;
    end

    if isstruct(varargin{1})
        newSettings = varargin{1};
        assert(isscalar(varargin), 'Too many input arguments');
    else
        newSettings = struct(varargin{:});
    end

    currentFields = fieldnames(currentSettings);
    newFields = fieldnames(newSettings);

    % Check for unknown fields
    badIdx = find(~ismember(newFields, currentFields));
    if ~isempty(badIdx)
        error('Unrecognized solver setting ''%s''', newFields{badIdx(1)});
    end

    % Convert linsys_solver string to integer
    if ismember('linsys_solver', newFields)
        if ischar(newSettings.linsys_solver)
            newSettings.linsys_solver = string_to_linsys_solver(newSettings.linsys_solver);
        end
    end

    % Check setup-only fields
    if ~isInit
        badIdx = find(ismember(newFields, setupOnlyFields));
        for i = badIdx(:)'
            if ~isequal(newSettings.(newFields{i}), currentSettings.(newFields{i}))
                error('Setting ''%s'' can only be changed at setup.', newFields{i});
            end
        end
    end

    % Merge new settings into current
    for i = 1:length(newFields)
        currentSettings.(newFields{i}) = double(newSettings.(newFields{i}));
    end
end


function str = linsys_solver_to_string(val)
    % Convert linsys_solver integer to string
    DIRECT   = osqp.constant('OSQP_DIRECT_SOLVER');
    INDIRECT = osqp.constant('OSQP_INDIRECT_SOLVER');
    switch val
        case DIRECT
            str = 'direct';
        case INDIRECT
            str = 'indirect';
        otherwise
            str = 'unknown';
    end
end


function val = string_to_linsys_solver(str)
    str = lower(str);
    switch str
        case 'direct'
            val = osqp.constant('OSQP_DIRECT_SOLVER');
        case 'indirect'
            val = osqp.constant('OSQP_INDIRECT_SOLVER');
        case ''
            val = osqp.constant('OSQP_DIRECT_SOLVER');
        otherwise
            warning('Linear system solver ''%s'' not recognized. Using direct solver.', str);
            val = osqp.constant('OSQP_DIRECT_SOLVER');
    end
end


