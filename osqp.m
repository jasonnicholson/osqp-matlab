classdef osqp < handle
    % OSQP interface class for OSQP solver v1.0.0
    %
    % This class provides a complete interface to the C implementation
    % of the OSQP solver.
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

    methods(Static)

        function out = default_settings()
            % DEFAULT_SETTINGS get the default solver settings structure
            out = osqp_mex('default_settings', 'static');
            out.linsys_solver = linsys_solver_to_string(out.linsys_solver);
        end

        function out = constant(constant_name)
            % CONSTANT Return solver constant
            %   C = CONSTANT(CONSTANT_NAME) return constant called CONSTANT_NAME
            out = osqp_mex('constant', 'static', constant_name);
        end

        function out = version()
            % VERSION Return OSQP version string
            out = osqp_mex('version', 'static');
        end

        function out = capabilities()
            % CAPABILITIES Return capability bitmask
            out = osqp_mex('capabilities', 'static');
        end

    end

    methods

        %% Constructor
        function this = osqp(varargin)
            % OSQP Construct OSQP solver class
            this.objectHandle = osqp_mex('new', varargin{:});
        end

        %% Destructor
        function delete(this)
            % DELETE Destroy OSQP solver class
            osqp_mex('delete', this.objectHandle);
        end

        %% current_settings
        function out = current_settings(this)
            % CURRENT_SETTINGS get the current solver settings structure
            out = osqp_mex('current_settings', this.objectHandle);
            out.linsys_solver = linsys_solver_to_string(out.linsys_solver);
        end

        %% update_settings
        function update_settings(this, varargin)
            % UPDATE_SETTINGS update the current solver settings structure
            %
            %   update_settings('setting1', val1, 'setting2', val2, ...)
            %   update_settings(settings_struct)

            newSettings = validate_settings(this, false, varargin{:});
            osqp_mex('update_settings', this.objectHandle, newSettings);
        end

        %% get_dimensions
        function [n, m] = get_dimensions(this)
            % GET_DIMENSIONS get the number of variables and constraints
            [n, m] = osqp_mex('get_dimensions', this.objectHandle);
        end

        %% has_capability
        function out = has_capability(this, cap_name)
            % HAS_CAPABILITY check whether a specific capability is available
            %   tf = HAS_CAPABILITY('OSQP_CAPABILITY_CODEGEN')
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
            [out.x, out.y, out.prim_inf_cert, out.dual_inf_cert, out.info] = ...
                osqp_mex('solve', this.objectHandle);
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

            allowedFields = {'q', 'l', 'u', 'Px', 'Px_idx', 'Ax', 'Ax_idx'};

            if isempty(varargin)
                return;
            elseif length(varargin) == 1
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
            elseif length(varargin) == 1
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

            osqp_mex('warm_start', this.objectHandle, x, y);
        end

        %% cold_start
        function cold_start(this)
            % COLD_START reset solver iterates to zero
            osqp_mex('cold_start', this.objectHandle);
        end

        %% codegen
        function codegen(this, target_dir, varargin)
            % CODEGEN generate C code for the parametric problem
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
            codegen_src_dir = fullfile(osqp_path, 'build', '_deps', ...
                'osqp-build', 'codegen_src');
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
            %
            %   [dq, dl, du] = adjoint_derivative_get_vec()

            [dq, dl, du] = osqp_mex('adjoint_derivative_get_vec', this.objectHandle);
        end

    end
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
        assert(length(varargin) == 1, 'Too many input arguments');
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


