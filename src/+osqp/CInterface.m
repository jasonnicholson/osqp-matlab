classdef CInterface < handle
    % OSQP.CINTERFACE  MEX wrapper around the OSQP C library.
    %
    %   Thin handle class that owns an opaque C workspace pointer and
    %   delegates every operation to osqp_mex.
    %
    %   solver = osqp.CInterface();
    %   solver.setup(P, q, A, l, u, 'max_iter', 1000);
    %   results = solver.solve();

    properties (SetAccess = private, Hidden = true)
        objectHandle   % Opaque pointer to the C OSQPWorkspace
    end

    % =================================================================
    %  Construction / destruction
    % =================================================================
    methods
        function obj = CInterface(varargin)
            obj.objectHandle = osqp_mex('new', varargin{:});
        end

        function delete(obj)
            if ~isempty(obj.objectHandle)
                osqp_mex('delete', obj.objectHandle);
            end
        end
    end

    % =================================================================
    %  Static methods
    % =================================================================
    methods (Static)
        function out = isMexAvailable()
            % ISMEXAVAILABLE  Check whether osqp_mex is on the path.
            out = exist('osqp_mex', 'file') == 3;
        end

        function out = default_settings()
            % DEFAULT_SETTINGS  Get the default settings structure from C.
            out = osqp_mex('default_settings', 'static');
            out.linsys_solver = osqp.CInterface.linsys_solver_to_string(out.linsys_solver);
        end

        function out = constant(constant_name)
            % CONSTANT  Return an OSQP constant by name.
            out = osqp_mex('constant', 'static', constant_name);
        end

        function out = version()
            % VERSION  Return the linked OSQP library version string.
            out = osqp_mex('version', 'static');
        end

        function out = capabilities()
            % CAPABILITIES  Return capability bitmask from the C library.
            out = osqp_mex('capabilities', 'static');
        end
    end

    % =================================================================
    %  Public API
    % =================================================================
    methods
        function setup(obj, P, q, A, l, u, varargin)
            % SETUP  Configure solver with problem data.
            %
            %   setup(P, q, A, l, u)
            %   setup(P, q, A, l, u, 'Name', Value, ...)
            %   setup(P, q, A, l, u, settings_struct)

            INFTY = 1e30;

            % Determine dimensions
            if isempty(P)
                if ~isempty(q)
                    n = numel(q);
                elseif ~isempty(A)
                    n = size(A, 2);
                else
                    error('OSQP:setup', 'Problem has no variables.');
                end
            else
                n = size(P, 1);
            end
            if isempty(A), m = 0; else, m = size(A, 1); end

            % Defaults
            if isempty(P), P = sparse(n, n); else, P = sparse(P); end
            if ~istriu(P), P = triu(P); end
            if isempty(q), q = zeros(n, 1); else, q = full(q(:)); end

            if isempty(A)
                A = sparse(0, n); l = zeros(0, 1); u = zeros(0, 1);
            else
                A = sparse(A);
                if isempty(l), l = -inf(m, 1); end
                if isempty(u), u =  inf(m, 1); end
                l = full(l(:)); u = full(u(:));
            end

            l = max(l, -INFTY);
            u = min(u,  INFTY);

            % Build settings struct for MEX
            theSettings = obj.validateSettings(true, varargin{:});

            osqp_mex('setup', obj.objectHandle, ...
                n, m, P, q, A, l, u, theSettings);
        end

        function out = solve(obj)
            % SOLVE  Solve the QP. Returns struct with x, y, info.
            [out.x, out.y, out.prim_inf_cert, out.dual_inf_cert, out.info] = ...
                osqp_mex('solve', obj.objectHandle);
        end

        function update(obj, varargin)
            % UPDATE  Modify problem vectors and/or matrices.
            %
            %   update('q', q_new, 'l', l_new, 'u', u_new)
            %   update('Px', Px_new, 'Px_idx', Px_idx_new)
            %   update('Ax', Ax_new, 'Ax_idx', Ax_idx_new)
            %   update(struct_with_fields)

            INFTY = 1e30;
            allowedFields = {'q', 'l', 'u', 'Px', 'Px_idx', 'Ax', 'Ax_idx'};

            if isempty(varargin), return; end

            if isscalar(varargin)
                if ~isstruct(varargin{1})
                    error('Single input should be a structure with new problem data');
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

            try q = double(full(newData.q(:)));         catch, q = [];      end
            try l = double(full(newData.l(:)));         catch, l = [];      end
            try u = double(full(newData.u(:)));         catch, u = [];      end
            try Px = double(full(newData.Px(:)));       catch, Px = [];     end
            try Px_idx = double(full(newData.Px_idx(:)));catch, Px_idx = []; end
            try Ax = double(full(newData.Ax(:)));       catch, Ax = [];     end
            try Ax_idx = double(full(newData.Ax_idx(:)));catch, Ax_idx = []; end

            [n, m] = obj.get_dimensions();

            assert(isempty(q) || length(q) == n, 'input ''q'' is the wrong size');
            assert(isempty(l) || length(l) == m, 'input ''l'' is the wrong size');
            assert(isempty(u) || length(u) == m, 'input ''u'' is the wrong size');
            assert(isempty(Px) || isempty(Px_idx) || length(Px) == length(Px_idx), ...
                'inputs ''Px'' and ''Px_idx'' must be the same size');
            assert(isempty(Ax) || isempty(Ax_idx) || length(Ax) == length(Ax_idx), ...
                'inputs ''Ax'' and ''Ax_idx'' must be the same size');

            % Convert to 0-based indexing for MEX
            if ~isempty(Px_idx), Px_idx = Px_idx - 1; end
            if ~isempty(Ax_idx), Ax_idx = Ax_idx - 1; end

            if ~isempty(u), u = min(u,  INFTY); end
            if ~isempty(l), l = max(l, -INFTY); end

            osqp_mex('update', obj.objectHandle, ...
                q, l, u, Px, Px_idx, length(Px), Ax, Ax_idx, length(Ax));
        end

        function warm_start(obj, varargin)
            % WARM_START  Set warm starting variables x and/or y.
            [n, m] = obj.get_dimensions();

            if isempty(varargin), return; end

            if isscalar(varargin) && isstruct(varargin{1})
                newData = varargin{1};
            else
                newData = struct(varargin{:});
            end

            try x = double(full(newData.x(:))); catch, x = []; end
            try y = double(full(newData.y(:))); catch, y = []; end

            assert(isempty(x) || length(x) == n, 'input ''x'' is the wrong size');
            assert(isempty(y) || length(y) == m, 'input ''y'' is the wrong size');

            osqp_mex('warm_start', obj.objectHandle, x, y);
        end

        function cold_start(obj)
            % COLD_START  Reset solver iterates to zero.
            osqp_mex('cold_start', obj.objectHandle);
        end

        function out = current_settings(obj)
            % CURRENT_SETTINGS  Get the current settings structure from C.
            out = osqp_mex('current_settings', obj.objectHandle);
            out.linsys_solver = osqp.CInterface.linsys_solver_to_string(out.linsys_solver);
        end

        function update_settings(obj, varargin)
            % UPDATE_SETTINGS  Update solver settings in the C workspace.
            newSettings = obj.validateSettings(false, varargin{:});
            osqp_mex('update_settings', obj.objectHandle, newSettings);
        end

        function [n, m] = get_dimensions(obj)
            % GET_DIMENSIONS  Get the number of variables and constraints.
            [n, m] = osqp_mex('get_dimensions', obj.objectHandle);
        end

        function out = has_capability(~, cap_name)
            % HAS_CAPABILITY  Check whether a specific capability is available.
            cap_val   = osqp.CInterface.constant(cap_name);
            cap_flags = osqp.CInterface.capabilities();
            out = bitand(uint32(cap_flags), uint32(cap_val)) ~= 0;
        end

        function codegen(obj, target_dir, varargin)
            % CODEGEN  Generate embeddable C code for the problem.
            %
            %   codegen(target_dir, 'Name', Value, ...)
            %
            % Options:
            %   'prefix'             - '' (default)
            %   'parameters'         - 'vectors' | 'matrices'
            %   'force_rewrite'      - false
            %   'float_type'         - false
            %   'printing_enable'    - false
            %   'profiling_enable'   - false
            %   'interrupt_enable'   - false
            %   'derivatives_enable' - false

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

            if exist(target_dir, 'dir')
                if p.Results.force_rewrite
                    rmdir(target_dir, 's');
                else
                    while true
                        prompt = sprintf( ...
                            'Directory "%s" already exists. Replace? y/n [y]: ', target_dir);
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

            defines.embedded_mode      = embedded_mode;
            defines.float_type         = p.Results.float_type;
            defines.printing_enable    = p.Results.printing_enable;
            defines.profiling_enable   = p.Results.profiling_enable;
            defines.interrupt_enable   = p.Results.interrupt_enable;
            defines.derivatives_enable = p.Results.derivatives_enable;

            % Ensure trailing filesep so C osqp_codegen() places files
            % inside target_dir (it concatenates dir + prefix + name).
            codegen_dir = target_dir;
            if codegen_dir(end) ~= filesep
                codegen_dir(end+1) = filesep;
            end

            osqp_mex('codegen', obj.objectHandle, ...
                codegen_dir, p.Results.prefix, defines);

            % Copy source files
            [osqp_path, ~, ~] = fileparts(which('osqp.m'));
            osqp_root = fileparts(osqp_path);
            codegen_src_dir = fullfile(osqp_root, 'build', 'codegen_src');
            if ~exist(codegen_src_dir, 'dir')
                codegen_src_dir = fullfile(osqp_root, 'build', '_deps', ...
                    'osqp-build', 'codegen_src');
            end
            if ~exist(codegen_src_dir, 'dir')
                codegen_src_dir = fullfile(osqp_path, 'codegen', 'sources');
            end
            if exist(codegen_src_dir, 'dir')
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

        function adjoint_derivative_compute(obj, dx, dy)
            % ADJOINT_DERIVATIVE_COMPUTE  Compute adjoint derivatives.
            [n, m] = obj.get_dimensions();
            assert(length(dx) == n, 'dx has wrong size');
            assert(length(dy) == m, 'dy has wrong size');
            osqp_mex('adjoint_derivative_compute', obj.objectHandle, ...
                double(full(dx(:))), double(full(dy(:))));
        end

        function [dP, dA] = adjoint_derivative_get_mat(obj)
            % ADJOINT_DERIVATIVE_GET_MAT  Retrieve derivative matrices.
            [dP, dA] = osqp_mex('adjoint_derivative_get_mat', obj.objectHandle);
        end

        function [dq, dl, du] = adjoint_derivative_get_vec(obj)
            % ADJOINT_DERIVATIVE_GET_VEC  Retrieve derivative vectors.
            [dq, dl, du] = osqp_mex('adjoint_derivative_get_vec', obj.objectHandle);
        end
    end

    % =================================================================
    %  Private helpers
    % =================================================================
    methods (Access = private)
        function theSettings = validateSettings(~, isSetup, varargin)
            % VALIDATESETTINGS  Parse name/value or struct settings for MEX.
            %
            % At setup time every field is accepted.
            % After setup only updatable fields are accepted.

            if isempty(varargin)
                theSettings = struct();
                return;
            end

            if isstruct(varargin{1})
                theSettings = varargin{1};
            else
                theSettings = struct(varargin{:});
            end

            % Convert string linsys_solver to integer for MEX
            if isfield(theSettings, 'linsys_solver')
                theSettings.linsys_solver = ...
                    osqp.CInterface.string_to_linsys_solver(theSettings.linsys_solver);
            end

            % Strip MATLAB-only fields
            matlabOnly = {'linear_solver'};
            for k = 1:numel(matlabOnly)
                if isfield(theSettings, matlabOnly{k})
                    theSettings = rmfield(theSettings, matlabOnly{k});
                end
            end

            if ~isSetup
                % After setup, only runtime-updatable fields are allowed.
                setupOnly = { ...
                    'linsys_solver', 'scaling', 'adaptive_rho', ...
                    'adaptive_rho_interval', 'adaptive_rho_fraction', ...
                    'polishing', 'verbose', ...
                    'scaled_termination', 'rho_is_vec'};
                fnames = fieldnames(theSettings);
                for k = 1:numel(fnames)
                    if ismember(fnames{k}, setupOnly)
                        error('OSQP:update_settings', ...
                            'Setting ''%s'' cannot be updated after setup.', fnames{k});
                    end
                end
            end
        end
    end

    methods (Static, Access = private)
        function str = linsys_solver_to_string(val)
            switch val
                case 0, str = 'qdldl';
                case 1, str = 'mkl pardiso';
                case 2, str = 'cuda pcg';
                otherwise, str = 'unknown';
            end
        end

        function val = string_to_linsys_solver(str)
            switch lower(str)
                case 'qdldl',      val = 0;
                case 'mkl pardiso', val = 1;
                case 'cuda pcg',   val = 2;
                otherwise
                    error('OSQP:settings', 'Unknown linsys_solver: %s', str);
            end
        end
    end
end
