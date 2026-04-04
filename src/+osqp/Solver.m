classdef Solver < handle
    % OSQP.SOLVER  Pure-MATLAB ADMM solver for quadratic programs.
    %
    %   Solves:  minimize   0.5 * x' * P * x + q' * x
    %            subject to l <= A * x <= u
    %
    %   solver = osqp.Solver();
    %   solver.setup(P, q, A, l, u);
    %   results = solver.solve();
    %
    %   solver = osqp.Solver();
    %   solver.setup(P, q, A, l, u, 'max_iter', 1000, 'verbose', false);

    % -----------------------------------------------------------------
    %  Problem data (unscaled)
    % -----------------------------------------------------------------
    properties (Access = private)
        P_triu      % upper triangular Hessian (sparse, n×n)
        q_          % linear cost (n×1)
        A_          % constraint matrix (sparse, m×n)
        l_          % lower bounds (m×1)
        u_          % upper bounds (m×1)
        n_ = 0      % number of variables
        m_ = 0      % number of constraints
    end

    % -----------------------------------------------------------------
    %  Settings and state
    % -----------------------------------------------------------------
    properties (Access = private)
        opts        % osqp.Options

        % Scaling
        scl         % struct: D, E, Dinv, Einv, c, cinv

        % Scaled problem
        Ps          % scaled P (upper tri)
        qs          % scaled q
        As          % scaled A
        ls          % scaled l
        us          % scaled u

        % ADMM iterates
        x_          % primal (n×1)
        z_          % slack  (m×1)
        y_          % dual   (m×1)

        % Rho vectors
        rho_vec     % per-constraint rho (m×1)
        rho_inv_vec % 1./rho_vec (m×1)
        constr_type % 0=inequality, 1=equality (m×1)

        % Linear system solver
        kkt_factor  % osqp.LinearSolver subclass

        % Flags / timing
        isSetup = false
        non_convex = false
        setup_time = 0
        solve_time = 0
        update_time = 0
        polish_time = 0

        % Cached algorithm constants (populated once during setup)
        C_          % struct: all constants from osqp.constant + derived values
    end

    % =====================================================================
    %  Public API
    % =====================================================================
    methods
        function setup(obj, P, q, A, l, u, varargin)
            % SETUP  Configure solver with problem data.
            %
            %   setup(P, q, A, l, u)
            %   setup(P, q, A, l, u, 'Name', Value, ...)
            %   setup(P, q, A, l, u, opts_struct)
            %   setup(P, q, A, l, u, osqp.Options)
            t_start = tic;

            % --- Parse settings ---
            obj.opts = osqp.Options();
            if ~isempty(varargin)
                if isa(varargin{1}, 'osqp.Options')
                    obj.opts = varargin{1};
                elseif isstruct(varargin{1})
                    obj.opts = osqp.Options.fromStruct(varargin{1});
                else
                    s = struct(varargin{:});
                    obj.opts = osqp.Options.fromStruct(s);
                end
            end

            % --- Cache algorithm constants ---
            obj.cacheConstants();

            % --- Determine dimensions ---
            if isempty(P)
                if ~isempty(q)
                    obj.n_ = numel(q);
                elseif ~isempty(A)
                    obj.n_ = size(A, 2);
                else
                    error('OSQP:setup', 'Problem has no variables.');
                end
            else
                obj.n_ = size(P, 1);
            end
            if isempty(A)
                obj.m_ = 0;
            else
                obj.m_ = size(A, 1);
            end
            n = obj.n_;
            m = obj.m_;

            % --- Default missing data ---
            if isempty(P), P = sparse(n, n); end
            if isempty(q), q = zeros(n, 1); end
            if isempty(A), A = sparse(0, n); l = zeros(0, 1); u = zeros(0, 1); end
            if isempty(l), l = -inf(m, 1); end
            if isempty(u), u = inf(m, 1); end

            % --- Validate and store ---
            q = double(full(q(:)));
            l = double(full(l(:)));
            u = double(full(u(:)));
            P = sparse(P);
            if ~istriu(P), P = triu(P); end

            obj.P_triu = P;
            obj.q_ = q;
            obj.A_ = sparse(A);
            obj.l_ = max(l, -obj.C_.OSQP_INFTY);
            obj.u_ = min(u,  obj.C_.OSQP_INFTY);

            % --- Convexity check (matching C: LDL-based inertia) ---
            %   C code factors KKT (with P+sigma*I) using QDLDL and
            %   verifies that D has exactly n positive entries.  We match
            %   by checking P+sigma*I is PD, then using LDL on P alone
            %   to detect indefiniteness (PSD is accepted, unlike Cholesky).
            Pfull = obj.P_triu + obj.P_triu' - diag(diag(obj.P_triu));
            obj.non_convex = false;
            if n > 0 && nnz(Pfull) > 0
                % Step 1: P + sigma*I must be PD (required for KKT to work)
                Psig = Pfull + obj.opts.sigma * speye(n);
                [~, flag_sig] = chol(Psig, 'lower');
                if flag_sig ~= 0
                    error('OSQP:NonConvex', ...
                        'P is non-convex and sigma is too small.');
                end
                % Step 2: Check if P itself is PSD using LDL factorization.
                % LDL handles PSD (zero eigenvalues) without failure,
                % unlike Cholesky which requires strict positive definiteness.
                obj.non_convex = ~osqp.Solver.isLDLpsd(Pfull);
            end

            % --- Initialize iterates ---
            obj.x_ = zeros(n, 1);
            obj.z_ = zeros(m, 1);
            obj.y_ = zeros(m, 1);

            % --- Classify constraints ---
            obj.constr_type = osqp.Solver.classifyConstraints(obj.l_, obj.u_, obj.C_);

            % --- Scale problem ---
            [obj.scl, obj.Ps, obj.qs, obj.As, obj.ls, obj.us] = ...
                osqp.Solver.scaleProblem(obj.P_triu, obj.q_, obj.A_, ...
                obj.l_, obj.u_, obj.opts, obj.C_);

            % --- Build rho vector ---
            [obj.rho_vec, obj.rho_inv_vec] = osqp.Solver.makeRhoVec( ...
                obj.constr_type, obj.opts.rho, m, obj.C_);

            % --- Factorize KKT ---
            obj.factorizeKKT();

            obj.isSetup = true;
            obj.setup_time = toc(t_start);
            obj.update_time = 0;
            obj.polish_time = 0;
        end

        function results = solve(obj)
            % SOLVE  Solve the QP. Returns struct with x, y, info.
            if ~obj.isSetup
                error('OSQP:solve', 'Solver not set up. Call setup() first.');
            end

            n = obj.n_;
            m = obj.m_;
            s = obj.opts;
            C = obj.C_;

            results = osqp.Solver.emptyResults(n, m, C);

            % Non-convex check
            if obj.non_convex
                results.info.status_val = C.OSQP_NON_CONVEX;
                results.info.status     = 'non_convex';
                results.info.obj_val    = nan;
                return;
            end

            t_solve = tic;

            % Apply scaling to initial iterates
            if s.warm_starting
                xs = obj.scl.Dinv .* obj.x_;
                zs = obj.scl.E    .* obj.z_;
                ys = obj.scl.Einv .* obj.y_ * obj.scl.c;
            else
                xs = zeros(n, 1);
                zs = zeros(m, 1);
                ys = zeros(m, 1);
            end

            status_val = C.OSQP_UNSOLVED;
            rho_updates = 0;
            adaptive_rho_interval = 0;

            if s.adaptive_rho && s.adaptive_rho_interval > 0
                adaptive_rho_interval = s.adaptive_rho_interval;
            end

            iter = 0; %#ok<NASGU> default for max_iter=0

            % ============================================================
            % ADMM iteration loop (C: osqp_solve main loop)
            % ============================================================
            for iter = 1:s.max_iter
                xs_prev = xs;
                zs_prev = zs;
                ys_prev = ys;

                % Step 1: KKT solve for (x_tilde, z_tilde)
                %   Solves [P+sigma*I  A'; A  -diag(1/rho)] * [xt; nut]
                %        = [sigma*x_prev - q; z_prev - rho_inv.*y]
                %   C: compute_rhs() + linsys_solver->solve()
                rhs_x = s.sigma * xs_prev - obj.qs;
                if m > 0
                    rhs_z = zs_prev - obj.rho_inv_vec .* ys;
                    rhs   = [rhs_x; rhs_z];
                else
                    rhs = rhs_x;
                end

                sol = obj.kkt_factor \ rhs;
                xtilde = sol(1:n);

                % C: update_xz_tilde() — recover z_tilde from KKT auxiliary
                if m > 0
                    ztilde = zs_prev + obj.rho_inv_vec .* (sol(n+1:end) - ys);
                else
                    ztilde = zeros(0, 1);
                end

                % Step 2: Over-relaxation + projection onto [l, u]
                %   x = alpha * xtilde + (1-alpha) * x_prev
                %   z = clip(alpha*ztilde + (1-alpha)*z_prev + rho_inv.*y, l, u)
                %   C: update_x() + update_z()
                xs = s.alpha * xtilde + (1 - s.alpha) * xs_prev;
                if m > 0
                    zs_relaxed = s.alpha * ztilde + (1 - s.alpha) * zs_prev;
                    zs = min(max(zs_relaxed + obj.rho_inv_vec .* ys, ...
                        obj.ls), obj.us);
                end

                % Step 3: Dual variable update
                %   y = y + rho .* (alpha*ztilde + (1-alpha)*z_prev - z)
                %   C: update_y()
                if m > 0
                    ys = ys + obj.rho_vec .* ( ...
                        s.alpha * ztilde + (1 - s.alpha) * zs_prev - zs);
                end

                % Periodic convergence check (C: check_termination)
                if s.check_termination > 0 && mod(iter, s.check_termination) == 0
                    [converged, status_val] = osqp.Solver.checkConvergence( ...
                        xs, zs, ys, xs_prev, ys_prev, ...
                        obj.Ps, obj.qs, obj.As, obj.ls, obj.us, ...
                        n, m, s, obj.scl, C);
                    if converged
                        break;
                    end
                end

                % Time limit (C: time_limit branch)
                if s.time_limit > 0 && toc(t_solve) >= s.time_limit
                    % Try approximate (10x tolerances) check before declaring time limit
                    [conv_tl, sv_tl] = osqp.Solver.checkConvergence( ...
                        xs, zs, ys, xs_prev, ys_prev, ...
                        obj.Ps, obj.qs, obj.As, obj.ls, obj.us, ...
                        n, m, s, obj.scl, C, true);  % approximate = true
                    if conv_tl
                        status_val = sv_tl;
                    else
                        status_val = C.OSQP_TIME_LIMIT_REACHED;
                    end
                    break;
                end

                % Adaptive rho update (C: adapt_rho)
                %   Recomputes rho based on primal/dual residual ratio and
                %   refactorizes KKT if rho changed significantly.
                if s.adaptive_rho
                    if adaptive_rho_interval == 0
                        % v1.0: iteration-count-based interval
                        adaptive_rho_interval = ...
                            C.OSQP_ADAPTIVE_RHO_MULTIPLE_TERMINATION * ...
                            max(s.check_termination, 1);
                    end
                    if adaptive_rho_interval > 0 && mod(iter, adaptive_rho_interval) == 0
                        new_rho = osqp.Solver.computeNewRho( ...
                            xs, zs, ys, obj.Ps, obj.qs, obj.As, ...
                            obj.opts.rho, s.adaptive_rho_tolerance, n, m, C);
                        if new_rho ~= obj.opts.rho
                            obj.opts.rho = new_rho;
                            [obj.rho_vec, obj.rho_inv_vec] = osqp.Solver.makeRhoVec( ...
                                obj.constr_type, new_rho, m, C);
                            obj.factorizeKKT();
                            rho_updates = rho_updates + 1;
                        end
                    end
                end
            end % for iter

            obj.solve_time = toc(t_solve);

            % Post-loop convergence check (C: check_termination after loop)
            %   First exact, then approximate (10x tolerances), then max_iter.
            if status_val == C.OSQP_UNSOLVED
                [~, status_val] = osqp.Solver.checkConvergence( ...
                    xs, zs, ys, xs_prev, ys_prev, ...
                    obj.Ps, obj.qs, obj.As, obj.ls, obj.us, ...
                    n, m, s, obj.scl, C);
                % If still unsolved, try approximate check (10x tolerances)
                % matching C code's post-loop check_termination(solver, 1)
                if status_val == C.OSQP_UNSOLVED
                    [~, status_val] = osqp.Solver.checkConvergence( ...
                        xs, zs, ys, xs_prev, ys_prev, ...
                        obj.Ps, obj.qs, obj.As, obj.ls, obj.us, ...
                        n, m, s, obj.scl, C, true);  % approximate = true
                end
                if status_val == C.OSQP_UNSOLVED
                    status_val = C.OSQP_MAX_ITER_REACHED;
                end
            end

            % Unscale solution back to original problem space
            %   x = D * xs,  z = Einv * zs,  y = E * ys / c
            %   C: unscale_solution()
            x_out = obj.scl.D    .* xs;
            z_out = obj.scl.Einv .* zs;
            y_out = obj.scl.E    .* ys * obj.scl.cinv;

            % Store for warm starting
            obj.x_ = x_out;
            obj.z_ = z_out;
            obj.y_ = y_out;

            % Fill results
            results.info.iter        = iter;
            results.info.status_val  = status_val;
            results.info.status      = osqp.Solver.statusStr(status_val);
            results.info.rho_updates = rho_updates;
            results.info.rho_estimate = obj.opts.rho;
            results.info.setup_time  = obj.setup_time;
            results.info.solve_time  = obj.solve_time;
            results.info.update_time = obj.update_time;

            % Build result struct based on termination status
            %   C: store_solution() and update_info()
            solution_present = ismember(status_val, [ ...
                C.OSQP_SOLVED, ...
                C.OSQP_SOLVED_INACCURATE, ...
                C.OSQP_MAX_ITER_REACHED]);

            if solution_present
                results.x = x_out;
                results.y = y_out;
                results.prim_inf_cert = nan(m, 1);
                results.dual_inf_cert = nan(n, 1);
                Pfull = obj.P_triu + obj.P_triu' - diag(diag(obj.P_triu));
                results.info.obj_val = 0.5 * (x_out' * (Pfull * x_out)) + ...
                    obj.q_' * x_out;
                results.info.pri_res = norm(obj.A_ * x_out - z_out, inf);
                results.info.dua_res = norm(Pfull * x_out + obj.q_ + obj.A_' * y_out, inf);
            elseif status_val == C.OSQP_PRIMAL_INFEASIBLE || ...
                    status_val == C.OSQP_PRIMAL_INFEASIBLE_INACCURATE
                results.x = nan(n, 1);
                results.y = nan(m, 1);
                delta_y = obj.scl.E .* (ys - ys_prev);
                if norm(delta_y, inf) > 0
                    delta_y = delta_y / norm(delta_y, inf);
                end
                results.prim_inf_cert = delta_y;
                results.dual_inf_cert = nan(n, 1);
                results.info.obj_val = inf;
                results.info.pri_res = inf;
                results.info.dua_res = inf;
            elseif status_val == C.OSQP_DUAL_INFEASIBLE || ...
                    status_val == C.OSQP_DUAL_INFEASIBLE_INACCURATE
                results.x = nan(n, 1);
                results.y = nan(m, 1);
                results.prim_inf_cert = nan(m, 1);
                delta_x = obj.scl.D .* (xs - xs_prev);
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

            % Solution polishing (C: polish.c)
            %   Identifies active constraints, forms a reduced KKT system,
            %   solves it with iterative refinement, then applies Moreau
            %   decomposition to ensure complementarity.
            results.info.status_polish = 0;
            if s.polishing && status_val == C.OSQP_SOLVED
                t_polish = tic;
                results = osqp.Solver.polishSolution(results, ...
                    obj.P_triu, obj.q_, obj.A_, obj.l_, obj.u_, s, C);
                obj.polish_time = toc(t_polish);
            end
            results.info.polish_time = obj.polish_time;
            results.info.run_time = obj.setup_time + obj.update_time + ...
                obj.solve_time + obj.polish_time;

            if s.verbose
                fprintf('OSQP(matlab): %s | iter=%d | obj=%.4e | pri=%.2e | dua=%.2e\n', ...
                    results.info.status, results.info.iter, ...
                    results.info.obj_val, results.info.pri_res, results.info.dua_res);
            end
        end

        function update(obj, varargin)
            % UPDATE  Modify problem vectors and/or matrices.
            %
            %   update('q', q_new, 'l', l_new, 'u', u_new)
            %   update('Px', Px_new, 'Px_idx', Px_idx_new)
            %   update('Ax', Ax_new, 'Ax_idx', Ax_idx_new)
            %   update(struct_with_fields)
            if ~obj.isSetup
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
                if numel(q_new) ~= obj.n_
                    error('OSQP:update', 'q must have length n=%d', obj.n_);
                end
                obj.q_ = q_new;
            end

            % l and u
            if isfield(newData, 'l') || isfield(newData, 'u')
                if isfield(newData, 'l')
                    l_new = double(full(newData.l(:)));
                    if numel(l_new) ~= obj.m_
                        error('OSQP:update', 'l must have length m=%d', obj.m_);
                    end
                    obj.l_ = max(l_new, -obj.C_.OSQP_INFTY);
                end
                if isfield(newData, 'u')
                    u_new = double(full(newData.u(:)));
                    if numel(u_new) ~= obj.m_
                        error('OSQP:update', 'u must have length m=%d', obj.m_);
                    end
                    obj.u_ = min(u_new, obj.C_.OSQP_INFTY);
                end
                obj.constr_type = osqp.Solver.classifyConstraints(obj.l_, obj.u_, obj.C_);
                [obj.rho_vec, obj.rho_inv_vec] = osqp.Solver.makeRhoVec( ...
                    obj.constr_type, obj.opts.rho, obj.m_, obj.C_);
                refactor = true;
            end

            % Px
            if isfield(newData, 'Px')
                Px_new = double(full(newData.Px(:)));
                [ri, ci, ~] = find(obj.P_triu);
                nz_orig = nnz(obj.P_triu);
                if isfield(newData, 'Px_idx')
                    idx = double(full(newData.Px_idx(:)));
                    if numel(Px_new) ~= numel(idx)
                        error('OSQP:update', 'Px and Px_idx must have same length');
                    end
                    old_vals = nonzeros(obj.P_triu);
                    old_vals(idx) = Px_new;
                    obj.P_triu = sparse(ri, ci, old_vals, obj.n_, obj.n_);
                else
                    if numel(Px_new) ~= nz_orig
                        error('OSQP:update', 'Px must have %d elements', nz_orig);
                    end
                    obj.P_triu = sparse(ri, ci, Px_new, obj.n_, obj.n_);
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
                    [ri, ci, vals] = find(obj.A_);
                    vals(idx) = Ax_new;
                    obj.A_ = sparse(ri, ci, vals, obj.m_, obj.n_);
                else
                    nz_orig = nnz(obj.A_);
                    if numel(Ax_new) ~= nz_orig
                        error('OSQP:update', 'Ax must have %d elements', nz_orig);
                    end
                    [ri, ci] = find(obj.A_);
                    obj.A_ = sparse(ri, ci, Ax_new, obj.m_, obj.n_);
                end
                refactor = true;
            end

            % Rescale and refactorize if needed
            if refactor
                [obj.scl, obj.Ps, obj.qs, obj.As, obj.ls, obj.us] = ...
                    osqp.Solver.scaleProblem(obj.P_triu, obj.q_, obj.A_, ...
                    obj.l_, obj.u_, obj.opts, obj.C_);
                obj.factorizeKKT();
            else
                % Only q changed; recompute scaled q
                obj.qs = obj.scl.c * (obj.scl.D .* obj.q_);
            end

            obj.update_time = obj.update_time + toc(t_update);
        end

        function warm_start(obj, varargin)
            % WARM_START  Set warm starting variables x and/or y.
            %
            %   warm_start('x', x, 'y', y)
            %   warm_start('x', x)
            %   warm_start('y', y)
            [n, m] = obj.get_dimensions();

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

            if ~isempty(x)
                obj.x_ = x;
                obj.z_ = obj.A_ * x;
            end
            if ~isempty(y)
                obj.y_ = y;
            end
        end

        function cold_start(obj)
            % COLD_START  Reset solver iterates to zero.
            obj.x_ = zeros(obj.n_, 1);
            obj.z_ = zeros(obj.m_, 1);
            obj.y_ = zeros(obj.m_, 1);
        end

        function update_settings(obj, varargin)
            % UPDATE_SETTINGS  Update solver settings after setup.
            %
            %   update_settings('max_iter', 1000, 'verbose', false)
            %   update_settings(settings_struct)
            if ~obj.isSetup
                error('OSQP:update_settings', 'Call setup() before update_settings().');
            end

            rho_updated = false;
            if isscalar(varargin) && isstruct(varargin{1})
                s = varargin{1};
                fnames = fieldnames(s);
                for k = 1:numel(fnames)
                    name = fnames{k};
                    if ~obj.opts.isUpdatable(name)
                        error('OSQP:update_settings', ...
                            'Setting ''%s'' cannot be updated after setup.', name);
                    end
                    obj.opts.(name) = s.(name);
                    if strcmp(name, 'rho'), rho_updated = true; end
                end
            else
                if mod(numel(varargin), 2) ~= 0
                    error('OSQP:update_settings', 'Arguments must be name/value pairs.');
                end
                for k = 1:2:numel(varargin)
                    name = varargin{k};
                    if ~obj.opts.isUpdatable(name)
                        error('OSQP:update_settings', ...
                            'Setting ''%s'' cannot be updated after setup.', name);
                    end
                    obj.opts.(name) = varargin{k+1};
                    if strcmp(name, 'rho'), rho_updated = true; end
                end
            end

            if rho_updated
                [obj.rho_vec, obj.rho_inv_vec] = osqp.Solver.makeRhoVec( ...
                    obj.constr_type, obj.opts.rho, obj.m_, obj.C_);
                obj.factorizeKKT();
            end
        end

        function out = current_settings(obj)
            % CURRENT_SETTINGS  Get the current solver settings as a struct.
            out = obj.opts.toStruct();
        end

        function [n, m] = get_dimensions(obj)
            % GET_DIMENSIONS  Get the number of variables and constraints.
            n = obj.n_;
            m = obj.m_;
        end
    end

    % =====================================================================
    %  Private helpers
    % =====================================================================
    methods (Access = private)
        function cacheConstants(obj)
            % CACHECONSTANTS  Read base constants from osqp.constant and
            % cache them (with derived values) for the solve loop.
            % Called once at the beginning of setup().  Derived values
            % (INFTY_THRESH, DIVISION_TOL) are computed from base
            % constants rather than stored as literals.
            C.OSQP_INFTY        = osqp.constant.OSQP_INFTY;
            C.OSQP_INFTY_THRESH = osqp.constant.OSQP_INFTY * osqp.constant.OSQP_MIN_SCALING;
            C.OSQP_DIVISION_TOL = 1 / osqp.constant.OSQP_INFTY;
            C.OSQP_RHO_MIN      = osqp.constant.OSQP_RHO_MIN;
            C.OSQP_RHO_MAX      = osqp.constant.OSQP_RHO_MAX;
            C.OSQP_RHO_EQ_OVER_RHO_INEQ = osqp.constant.OSQP_RHO_EQ_OVER_RHO_INEQ;
            C.OSQP_RHO_TOL      = osqp.constant.OSQP_RHO_TOL;
            C.OSQP_MIN_SCALING   = osqp.constant.OSQP_MIN_SCALING;
            C.OSQP_MAX_SCALING   = osqp.constant.OSQP_MAX_SCALING;
            C.OSQP_ADAPTIVE_RHO_MULTIPLE_TERMINATION = ...
                osqp.constant.OSQP_ADAPTIVE_RHO_MULTIPLE_TERMINATION;

            % Status codes
            C.OSQP_SOLVED                       = osqp.constant.OSQP_SOLVED;
            C.OSQP_SOLVED_INACCURATE            = osqp.constant.OSQP_SOLVED_INACCURATE;
            C.OSQP_PRIMAL_INFEASIBLE            = osqp.constant.OSQP_PRIMAL_INFEASIBLE;
            C.OSQP_PRIMAL_INFEASIBLE_INACCURATE = osqp.constant.OSQP_PRIMAL_INFEASIBLE_INACCURATE;
            C.OSQP_DUAL_INFEASIBLE              = osqp.constant.OSQP_DUAL_INFEASIBLE;
            C.OSQP_DUAL_INFEASIBLE_INACCURATE   = osqp.constant.OSQP_DUAL_INFEASIBLE_INACCURATE;
            C.OSQP_MAX_ITER_REACHED             = osqp.constant.OSQP_MAX_ITER_REACHED;
            C.OSQP_TIME_LIMIT_REACHED           = osqp.constant.OSQP_TIME_LIMIT_REACHED;
            C.OSQP_NON_CONVEX                   = osqp.constant.OSQP_NON_CONVEX;
            C.OSQP_UNSOLVED                     = osqp.constant.OSQP_UNSOLVED;

            obj.C_ = C;
        end

        function factorizeKKT(obj)
            % FACTORIZEKKT  Build and factorize the KKT matrix.
            K = osqp.LinearSolver.buildKKT( ...
                obj.Ps, obj.As, obj.rho_vec, obj.opts.sigma, obj.n_, obj.m_);

            switch obj.opts.linear_solver
                case 'matlab_ldl'
                    if isempty(obj.kkt_factor) || ...
                            ~isa(obj.kkt_factor, 'osqp.linsys.MatlabLDLSolver')
                        obj.kkt_factor = osqp.linsys.MatlabLDLSolver(K);
                    else
                        obj.kkt_factor.refactorize(K);
                    end
                case 'qdldl'
                    if isempty(obj.kkt_factor) || ...
                            ~isa(obj.kkt_factor, 'osqp.linsys.QDLDLSolver')
                        obj.kkt_factor = osqp.linsys.QDLDLSolver(K);
                    else
                        obj.kkt_factor.refactorize(K);
                    end
                case 'qdldl_c'
                    if isempty(obj.kkt_factor) || ...
                            ~isa(obj.kkt_factor, 'osqp.linsys.QDLDLCSolver')
                        obj.kkt_factor = osqp.linsys.QDLDLCSolver(K);
                    else
                        obj.kkt_factor.refactorize(K);
                    end
                otherwise
                    error('OSQP:setup', 'Unknown linear solver: %s', ...
                        obj.opts.linear_solver);
            end
        end
    end

    % =====================================================================
    %  Static helpers
    % =====================================================================
    methods (Static, Access = private)

        function str = statusStr(val)
            switch val
                case  1,  str = 'solved';
                case  2,  str = 'solved_inaccurate';
                case  3,  str = 'primal_infeasible';
                case  4,  str = 'primal_infeasible_inaccurate';
                case  5,  str = 'dual_infeasible';
                case  6,  str = 'dual_infeasible_inaccurate';
                case  7,  str = 'maximum_iterations_reached';
                case  8,  str = 'time_limit_reached';
                case  9,  str = 'non_convex';
                case 10,  str = 'interrupted';
                case 11,  str = 'unsolved';
                otherwise, str = 'unsolved';
            end
        end

        function tf = isLDLpsd(P)
            % ISLDLPSD  Check if symmetric P is positive semidefinite via LDL'.
            %   MATLAB's ldl() produces a block-diagonal D with 1x1 and 2x2
            %   blocks.  P is PSD iff every eigenvalue of every block is
            %   nonneg.  Unlike Cholesky, LDL handles exact zeros in D.
            %   Matches C behavior: QDLDL counts strictly positive D entries;
            %   zero entries (from PSD P+sigma*I) still yield n positive
            %   entries because sigma shifts them above zero.  Here we check
            %   P alone (without sigma) to detect indefiniteness early.
            [~, D, ~] = ldl(sparse(P));
            n = size(D, 1);
            if n == 0, tf = true; return; end
            d = diag(D);
            if n > 1
                sub = diag(D, -1);
            else
                sub = zeros(0, 1);
            end
            tf = true;
            k = 1;
            while k <= n
                if k < n && sub(k) ~= 0
                    % 2x2 block: PSD iff trace >= 0 and det >= 0
                    tr = d(k) + d(k+1);
                    dt = d(k) * d(k+1) - sub(k)^2;
                    if tr < 0 || dt < 0
                        tf = false; return;
                    end
                    k = k + 2;
                else
                    if d(k) < 0
                        tf = false; return;
                    end
                    k = k + 1;
                end
            end
        end

        function r = emptyResults(n, m, C)
            r.x = nan(n, 1);
            r.y = nan(m, 1);
            r.prim_inf_cert = nan(m, 1);
            r.dual_inf_cert = nan(n, 1);
            r.info.iter          = 0;
            r.info.status        = 'unsolved';
            r.info.status_val    = C.OSQP_UNSOLVED;
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

        function ctype = classifyConstraints(l, u, C)
            % CLASSIFYCONSTRAINTS  Constraint type for per-row rho.
            %   Matches C code ew_bounds_type() in auxil.c:
            %   -1 = loose (both bounds >= INFTY_THRESH) → uses RHO_MIN
            %    0 = inequality (l != u)                 → uses rho
            %    1 = equality  (|l - u| < RHO_TOL)       → uses RHO_EQ * rho
            m = numel(l);
            thr = C.OSQP_INFTY_THRESH;
            ctype = zeros(m, 1);             % default: inequality (0)
            loose = (l <= -thr) & (u >= thr);
            ctype(loose) = -1;
            eq = abs(l - u) < C.OSQP_RHO_TOL;
            ctype(eq) = 1;
        end

        function [rho_vec, rho_inv_vec] = makeRhoVec(ctype, rho, m, C)
            % MAKERHOVEX  Per-constraint rho vector (C: set_rho_vec).
            %   Assigns different rho values by constraint type:
            %   ctype -1 (loose):     RHO_MIN  (1e-6)
            %   ctype  0 (inequality): rho
            %   ctype  1 (equality):   RHO_EQ_OVER_RHO_INEQ * rho  (1e3 * rho)
            if m == 0
                rho_vec = zeros(0, 1); rho_inv_vec = zeros(0, 1); return;
            end
            rho = max(min(rho, C.OSQP_RHO_MAX), C.OSQP_RHO_MIN);
            rho_vec = rho * ones(m, 1);
            rho_vec(ctype == -1) = C.OSQP_RHO_MIN;
            rho_vec(ctype == 1)  = C.OSQP_RHO_EQ_OVER_RHO_INEQ * rho;
            rho_inv_vec = 1 ./ rho_vec;
        end

        function [scl, Ps, qs, As, ls, us] = scaleProblem(P_triu, q, A, l, u, opts, C)
            % SCALEPROBLEM  Ruiz equilibration (C: scaling.c / scale_data).
            %   Each iteration: compute column inf-norms of the KKT matrix
            %   [P A'; A 0], take their sqrt-inverse as scaling factors D/E,
            %   then apply to P, A, q.  Finally normalize cost by c.
            %   Accumulates D, E, c over opts.scaling iterations.
            n = size(P_triu, 1);
            m = size(A, 1);
            num_iter = opts.scaling;
            INF = C.OSQP_INFTY;
            MIN_S = C.OSQP_MIN_SCALING;
            MAX_S = C.OSQP_MAX_SCALING;

            D = ones(n, 1);
            E = ones(m, 1);
            c = 1.0;

            if num_iter == 0
                Ps = P_triu; qs = q; As = A; ls = l; us = u;
                scl.D = D; scl.E = E; scl.Dinv = D; scl.Einv = E;
                scl.c = 1; scl.cinv = 1;
                return;
            end

            % Work on full symmetric P for norm computation; modify in-place
            Pw = P_triu + P_triu' - diag(diag(P_triu));
            Aw = A;
            qw = q;

            for k = 1:num_iter
                % --- Scaling iteration k ---
                % Compute column inf-norms of [P; A] (the KKT top block).
                % D_temp scales variables (columns), E_temp scales constraints (rows).
                D_temp = full(max(abs(Pw)))';           % col inf-norms of P (n×1)
                if m > 0
                    D_temp = max(D_temp, full(max(abs(Aw)))');
                end

                % E_temp = row inf-norms of A_current
                if m > 0
                    E_temp = full(max(abs(Aw), [], 2)); % row inf-norms (m×1)
                else
                    E_temp = ones(0, 1);
                end

                % Clamp scaling factors to [MIN_SCALING, MAX_SCALING]
                % (C: limit_scaling_vector)
                D_temp(D_temp < MIN_S) = 1.0;
                D_temp = min(D_temp, MAX_S);
                if m > 0
                    E_temp(E_temp < MIN_S) = 1.0;
                    E_temp = min(E_temp, MAX_S);
                end

                % Sqrt and invert
                D_temp = 1 ./ sqrt(D_temp);
                if m > 0
                    E_temp = 1 ./ sqrt(E_temp);
                end

                % Scale P, A, q in-place (matching C)
                Dsp = spdiags(D_temp, 0, n, n);
                Pw = Dsp * Pw * Dsp;
                if m > 0
                    Esp = spdiags(E_temp, 0, m, m);
                    Aw = Esp * Aw * Dsp;
                end
                qw = D_temp .* qw;

                % Accumulate
                D = D .* D_temp;
                if m > 0
                    E = E .* E_temp;
                end

                % --- Cost normalization ---
                % Normalize objective by the average column inf-norm of P
                % combined with the inf-norm of q (C: compute_inf_norm_cols_P_scaled)
                col_norms_P = full(max(abs(Pw)));
                avg_col_norm = sum(col_norms_P, 2) / max(n, 1);
                inf_norm_q = norm(qw, inf);

                cost_measure = max(avg_col_norm, inf_norm_q);
                % limit_scaling_scalar
                if cost_measure < MIN_S
                    cost_measure = 1.0;
                end
                cost_measure = min(cost_measure, MAX_S);

                c_temp = 1 / cost_measure;
                Pw = c_temp * Pw;
                qw = c_temp * qw;
                c = c * c_temp;
            end

            Dinv = 1 ./ D;
            Einv = 1 ./ E;
            cinv = 1 / c;

            Ps = triu(Pw);
            qs = qw;
            As = Aw;
            ls = max(E .* l, -INF);
            us = min(E .* u,  INF);

            scl.D = D; scl.E = E; scl.Dinv = Dinv; scl.Einv = Einv;
            scl.c = c; scl.cinv = cinv;
        end

        function [converged, status_val] = checkConvergence( ...
                xs, zs, ys, xs_prev, ys_prev, ...
                Ps, qs, As, ls, us, ~, m, s, scl, C, approximate)
            % CHECKCONVERGENCE  Termination checks matching C code.
            %
            %   When s.scaled_termination is true (or scaling is off),
            %   all comparisons use scaled quantities directly.
            %
            %   When s.scaled_termination is false (C default) and
            %   scaling is active, residuals and tolerances are unscaled
            %   using D, Dinv, E, Einv, c, cinv from scl.
            %
            %   When approximate=true, uses 10x tolerances and returns
            %   INACCURATE status variants.
            if nargin < 16, approximate = false; end
            status_val = C.OSQP_UNSOLVED;
            converged  = false;

            % Determine whether to unscale
            do_unscale = s.scaling > 0 && ~s.scaled_termination;

            if approximate
                tol_mult = 10;
            else
                tol_mult = 1;
            end

            Pfull = Ps + Ps' - diag(diag(Ps));
            Pxs = Pfull * xs;

            % --- Primal residual: ||Ax - z||_inf ---
            %   (optionally unscaled: ||Einv .* (Ax - z)||_inf)
            if m > 0
                prim_resid_vec = As * xs - zs;
                Axs = As * xs;
            else
                prim_resid_vec = zeros(0, 1);
                Axs = zeros(0, 1);
            end

            if m > 0
                if do_unscale
                    prim_res = norm(scl.Einv .* prim_resid_vec, inf);
                else
                    prim_res = norm(prim_resid_vec, inf);
                end
            else
                prim_res = 0;
            end

            % --- Dual residual: ||Px + q + A'y||_inf ---
            %   (optionally unscaled: cinv * ||Dinv .* (Px + q + A'y)||_inf)
            if m > 0
                Atys = As' * ys;
                dual_resid_vec = Pxs + qs + Atys;
            else
                Atys = zeros(size(qs));
                dual_resid_vec = Pxs + qs;
            end
            if do_unscale
                dual_res = scl.cinv * norm(scl.Dinv .* dual_resid_vec, inf);
            else
                dual_res = norm(dual_resid_vec, inf);
            end

            % --- Check for divergence (non-convexity) ---
            if prim_res > C.OSQP_INFTY || ...
                    dual_res > C.OSQP_INFTY
                status_val = C.OSQP_NON_CONVEX;
                converged = true;
                return;
            end

            % --- Primal tolerance: eps_abs + eps_rel * max(||Ax||, ||z||) ---
            if m > 0
                if do_unscale
                    max_rel_prim = max( ...
                        norm(scl.Einv .* zs, inf), ...
                        norm(scl.Einv .* Axs, inf));
                else
                    max_rel_prim = max(norm(zs, inf), norm(Axs, inf));
                end
                eps_prim = tol_mult * (s.eps_abs + s.eps_rel * max_rel_prim);
            else
                eps_prim = tol_mult * s.eps_abs;
            end

            % --- Dual tolerance: eps_abs + eps_rel * max(||Px||, ||q||, ||A'y||) ---
            if do_unscale
                max_rel_dual = norm(scl.Dinv .* qs, inf);
                if m > 0
                    max_rel_dual = max(max_rel_dual, ...
                        norm(scl.Dinv .* Atys, inf));
                end
                max_rel_dual = max(max_rel_dual, ...
                    norm(scl.Dinv .* Pxs, inf));
                max_rel_dual = max_rel_dual * scl.cinv;
            else
                max_rel_dual = norm(qs, inf);
                if m > 0
                    max_rel_dual = max(max_rel_dual, norm(Atys, inf));
                end
                max_rel_dual = max(max_rel_dual, norm(Pxs, inf));
            end
            eps_dual = tol_mult * (s.eps_abs + s.eps_rel * max_rel_dual);

            % --- Optimality ---
            if prim_res <= eps_prim && dual_res <= eps_dual
                if approximate
                    status_val = C.OSQP_SOLVED_INACCURATE;
                else
                    status_val = C.OSQP_SOLVED;
                end
                converged = true;
                return;
            end

            % --- Primal infeasibility certificate (C: is_primal_infeasible) ---
            %   If ||A' * delta_y|| < eps * ||delta_y|| and
            %   support_fn(delta_y, l, u) < 0, problem is primal infeasible.
            if m > 0
                dy  = ys - ys_prev;
                eps_pi = tol_mult * s.eps_prim_inf;

                % Project delta_y onto polar of recession cone of [l, u]
                INF_THRESH = C.OSQP_INFTY_THRESH;
                for i = 1:m
                    li_fin = abs(ls(i)) < INF_THRESH;
                    ui_fin = abs(us(i)) < INF_THRESH;
                    if ~li_fin && ~ui_fin
                        dy(i) = 0;
                    elseif ~li_fin
                        dy(i) = min(dy(i), 0);
                    elseif ~ui_fin
                        dy(i) = max(dy(i), 0);
                    end
                end

                % Compute norm (unscale if needed)
                if do_unscale
                    norm_dy = norm(scl.E .* dy, inf);
                else
                    norm_dy = norm(dy, inf);
                end

                if norm_dy > eps_pi
                    % Support function check: u'*max(dy,0) + l'*min(dy,0)
                    sup_val = 0;
                    for i = 1:m
                        if abs(us(i)) < INF_THRESH
                            sup_val = sup_val + us(i) * max(dy(i), 0);
                        end
                        if abs(ls(i)) < INF_THRESH
                            sup_val = sup_val + ls(i) * min(dy(i), 0);
                        end
                    end

                    if sup_val < 0
                        ATdy = As' * dy;
                        if do_unscale
                            ATdy = scl.Dinv .* ATdy;
                        end
                        if norm(ATdy, inf) < eps_pi * norm_dy
                            if approximate
                                status_val = C.OSQP_PRIMAL_INFEASIBLE_INACCURATE;
                            else
                                status_val = C.OSQP_PRIMAL_INFEASIBLE;
                            end
                            converged = true;
                            return;
                        end
                    end
                end
            end

            % --- Dual infeasibility certificate (C: is_dual_infeasible) ---
            %   If q' * delta_x < 0 and ||P * delta_x|| < eps * ||delta_x||
            %   and A * delta_x in recession cone of [l, u], problem is
            %   dual infeasible (unbounded below).
            dx = xs - xs_prev;
            eps_di = tol_mult * s.eps_dual_inf;

            % Compute norm and cost_scaling (unscale if needed)
            if do_unscale
                norm_dx = norm(scl.D .* dx, inf);
                cost_scaling = scl.c;
            else
                norm_dx = norm(dx, inf);
                cost_scaling = 1.0;
            end

            if norm_dx > C.OSQP_DIVISION_TOL
                % Check q' * delta_x < 0
                qdx = qs' * dx;
                if qdx < 0
                    % Check ||P * delta_x|| < cost_scaling * eps * ||delta_x||
                    Pdx = Pfull * dx;
                    if do_unscale
                        Pdx = scl.Dinv .* Pdx;
                    end
                    if norm(Pdx, inf) < cost_scaling * eps_di * norm_dx
                        % Check A * delta_x in recession cone of [l, u]
                        if m > 0
                            Adx = As * dx;
                            if do_unscale
                                Adx = scl.Einv .* Adx;
                            end
                            ok = osqp.Solver.dualInfCheck( ...
                                Adx, ls, us, eps_di * norm_dx, C);
                        else
                            ok = true;
                        end
                        if ok
                            if approximate
                                status_val = C.OSQP_DUAL_INFEASIBLE_INACCURATE;
                            else
                                status_val = C.OSQP_DUAL_INFEASIBLE;
                            end
                            converged = true;
                            return;
                        end
                    end
                end
            end
        end

        function ok = primalInfCheck(v, l, u, eps, C)
            % Support function of polar recession cone.
            % Used only as a helper — the main primal inf check
            % is now inline in checkConvergence with proper unscaling.
            INF_THRESH = C.OSQP_INFTY_THRESH;
            vpos = max(v, 0);
            vneg = min(v, 0);
            val  = 0;
            for i = 1:numel(l)
                if abs(u(i)) < INF_THRESH, val = val + u(i) * vpos(i); end
                if abs(l(i)) < INF_THRESH, val = val + l(i) * vneg(i); end
            end
            ok = val < -eps;
        end

        function ok = dualInfCheck(Adx, l, u, eps, C)
            % Check whether Adx lies in the recession cone of [l, u].
            % eps is pre-multiplied by norm_delta_x by the caller,
            % matching C code's in_reccone(Adx, l, u, thr, eps*norm).
            INF_THRESH = C.OSQP_INFTY_THRESH;
            ok = true;
            for i = 1:numel(l)
                li = l(i); ui = u(i);
                li_fin = abs(li) < INF_THRESH;
                ui_fin = abs(ui) < INF_THRESH;
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

        function new_rho = computeNewRho(xs, zs, ys, Ps, qs, As, rho, tol, ~, m, C)
            % COMPUTENEWRHO  Adaptive rho estimate (C: compute_rho_estimate).
            %   rho_new = rho * sqrt(prim_res_normalized / dual_res_normalized)
            %   Only updates if the ratio exceeds adaptive_rho_tolerance.
            DTOL = C.OSQP_DIVISION_TOL;
            Pfull = Ps + Ps' - diag(diag(Ps));

            % Primal residual
            if m > 0
                prim_res = norm(As * xs - zs, inf);
            else
                prim_res = 0;
            end

            % Dual residual
            Pxs = Pfull * xs;
            if m > 0
                dual_res = norm(Pxs + qs + As' * ys, inf);
            else
                dual_res = norm(Pxs + qs, inf);
            end

            if prim_res < DTOL || dual_res < DTOL
                new_rho = rho; return;
            end

            % Normalize primal residual: max(||z||, ||Ax||)
            prim_norm = norm(zs, inf);
            if m > 0
                prim_norm = max(prim_norm, norm(As * xs, inf));
            end
            prim_res_n = prim_res / (prim_norm + DTOL);

            % Normalize dual residual: max(||q||, ||A'y||, ||Px||)
            dual_norm = norm(qs, inf);
            if m > 0
                dual_norm = max(dual_norm, norm(As' * ys, inf));
            end
            dual_norm = max(dual_norm, norm(Pxs, inf));
            dual_res_n = dual_res / (dual_norm + DTOL);

            rho_estimate = rho * sqrt(prim_res_n / dual_res_n);
            rho_estimate = max(min(rho_estimate, C.OSQP_RHO_MAX), ...
                               C.OSQP_RHO_MIN);

            if rho_estimate > rho * tol || rho_estimate < rho / tol
                new_rho = rho_estimate;
            else
                new_rho = rho;
            end
        end

        function results = polishSolution(results, P_triu, q, A, l, u, s, C)
            % POLISHSOLUTION  Solution polishing (C: polish.c).
            %   1. Identify active constraints from dual sign and primal slack.
            %   2. If no active set: solve unconstrained with delta-regularized P.
            %   3. Otherwise: form reduced KKT with active rows of A,
            %      solve with iterative refinement, then apply Moreau
            %      decomposition z = clip(y + Ax, l, u) for complementarity.
            n = numel(q);
            m = size(A, 1);
            x = results.x;
            y = results.y;
            Pfull = P_triu + P_triu' - diag(diag(P_triu));
            INF = C.OSQP_INFTY;

            % Determine active set from sign of dual and primal distance
            Ax = A * x;
            tol_act = max(s.delta, s.eps_abs);
            act_l = (y < -s.delta & abs(l) < INF) | ...
                (abs(l) < INF & abs(Ax - l) <= tol_act);
            act_u = (y >  s.delta & abs(u) < INF) | ...
                (abs(u) < INF & abs(Ax - u) <= tol_act);
            active_idx = find(act_l | act_u);
            n_act = numel(active_idx);

            if n_act == 0
                % No active constraints → solve unconstrained regularized QP
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

            % Form reduced KKT system with active constraints (C: form_Ared)
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
                % Solve delta-regularized KKT system
                %   [P + delta*I   A_act'] [x]   [-q   ]
                %   [A_act        -delta*I] [y] = [b_act]
                Kreg = [Pfull + s.delta*speye(n), A_act'; ...
                    A_act, -s.delta*speye(n_act)];
                rhs0 = [-q; b_act];
                sol0 = Kreg \ rhs0;
                x_pol = sol0(1:n);
                y_act = sol0(n+1:end);

                % Iterative refinement (C: polish_refine_iter loops)
                for refine = 1:s.polish_refine_iter
                    res = rhs0 - Kreg * [x_pol; y_act];
                    delta_sol = Kreg \ res;
                    x_pol = x_pol + delta_sol(1:n);
                    y_act = y_act + delta_sol(n+1:end);
                end

                y_pol = zeros(m, 1);
                y_pol(active_idx) = y_act;

                Ax_pol = A * x_pol;

                % Moreau decomposition: project (y + Ax) onto [l,u] to
                % get z, then y_pol = (y + Ax) - z ensures z ∈ C, y ∈ N_C(z)
                y_pol = y_pol + Ax_pol;
                z_pol = min(max(y_pol, l), u);
                y_pol = y_pol - z_pol;

                feas_l = all(z_pol >= l - tol_act | abs(l) >= INF);
                feas_u = all(z_pol <= u + tol_act | abs(u) >= INF);
                if feas_l && feas_u
                    obj_pol = 0.5 * x_pol' * Pfull * x_pol + q' * x_pol;
                    pri_res_pol = norm(Ax_pol - z_pol, inf);
                    dua_res_pol = norm(Pfull * x_pol + q + A' * y_pol, inf);
                    pri_res_cur = results.info.pri_res;
                    dua_res_cur = results.info.dua_res;
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
