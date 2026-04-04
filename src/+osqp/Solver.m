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
    end

    % -----------------------------------------------------------------
    %  Status constants
    % -----------------------------------------------------------------
    properties (Constant, Hidden = true)
        OSQP_INFTY = 1e30

        % Rho management constants (matching C code)
        RHO_MIN              = 1e-6
        RHO_MAX              = 1e6
        RHO_EQ_OVER_RHO_INEQ = 1e3
        RHO_TOL              = 1e-4
        MIN_SCALING          = 1e-4
        MAX_SCALING          = 1e4
        DIVISION_TOL         = 1e-10

        STATUS_SOLVED                       = 1
        STATUS_SOLVED_INACCURATE            = 2
        STATUS_PRIMAL_INFEASIBLE            = 3
        STATUS_PRIMAL_INFEASIBLE_INACCURATE = 4
        STATUS_DUAL_INFEASIBLE              = 5
        STATUS_DUAL_INFEASIBLE_INACCURATE   = 6
        STATUS_MAX_ITER_REACHED             = 7
        STATUS_TIME_LIMIT_REACHED           = 8
        STATUS_NON_CONVEX                   = 9
        STATUS_SIGINT                       = 10
        STATUS_UNSOLVED                     = 11
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
            obj.l_ = max(l, -osqp.Solver.OSQP_INFTY);
            obj.u_ = min(u,  osqp.Solver.OSQP_INFTY);

            % --- Convexity check ---
            Pfull = obj.P_triu + obj.P_triu' - diag(diag(obj.P_triu));
            obj.non_convex = false;
            if n > 0 && nnz(Pfull) > 0
                thresh = 1e-7;
                [~, flag] = chol(Pfull + thresh * speye(n), 'lower');
                if flag ~= 0
                    Psig = Pfull + obj.opts.sigma * speye(n);
                    [~, flag2] = chol(Psig + thresh * speye(n), 'lower');
                    if flag2 ~= 0
                        error('OSQP:NonConvex', ...
                            'P is non-convex and sigma is too small.');
                    end
                    obj.non_convex = true;
                end
            end

            % --- Initialize iterates ---
            obj.x_ = zeros(n, 1);
            obj.z_ = zeros(m, 1);
            obj.y_ = zeros(m, 1);

            % --- Classify constraints ---
            obj.constr_type = osqp.Solver.classifyConstraints(obj.l_, obj.u_);

            % --- Scale problem ---
            [obj.scl, obj.Ps, obj.qs, obj.As, obj.ls, obj.us] = ...
                osqp.Solver.scaleProblem(obj.P_triu, obj.q_, obj.A_, ...
                obj.l_, obj.u_, obj.opts);

            % --- Build rho vector ---
            [obj.rho_vec, obj.rho_inv_vec] = osqp.Solver.makeRhoVec( ...
                obj.constr_type, obj.opts.rho, m);

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

            results = osqp.Solver.emptyResults(n, m);

            % Non-convex check
            if obj.non_convex
                results.info.status_val = osqp.Solver.STATUS_NON_CONVEX;
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

            status_val = osqp.Solver.STATUS_UNSOLVED;
            rho_updates = 0;
            adaptive_rho_interval = 0;

            if s.adaptive_rho && s.adaptive_rho_interval > 0
                adaptive_rho_interval = s.adaptive_rho_interval;
            end

            iter = 0; %#ok<NASGU> default for max_iter=0

            for iter = 1:s.max_iter
                xs_prev = xs;
                zs_prev = zs;
                ys_prev = ys;

                % Step 1: Solve KKT for x_tilde, z_tilde
                rhs_x = s.sigma * xs_prev - obj.qs;
                if m > 0
                    rhs_z = zs_prev - obj.rho_inv_vec .* ys;
                    rhs   = [rhs_x; rhs_z];
                else
                    rhs = rhs_x;
                end

                sol = obj.kkt_factor \ rhs;
                xtilde = sol(1:n);

                if m > 0
                    ztilde = zs_prev + obj.rho_inv_vec .* (sol(n+1:end) - ys);
                else
                    ztilde = zeros(0, 1);
                end

                % Step 2: Relaxation + projection
                xs = s.alpha * xtilde + (1 - s.alpha) * xs_prev;
                if m > 0
                    zs_relaxed = s.alpha * ztilde + (1 - s.alpha) * zs_prev;
                    zs = min(max(zs_relaxed + obj.rho_inv_vec .* ys, ...
                        obj.ls), obj.us);
                end

                % Step 3: Dual update
                if m > 0
                    ys = ys + obj.rho_vec .* ( ...
                        s.alpha * ztilde + (1 - s.alpha) * zs_prev - zs);
                end

                % Convergence check
                if s.check_termination > 0 && mod(iter, s.check_termination) == 0
                    [converged, status_val] = osqp.Solver.checkConvergence( ...
                        xs, zs, ys, xs_prev, ys_prev, ...
                        obj.Ps, obj.qs, obj.As, obj.ls, obj.us, ...
                        n, m, s);
                    if converged
                        break;
                    end
                end

                % Time limit
                if s.time_limit > 0 && toc(t_solve) >= s.time_limit
                    % Try approximate check before declaring time limit
                    [conv_tl, sv_tl] = osqp.Solver.checkConvergence( ...
                        xs, zs, ys, xs_prev, ys_prev, ...
                        obj.Ps, obj.qs, obj.As, obj.ls, obj.us, ...
                        n, m, s, true);  % approximate = true
                    if conv_tl
                        status_val = sv_tl;
                    else
                        status_val = osqp.Solver.STATUS_TIME_LIMIT_REACHED;
                    end
                    break;
                end

                % Adaptive rho
                if s.adaptive_rho
                    if adaptive_rho_interval == 0
                        if iter == 1
                            t_iter_start = tic;
                        elseif iter == 2
                            t_one_iter = toc(t_iter_start);
                            t_setup_val = max(obj.setup_time, 1e-10);
                            if t_one_iter > 0
                                adaptive_rho_interval = max(1, ...
                                    round(s.adaptive_rho_fraction * t_setup_val / t_one_iter));
                            else
                                adaptive_rho_interval = 25;
                            end
                        end
                    end
                    if adaptive_rho_interval > 0 && mod(iter, adaptive_rho_interval) == 0
                        new_rho = osqp.Solver.computeNewRho( ...
                            xs, zs, ys, obj.Ps, obj.qs, obj.As, ...
                            obj.opts.rho, s.adaptive_rho_tolerance, n, m);
                        if new_rho ~= obj.opts.rho
                            obj.opts.rho = new_rho;
                            [obj.rho_vec, obj.rho_inv_vec] = osqp.Solver.makeRhoVec( ...
                                obj.constr_type, new_rho, m);
                            obj.factorizeKKT();
                            rho_updates = rho_updates + 1;
                        end
                    end
                end
            end % for iter

            obj.solve_time = toc(t_solve);

            % Final convergence test if loop ended without converging
            if status_val == osqp.Solver.STATUS_UNSOLVED
                [~, status_val] = osqp.Solver.checkConvergence( ...
                    xs, zs, ys, xs_prev, ys_prev, ...
                    obj.Ps, obj.qs, obj.As, obj.ls, obj.us, ...
                    n, m, s);
                % If still unsolved, try approximate check (10x tolerances)
                % matching C code's post-loop check_termination(solver, 1)
                if status_val == osqp.Solver.STATUS_UNSOLVED
                    [~, status_val] = osqp.Solver.checkConvergence( ...
                        xs, zs, ys, xs_prev, ys_prev, ...
                        obj.Ps, obj.qs, obj.As, obj.ls, obj.us, ...
                        n, m, s, true);  % approximate = true
                end
                if status_val == osqp.Solver.STATUS_UNSOLVED
                    status_val = osqp.Solver.STATUS_MAX_ITER_REACHED;
                end
            end

            % Unscale solution
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

            solution_present = ismember(status_val, [ ...
                osqp.Solver.STATUS_SOLVED, ...
                osqp.Solver.STATUS_SOLVED_INACCURATE, ...
                osqp.Solver.STATUS_MAX_ITER_REACHED]);

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
            elseif status_val == osqp.Solver.STATUS_PRIMAL_INFEASIBLE || ...
                    status_val == osqp.Solver.STATUS_PRIMAL_INFEASIBLE_INACCURATE
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
            elseif status_val == osqp.Solver.STATUS_DUAL_INFEASIBLE || ...
                    status_val == osqp.Solver.STATUS_DUAL_INFEASIBLE_INACCURATE
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

            % Polishing (C code only polishes when status == OSQP_SOLVED)
            results.info.status_polish = 0;
            if s.polishing && status_val == osqp.Solver.STATUS_SOLVED
                t_polish = tic;
                results = osqp.Solver.polishSolution(results, ...
                    obj.P_triu, obj.q_, obj.A_, obj.l_, obj.u_, s);
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
                    obj.l_ = max(l_new, -osqp.Solver.OSQP_INFTY);
                end
                if isfield(newData, 'u')
                    u_new = double(full(newData.u(:)));
                    if numel(u_new) ~= obj.m_
                        error('OSQP:update', 'u must have length m=%d', obj.m_);
                    end
                    obj.u_ = min(u_new, osqp.Solver.OSQP_INFTY);
                end
                obj.constr_type = osqp.Solver.classifyConstraints(obj.l_, obj.u_);
                [obj.rho_vec, obj.rho_inv_vec] = osqp.Solver.makeRhoVec( ...
                    obj.constr_type, obj.opts.rho, obj.m_);
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
                    obj.l_, obj.u_, obj.opts);
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

            x_updated = false;
            y_updated = false;
            if ~isempty(x)
                obj.x_ = x;
                obj.z_ = obj.A_ * x;
                x_updated = true;
            end
            if ~isempty(y)
                obj.y_ = y;
                y_updated = true;
            end
            if x_updated && ~y_updated
                obj.y_ = zeros(m, 1);
            elseif ~x_updated && y_updated
                obj.x_ = zeros(n, 1);
                obj.z_ = zeros(m, 1);
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
                    obj.constr_type, obj.opts.rho, obj.m_);
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

        function r = emptyResults(n, m)
            r.x = nan(n, 1);
            r.y = nan(m, 1);
            r.prim_inf_cert = nan(m, 1);
            r.dual_inf_cert = nan(n, 1);
            r.info.iter          = 0;
            r.info.status        = 'unsolved';
            r.info.status_val    = osqp.Solver.STATUS_UNSOLVED;
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

        function ctype = classifyConstraints(l, u)
            % Classify constraints matching C code ew_bounds_type:
            %   -1 = loose (both bounds >= INFTY*MIN_SCALING) -> RHO_MIN
            %    0 = inequality (l != u)                      -> rho
            %    1 = equality  (|l - u| < RHO_TOL)            -> RHO_EQ * rho
            m = numel(l);
            INF = osqp.Solver.OSQP_INFTY;
            MIN_S = osqp.Solver.MIN_SCALING;
            thr = INF * MIN_S;
            ctype = zeros(m, 1);             % default: inequality (0)
            loose = (l <= -thr) & (u >= thr);
            ctype(loose) = -1;
            eq = abs(l - u) < osqp.Solver.RHO_TOL;
            ctype(eq) = 1;
        end

        function [rho_vec, rho_inv_vec] = makeRhoVec(ctype, rho, m)
            % Build per-constraint rho vector matching C code set_rho_vec:
            %   ctype -1 (loose):     RHO_MIN
            %   ctype  0 (inequality): rho
            %   ctype  1 (equality):   RHO_EQ_OVER_RHO_INEQ * rho
            if m == 0
                rho_vec = zeros(0, 1); rho_inv_vec = zeros(0, 1); return;
            end
            rho = max(min(rho, osqp.Solver.RHO_MAX), osqp.Solver.RHO_MIN);
            rho_vec = rho * ones(m, 1);
            rho_vec(ctype == -1) = osqp.Solver.RHO_MIN;
            rho_vec(ctype == 1)  = osqp.Solver.RHO_EQ_OVER_RHO_INEQ * rho;
            rho_inv_vec = 1 ./ rho_vec;
        end

        function [scl, Ps, qs, As, ls, us] = scaleProblem(P_triu, q, A, l, u, opts)
            % Ruiz equilibration matching C code scale_data / compute_inf_norm_cols_KKT.
            % Scales P, A, q in-place each iteration (as C does), accumulating D, E, c.
            n = size(P_triu, 1);
            m = size(A, 1);
            num_iter = opts.scaling;
            INF = osqp.Solver.OSQP_INFTY;
            MIN_S = osqp.Solver.MIN_SCALING;
            MAX_S = osqp.Solver.MAX_SCALING;

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
                % --- Column inf-norms of KKT top block ---
                % D_temp = max(col_inf_norm(P_current), col_inf_norm(A_current))
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

                % Limit scaling (matching limit_scaling_vector in C)
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

                % --- Cost normalization (matching C) ---
                % avg column inf-norm of current P
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
                Ps, qs, As, ls, us, ~, m, s, approximate)
            % Check convergence. When approximate=true, uses 10x tolerances
            % and returns INACCURATE status variants (matching C code).
            if nargin < 14, approximate = false; end
            status_val = osqp.Solver.STATUS_UNSOLVED;
            converged  = false;

            if approximate
                tol_mult = 10;
            else
                tol_mult = 1;
            end

            Pfull = Ps + Ps' - diag(diag(Ps));
            Pxs = Pfull * xs;

            % --- Check for divergence (non-convexity) ---
            if m > 0
                prim_res = norm(As * xs - zs, inf);
            else
                prim_res = 0;
            end
            if m > 0
                dual_res = norm(Pxs + qs + As' * ys, inf);
            else
                dual_res = norm(Pxs + qs, inf);
            end
            if prim_res > osqp.Solver.OSQP_INFTY || ...
                    dual_res > osqp.Solver.OSQP_INFTY
                status_val = osqp.Solver.STATUS_NON_CONVEX;
                converged = true;
                return;
            end

            % Adaptive tolerances
            if m > 0
                Axs_norm = norm(As * xs, inf);
                zs_norm  = norm(zs, inf);
                eps_prim = tol_mult * (s.eps_abs + s.eps_rel * max(Axs_norm, zs_norm));
            else
                eps_prim = tol_mult * s.eps_abs;
            end
            Pxs_norm  = norm(Pxs, inf);
            Atys_norm = 0;
            if m > 0
                Atys_norm = norm(As' * ys, inf);
            end
            qs_norm  = norm(qs, inf);
            eps_dual = tol_mult * (s.eps_abs + s.eps_rel * max([Pxs_norm; Atys_norm; qs_norm]));

            % Optimality
            if prim_res <= eps_prim && dual_res <= eps_dual
                if approximate
                    status_val = osqp.Solver.STATUS_SOLVED_INACCURATE;
                else
                    status_val = osqp.Solver.STATUS_SOLVED;
                end
                converged = true;
                return;
            end

            % Primal infeasibility (certificate: delta_y)
            if m > 0
                dy      = ys - ys_prev;
                norm_dy = norm(dy, inf);
                eps_pi  = tol_mult * s.eps_prim_inf;
                if norm_dy > eps_pi
                    dy_n  = dy / norm_dy;
                    ATdy  = norm(As' * dy_n, inf);
                    lu_ok = osqp.Solver.primalInfCheck(dy_n, ls, us, eps_pi);
                    if ATdy <= eps_pi && lu_ok
                        if approximate
                            status_val = osqp.Solver.STATUS_PRIMAL_INFEASIBLE_INACCURATE;
                        else
                            status_val = osqp.Solver.STATUS_PRIMAL_INFEASIBLE;
                        end
                        converged = true;
                        return;
                    end
                end
            end

            % Dual infeasibility (certificate: delta_x)
            dx      = xs - xs_prev;
            norm_dx = norm(dx, inf);
            eps_di  = tol_mult * s.eps_dual_inf;
            if norm_dx > eps_di
                dx_n = dx / norm_dx;
                Pdx  = norm(Pfull * dx_n, inf);
                qdx  = qs' * dx_n;
                if Pdx <= eps_di && qdx < -eps_di
                    if m > 0
                        Adx = As * dx_n;
                        ok  = osqp.Solver.dualInfCheck(Adx, ls, us, eps_di);
                    else
                        ok = true;
                    end
                    if ok
                        if approximate
                            status_val = osqp.Solver.STATUS_DUAL_INFEASIBLE_INACCURATE;
                        else
                            status_val = osqp.Solver.STATUS_DUAL_INFEASIBLE;
                        end
                        converged = true;
                        return;
                    end
                end
            end
        end

        function ok = primalInfCheck(v, l, u, eps)
            INF = osqp.Solver.OSQP_INFTY;
            vpos = max(v, 0);
            vneg = min(v, 0);
            val  = 0;
            for i = 1:numel(l)
                if abs(u(i)) < INF, val = val + u(i) * vpos(i); end
                if abs(l(i)) < INF, val = val + l(i) * vneg(i); end
            end
            ok = val < -eps;
        end

        function ok = dualInfCheck(Adx, l, u, eps)
            INF = osqp.Solver.OSQP_INFTY;
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

        function new_rho = computeNewRho(xs, zs, ys, Ps, qs, As, rho, tol, ~, m)
            % Compute rho estimate matching C code compute_rho_estimate.
            % Normalizes primal/dual residuals before computing ratio.
            DTOL = osqp.Solver.DIVISION_TOL;
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
            rho_estimate = max(min(rho_estimate, osqp.Solver.RHO_MAX), ...
                               osqp.Solver.RHO_MIN);

            if rho_estimate > rho * tol || rho_estimate < rho / tol
                new_rho = rho_estimate;
            else
                new_rho = rho;
            end
        end

        function results = polishSolution(results, P_triu, q, A, l, u, s)
            n = numel(q);
            m = size(A, 1);
            x = results.x;
            y = results.y;
            Pfull = P_triu + P_triu' - diag(diag(P_triu));
            INF = osqp.Solver.OSQP_INFTY;

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
                    pri_res_pol = norm(Ax_pol - min(max(Ax_pol, l), u), inf);
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
