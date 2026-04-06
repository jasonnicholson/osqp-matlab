classdef QDLDLCSolver < osqp.LinearSolver
    % OSQP.LINSYS.QDLDLCSOLVER  KKT solver using C-MEX QDLDL factorization.
    %
    %   Wraps qdldl_c_factor_mex (compiled C QDLDL) for fast LDL
    %   factorization of quasi-definite KKT matrices.
    %
    %   Requires qdldl_c_factor_mex MEX file to be on the MATLAB path.
    %
    %   Example:
    %     K = osqp.LinearSolver.buildKKT(P, A, rho_vec, sigma, n, m);
    %     solver = osqp.linsys.QDLDLCSolver(K);
    %     x = solver \ b;

    properties (Access = private)
        perm      % permutation vector (or empty)
        L_unit    % unit lower triangular factor (L + I)
        Dinv_vec  % inverse diagonal vector
    end

    methods
        function obj = QDLDLCSolver(K)
            % QDLDLCSOLVER  Construct from symmetric matrix K.
            arguments
                K (:,:) {mustBeNumeric}
            end
            if ~exist('qdldl_c_factor_mex', 'file')
                error('osqp:linsys:QDLDLCSolver', ...
                    'qdldl_c_factor_mex not found. Build or add to path.');
            end
            obj.factorize(K);
        end

        function x = solve(obj, b)
            % SOLVE  Solve K*x = b using cached C-QDLDL factorization.
            if isempty(obj.perm)
                tmp = b;
            else
                tmp = b(obj.perm);
            end

            tmp = obj.L_unit \ tmp;
            tmp = tmp .* obj.Dinv_vec;
            tmp = obj.L_unit' \ tmp;

            if isempty(obj.perm)
                x = tmp;
            else
                x = b;
                x(obj.perm) = tmp;
            end
        end

        function refactorize(obj, K)
            % REFACTORIZE  Rebuild the C-QDLDL factorization with new K.
            obj.factorize(K);
        end
    end

    methods (Access = private)
        function factorize(obj, K)
            n = size(K, 1);
            if ~issparse(K)
                K = sparse(K);
            end
            Ktriu = triu(K);

            % AMD permutation for fill-reducing ordering
            obj.perm = symamd(Ktriu);
            obj.perm = obj.perm(:);

            Kfull = Ktriu + Ktriu' - spdiags(diag(Ktriu), 0, n, n);
            Kperm = triu(Kfull(obj.perm, obj.perm));

            [L, obj.Dinv_vec] = qdldl_c_factor_mex(Kperm);
            obj.L_unit = L + speye(n);
        end
    end
end
