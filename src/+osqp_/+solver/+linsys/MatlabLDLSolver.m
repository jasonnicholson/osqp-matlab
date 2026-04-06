classdef MatlabLDLSolver < osqp.LinearSolver
    % OSQP.LINSYS.MATLABLDLSOLVER  KKT solver using MATLAB's decomposition(K,'ldl').
    %
    %   Uses MATLAB's built-in LDL factorization via the decomposition
    %   object.  Handles the indefinite (quasi-definite) KKT matrices that
    %   arise in OSQP.  The decomposition object caches the factorization
    %   for fast repeated solves.
    %
    %   Example:
    %     K = osqp.LinearSolver.buildKKT(P, A, rho_vec, sigma, n, m);
    %     solver = osqp.linsys.MatlabLDLSolver(K);
    %     x = solver \ b;

    properties (Access = private)
        dK   % decomposition object
    end

    methods
        function obj = MatlabLDLSolver(K)
            % MATLABLDLSOLVER  Construct from upper-triangular KKT matrix K.
            arguments
                K (:,:) {mustBeNumeric}
            end
            obj.dK = decomposition(K, 'ldl', 'upper');
        end

        function x = solve(obj, b)
            % SOLVE  Solve K*x = b using cached LDL decomposition.
            x = obj.dK \ b;
        end

        function refactorize(obj, K)
            % REFACTORIZE  Rebuild the decomposition with new K.
            obj.dK = decomposition(K, 'ldl', 'upper');
        end
    end
end
