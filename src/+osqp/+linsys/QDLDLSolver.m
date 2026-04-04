classdef QDLDLSolver < osqp.LinearSolver
    % OSQP.LINSYS.QDLDLSOLVER  KKT solver using pure-MATLAB QDLDL.
    %
    %   Wraps QDLDLFactorization.fromMatrix() from the qdldl.m package.
    %   Requires qdldl.m/src to be on the MATLAB path.
    %
    %   Example:
    %     K = osqp.LinearSolver.buildKKT(P, A, rho_vec, sigma, n, m);
    %     solver = osqp.linsys.QDLDLSolver(K);
    %     x = solver \ b;

    properties (Access = private)
        factor   % QDLDLFactorization handle
        K_       % cached matrix for refactor
    end

    methods
        function obj = QDLDLSolver(K)
            % QDLDLSOLVER  Construct from symmetric matrix K.
            arguments
                K (:,:) {mustBeNumeric}
            end
            if ~exist('QDLDLFactorization', 'class')
                error('osqp:linsys:QDLDLSolver', ...
                    'QDLDLFactorization not found. Add qdldl.m/src to path.');
            end
            obj.K_ = K;
            obj.factor = QDLDLFactorization.fromMatrix(K);
        end

        function x = solve(obj, b)
            % SOLVE  Solve K*x = b using cached QDLDL factorization.
            x = obj.factor \ b;
        end

        function refactorize(obj, K)
            % REFACTORIZE  Rebuild the QDLDL factorization with new K.
            obj.K_ = K;
            obj.factor = QDLDLFactorization.fromMatrix(K);
        end
    end
end
