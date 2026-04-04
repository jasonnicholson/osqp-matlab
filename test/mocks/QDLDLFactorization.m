classdef QDLDLFactorization
    % QDLDLFACTORIZATION  Mock for testing QDLDLSolver without qdldl.m.
    %
    %   Provides the same interface as the real QDLDLFactorization but uses
    %   MATLAB's built-in decomposition for the actual factorization.

    properties (Access = private)
        decomp  % MATLAB decomposition object
    end

    methods (Static)
        function obj = fromMatrix(K)
            obj = QDLDLFactorization();
            if ~issparse(K)
                K = sparse(K);
            end
            Kfull = K + triu(K, 1)';
            obj.decomp = decomposition(Kfull, 'ldl', 'upper');
        end
    end

    methods
        function x = mldivide(obj, b)
            x = obj.decomp \ b;
        end
    end
end
