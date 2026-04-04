classdef MATLABLDLSolver < handle
% MATLABLDLSolver  Thin wrapper around MATLAB's built-in ldl().
%
%   Provides a cached LDL factorization for repeated solves via mldivide.
%   Uses vector permutation output for reduced storage.
%
%   MATLAB's ldl() computes K(p,p) = L * D * L' where L is unit lower
%   triangular (sparse) and D is block-diagonal (1x1 or 2x2 Bunch-Kaufman
%   blocks for indefinite matrices such as the OSQP KKT system).

    properties (Access = private)
        L       % unit lower triangular sparse factor
        D       % block-diagonal factor (sparse)
        p       % permutation vector: K(p,p) = L*D*L'
    end

    methods
        function obj = MATLABLDLSolver(K)
            [obj.L, obj.D, obj.p] = ldl(K, 'vector');
        end

        function refactorize(obj, K)
            [obj.L, obj.D, obj.p] = ldl(K, 'vector');
        end

        function x = mldivide(obj, b)
            % Solve K*x = b  where  K(p,p) = L * D * L'
            c      = obj.L  \ b(obj.p);
            c      = obj.D  \ c;
            c      = obj.L' \ c;
            x      = b;
            x(obj.p) = c;
        end
    end
end
