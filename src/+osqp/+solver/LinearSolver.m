classdef (Abstract) LinearSolver < handle
    % OSQP.LINEARSOLVER  Abstract base class for KKT linear system solvers.
    %
    %   Subclasses must implement:
    %     solve(obj, b)       — solve K*x = b using cached factorization
    %     refactorize(obj, K) — update factorization with new K
    %
    %   The concrete mldivide method enables backslash syntax:
    %     x = solver \ b;
    %
    %   Static factory method buildKKT assembles the KKT matrix from
    %   problem data. Subclasses call this, then factorize.

    methods (Abstract)
        x = solve(obj, b)
        refactorize(obj, K)
    end

    methods
        function x = mldivide(obj, b)
            % MLDIVIDE  Enable backslash syntax: x = solver \ b
            x = obj.solve(b);
        end
    end

    methods (Static)
        function K = buildKKT(P_triu, A, rho_vec, sigma, n, m)
            % BUILDKKT  Assemble the quasi-definite KKT matrix (upper triangular).
            %
            %   K = triu([P + sigma*I,   A'        ])
            %           ([0,         -diag(1./rho) ])
            %
            %   Input P_triu is upper triangular; output K is upper triangular.
            %   This matches the C OSQP form_KKT which stores only the upper
            %   triangle. MATLAB decomposition('ldl') and QDLDL both read only
            %   the upper triangle, so this is correct for all backends.
            Ktl = P_triu + sigma * speye(n);
            if m > 0
                Kbr = -spdiags(1 ./ rho_vec, 0, m, m);
                K   = [Ktl, A'; sparse(m, n), Kbr];
            else
                K = Ktl;
            end
        end
    end
end
