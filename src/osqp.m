function solver = osqp(args)
% OSQP  Create an OSQP quadratic programming solver.
%
%   solver = osqp()
%       Returns an osqp.CInterface (C MEX backend). Errors if MEX is
%       unavailable.
%
%   solver = osqp(backend="matlab")
%       Returns an osqp.Solver (pure-MATLAB ADMM backend).
%
%   solver = osqp(backend="c")
%       Same as osqp() — explicit C backend.
%
%   The returned object supports: setup, solve, update, warm_start,
%   cold_start, current_settings, update_settings, get_dimensions.
%   CInterface additionally supports: codegen, adjoint_derivative_*,
%   has_capability, version, constant, capabilities.
%
%   See also: osqp.CInterface, osqp.Solver, osqp.Options

    arguments
        args.backend (1,:) char {mustBeMember(args.backend, {'c', 'matlab'})} = 'c';
    end

    switch args.backend
        case 'matlab'
            solver = osqp.Solver();
        case 'c'
            if ~osqp.CInterface.isMexAvailable()
                error('OSQP:osqp', 'C backend requested but osqp_mex is not available.');
            end
            solver = osqp.CInterface();
    end
end
