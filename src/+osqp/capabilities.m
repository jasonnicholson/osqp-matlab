function out = capabilities()
% OSQP.CAPABILITIES  Return the OSQP capability bitmask.
%
%   cap = osqp.capabilities()

    if osqp.CInterface.isMexAvailable()
        out = osqp.CInterface.capabilities();
    else
        out = 0;
    end
end
