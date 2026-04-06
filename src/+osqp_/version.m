function out = version()
% OSQP.VERSION  Return the OSQP version string.
%
%   v = osqp.version()

    if osqp.CInterface.isMexAvailable()
        out = osqp.CInterface.version();
    else
        out = '1.0.0-matlab';
    end
end
