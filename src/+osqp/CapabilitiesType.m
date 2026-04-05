classdef CapabilitiesType
    % OSQP.CAPABILITIESTYPE  Flags for compiled solver capabilities.
    %   Mirrors osqp_capabilities_type in osqp_api_constants.h.
    %

    enumeration
        OSQP_CAPABILITY_DIRECT_SOLVER   (1)   % 0x01 — direct linear solver present
        OSQP_CAPABILITY_INDIRECT_SOLVER (2)   % 0x02 — indirect linear solver present
        OSQP_CAPABILITY_CODEGEN         (4)   % 0x04 — code generation present
        OSQP_CAPABILITY_UPDATE_MATRICES (8)   % 0x08 — problem matrices can be updated
        OSQP_CAPABILITY_DERIVATIVES     (16)  % 0x10 — solution derivatives available
    end
end
