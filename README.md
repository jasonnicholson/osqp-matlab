# OSQP-MATLAB — A Fork of MATLAB Interface for OSQP - UNDER HEAVY CONSTRUCTION

I forked this off from the official [OSQP-MATLAB](https://github.com/osqp/osqp-matlab) repository to build an updated OSQP C interface and implement a pure MATLAB version.

**Don't expect anything to work right now. I am ripping it apart and rebuilding it from the ground up. 2026-May-04**


MATLAB wrapper for [OSQP](https://osqp.org/), the Operator Splitting QP Solver and a pure-MATLAB implementation for environments without a C compiler.

## Roadmap

Project roadmap: [docs/roadmap.md](docs/roadmap.md)

There are two backends available:
`solver = osqp(backend='c')`
`solver = osqp(backend='matlab')`

OSQP solves convex quadratic programs of the form:

```
minimize        0.5 x' P x + q' x

subject to      l <= A x <= u
```

where `x in R^n` is the optimization variable. `P in S^n_+` is a positive semidefinite matrix, `q in R^n`, `A in R^{m x n}`, and `l, u in R^m U {±inf}^m`.

## Prerequisites

- MATLAB R2018a or later
- A supported C/C++ compiler (see `mex -setup`)
- CMake 3.15 or later (must be on the system PATH)

### Optional

- Intel MKL — for the MKL algebra backend
- NVIDIA CUDA Toolkit — for the CUDA algebra backend

## Installation

Clone the repository:

```bash
git clone https://github.com/osqp/osqp-matlab.git
cd osqp-matlab
```

Set up the development path and build the MEX interface from within MATLAB:

```matlab
setupOSQPdevelopmentPath   % adds src/ and utils/ to the MATLAB path
make_osqp                  % default (builtin/QDLDL algebra backend)
```

The build system uses CMake `FetchContent` to automatically download OSQP v1.0.0 and its dependencies (QDLDL). No git submodules are needed.

Build options:

```matlab
make_osqp('algebra', 'mkl')     % Intel MKL algebra backend
make_osqp('algebra', 'cuda')    % NVIDIA CUDA algebra backend
make_osqp('-verbose')           % print full compiler output
make_osqp('clean')              % remove compiled MEX files
make_osqp('purge')              % clean + remove build directory
```

## Quick Start

```matlab
P = sparse([4 1; 1 2]);
q = [1; 1];
A = sparse([1 1; 1 0; 0 1]);
l = [1; 0; 0];
u = [1; 0.7; 0.7];

solver = osqp;                                  % C backend (default)
solver.setup(P, q, A, l, u, 'verbose', false);
res = solver.solve();

fprintf('x = [%.4f, %.4f]\n', res.x(1), res.x(2));
delete(solver);
```

### Pure-MATLAB Backend

A pure-MATLAB ADMM solver is available for environments without a C compiler:

```matlab
solver = osqp(backend="matlab");
solver.setup(P, q, A, l, u, 'verbose', false);
res = solver.solve();
```

## API Reference

### Core Methods

| Method | Description |
|--------|-------------|
| `setup(P, q, A, l, u, varargin)` | Configure solver with problem data and optional settings |
| `solve()` | Solve the QP, returns struct with `x`, `y`, `info` |
| `update('q', q_new, 'l', l_new, 'u', u_new)` | Update problem vectors |
| `update('Px', Px_new, 'Px_idx', idx)` | Update P matrix values |
| `update('Ax', Ax_new, 'Ax_idx', idx)` | Update A matrix values |
| `warm_start('x', x0, 'y', y0)` | Set warm-start primal/dual variables |
| `cold_start()` | Reset solver iterates to zero |

### Settings

| Method | Description |
|--------|-------------|
| `osqp.default_settings()` | Get default settings structure (static) |
| `current_settings()` | Get current solver settings |
| `update_settings(varargin)` | Update solver settings (key-value pairs or struct) |

### Code Generation

| Method | Description |
|--------|-------------|
| `codegen(dir, varargin)` | Generate embeddable C code for the current problem |

Code generation options: `'parameters'` (`'vectors'` or `'matrices'`), `'force_rewrite'` (overwrite), `'float_type'` (use single precision), `'printing_enable'` (enable printing), `'profiling_enable'` (enable profiling), `'interrupt_enable'` (enable interrupt), `'derivatives_enable'` (enable derivatives), `'prefix'` (symbol prefix).

### Adjoint Derivatives

| Method | Description |
|--------|-------------|
| `adjoint_derivative_compute(dx, dy)` | Compute adjoint derivatives given loss gradients |
| `adjoint_derivative_get_mat()` | Get derivative matrices `[dP, dA]` |
| `adjoint_derivative_get_vec()` | Get derivative vectors `[dq, dl, du]` |

### Utility

| Method | Description |
|--------|-------------|
| `osqp.version()` | OSQP version string |
| `osqp.constant(name)` | Query an OSQP constant |
| `osqp.capabilities()` | Capability bitmask |
| `has_capability(cap)` | Check for a specific capability constant |
| `get_dimensions()` | Returns `[n, m]` |

## Backends

| Feature | C (`osqp()`) | MATLAB (`osqp(backend="matlab")`) |
|---------|:---:|:---:|
| `setup` / `solve` / `update` | Yes | Yes |
| `warm_start` / `cold_start` | Yes | Yes |
| `update_settings` | Yes | Yes |
| `codegen` | Yes | No |
| `adjoint_derivative_*` | Yes | No |
| Linear solvers | QDLDL (C) | MATLAB LDL / QDLDL (MATLAB) |
| Requires MEX compilation | Yes | No |

## Key Settings

| Setting | Default | Description |
|---------|---------|-------------|
| `rho` | 0.1 | ADMM penalty parameter |
| `sigma` | 1e-6 | ADMM penalty parameter |
| `max_iter` | 4000 | Maximum iterations |
| `eps_abs` | 1e-3 | Absolute convergence tolerance |
| `eps_rel` | 1e-3 | Relative convergence tolerance |
| `eps_prim_inf` | 1e-4 | Primal infeasibility tolerance |
| `eps_dual_inf` | 1e-4 | Dual infeasibility tolerance |
| `polishing` | false | Enable solution polishing |
| `warm_starting` | true | Enable warm starting |
| `verbose` | true | Enable printing |
| `linsys_solver` | `'direct'` | Linear system solver (`'direct'`/`'indirect'`) |
| `scaling` | 10 | Number of scaling iterations (0 = disabled) |
| `scaled_termination` | false | Check convergence in scaled space (true) or unscaled (false, matches C default) |
| `alpha` | 1.6 | ADMM relaxation parameter |
| `check_termination` | 25 | Check termination every N iterations |

## Running Tests

```matlab
setupOSQPdevelopmentPath          % if not already done
run_osqp_unit_tests               % run all unit tests
run_osqp_unit_tests(true)         % run tests with HTML coverage report
run_osqp_cinterface_unit_tests    % run CInterface-focused unit tests only
run_osqp_matlab_unit_tests        % run pure-MATLAB solver unit tests only
```

## Running Examples

```matlab
setupOSQPdevelopmentPath          % if not already done
run_osqp_examples                 % run all example scripts
```

## Project Structure

```
osqp-matlab/
├── setupOSQPdevelopmentPath.m   % adds src/ and utils/ to the MATLAB path
├── src/
│   ├── osqp.m                   % factory function (selects backend)
│   ├── make_osqp.m              % build script for MEX interface
│   ├── osqp_mex.cpp/.hpp        % MEX gateway (C ↔ MATLAB bridge)
│   ├── CMakeLists.txt           % CMake build configuration
│   ├── +osqp/                   % OSQP package namespace
│   │   ├── CInterface.m         % C backend wrapper over osqp_mex
│   │   ├── Solver.m             % Pure-MATLAB ADMM backend
│   │   ├── SolverOptions.m      % Shared options object (both backends)
│   │   ├── capabilities.m       % C backend capabilities query
│   │   ├── default_settings.m   % Package-level default settings
│   │   ├── version.m            % OSQP version helper
│   │   ├── private/             % Shared validation/helpers (e.g., validateData)
│   │   └── +solver/             % MATLAB solver internals
│   │       ├── Constants.m      % Algorithm constants for MATLAB backend
│   │       ├── LinearSolver.m   % KKT assembly / abstract linear solver utilities
│   │       ├── StatusType.m     % Status enums for MATLAB backend
│   │       ├── ErrorType.m      % Error enums for MATLAB backend
│   │       ├── LinSysSolverType.m   % Linear solver enum definitions
│   │       ├── PrecondType.m    % Preconditioner enum definitions
│   │       ├── PolishStatusType.m   % Polishing status enums
│   │       └── +linsys/         % MATLAB backend linear solver implementations
│   │           ├── MatlabLDLSolver.m
│   │           ├── QDLDLSolver.m
│   │           └── QDLDLCSolver.m
│   └── codegen/                 % Code generation support
├── examples/                    % Example scripts
├── test/
│   ├── run_osqp_unit_tests.m            % all tests or backend-scoped tests
│   ├── run_osqp_cinterface_unit_tests.m % CInterface-focused unit suite
│   ├── run_osqp_matlab_unit_tests.m     % pure-MATLAB unit suite
│   ├── unit/                    % Unit test suite (matlab.unittest)
│   │   ├── basic_tests.m, ...   % CInterface-focused tests
│   │   └── solver_tests.m, ...  % MATLAB backend tests
│   ├── integration/             % Integration tests (script-based)
│   └── mocks/                   % Test doubles
├── utils/                       % Development utilities
└── build/                       % CMake build output (auto-generated)
```

### Backend-Oriented Layout

- CInterface backend (default via `osqp()`):
	- Entry point: `src/+osqp/CInterface.m`
	- Native bridge: `src/osqp_mex.cpp`, `src/osqp_mex.hpp`, `src/make_osqp.m`
	- C-only APIs: codegen, adjoint derivatives, capability/constant queries
- Pure-MATLAB backend (via `osqp(backend="matlab")`):
	- Entry point: `src/+osqp/Solver.m`
	- Core internals: `src/+osqp/+solver/*`
	- Linear solvers: `src/+osqp/+solver/+linsys/*`
- Shared across both backends:
	- Factory: `src/osqp.m`
	- Settings model: `src/+osqp/SolverOptions.m`
	- Common data validation: `src/+osqp/private/validateData.m`

## Migration from v0.6.x

Key breaking changes from OSQP ≤ 0.6.x:

- **Factory function**: `osqp` is now a function that returns a backend object (`osqp.CInterface` or `osqp.Solver`), not a class
- **Dual backends**: Use `osqp()` for C backend (default), `osqp(backend="matlab")` for pure-MATLAB
- **Package namespace**: All source files live under the `+osqp` package (`osqp.Solver`, `osqp.Options`, etc.)
- **Settings renamed**: `warm_start` → `warm_starting`, `polish` → `polishing`
- **Linear system solvers**: `'qdldl'` → `'direct'`, `'mkl pardiso'` → `'direct'` (select backend at build time)
- **Code generation**: Uses native `osqp_codegen()` instead of template rendering
- **No submodule**: OSQP is fetched automatically via CMake FetchContent
- **New features**: Adjoint derivatives, cold start, CUDA backend, capability queries

## License

Apache 2.0. See [LICENSE](LICENSE) for details.
