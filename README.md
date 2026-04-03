# OSQP-MATLAB — MATLAB Interface for OSQP v1.0.0

MATLAB wrapper for [OSQP](https://osqp.org/) v1.0.0: the Operator Splitting QP Solver.

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

solver = osqp;
solver.setup(P, q, A, l, u, 'verbose', false);
res = solver.solve();

fprintf('x = [%.4f, %.4f]\n', res.x(1), res.x(2));
delete(solver);
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

Code generation options: `'parameters'` (`'vectors'`, `'matrices'`, or `'none'`), `'force'` (overwrite), `'float'` (use single precision), `'printing'` (enable printing), `'profiling'` (enable profiling), `'interrupt'` (enable interrupt), `'derivatives'` (enable derivatives).

### Adjoint Derivatives

| Method | Description |
|--------|-------------|
| `adjoint_derivative_compute(dx, dy)` | Compute adjoint derivatives given loss gradients |
| `adjoint_derivative_get_mat()` | Get derivative matrices `[dP, dA]` |
| `adjoint_derivative_get_vec()` | Get derivative vectors `[dq, dl, du]` |

### Utility (Static Methods)

| Method | Description |
|--------|-------------|
| `osqp.version()` | OSQP version string |
| `osqp.constant(name)` | Query an OSQP constant |
| `osqp.capabilities()` | Capability bitmask |
| `has_capability(cap)` | Check for a specific capability constant |
| `get_dimensions()` | Returns `[m, n]` |

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
| `alpha` | 1.6 | ADMM relaxation parameter |
| `check_termination` | 25 | Check termination every N iterations |

## Running Tests

```matlab
setupOSQPdevelopmentPath   % if not already done
run_osqp_tests             % run all 37 unit tests
run_osqp_tests(true)       % run tests with HTML coverage report
```

## Running Examples

```matlab
setupOSQPdevelopmentPath   % if not already done
run_osqp_examples          % run all example scripts
```

## Project Structure

```
osqp-matlab/
├── setupOSQPdevelopmentPath.m   % adds src/ and utils/ to the MATLAB path
├── src/                         % packaged source files
│   ├── osqp.m                   % main OSQP class
│   ├── make_osqp.m              % build script
│   ├── osqp_mex.cpp/.hpp        % MEX gateway
│   ├── CMakeLists.txt           % CMake build configuration
│   └── codegen/                 % code generation support
├── examples/                    % example scripts
├── unittests/                   % unit test suite
└── utils/                       % development utilities
```

## Migration from v0.6.x

Key breaking changes from OSQP ≤ 0.6.x:

- **Settings renamed**: `warm_start` → `warm_starting`, `polish` → `polishing`
- **Linear system solvers**: `'qdldl'` → `'direct'`, `'mkl pardiso'` → `'direct'` (select backend at build time)
- **Code generation**: Uses native `osqp_codegen()` instead of template rendering
- **No submodule**: OSQP is fetched automatically via CMake FetchContent
- **New features**: Adjoint derivatives, cold start, CUDA backend, capability queries

## License

Apache 2.0. See [LICENSE](LICENSE) for details.
