# OSQP-MATLAB Roadmap

Date: 2026-05-24  
Owner: Jason Nicholson  
Status: Draft for discussion

## 1) Goal

Build a stable, testable, and deployment-ready OSQP-MATLAB fork that:

1. Has a reliable OSQP C API v1 interface.
2. Supports Simulink simulation and code generation workflows.
3. Supports embedded deployment experiments.
4. Has a reliable pure MATLAB solver path.

## 2) Project Principles

1. Correctness first, then performance, then convenience.
2. Keep C backend and MATLAB backend behavior aligned where practical.
3. Use explicit acceptance criteria for every milestone.
4. Minimize duplicate code and hard-coded constants.
5. Keep Simulink integration maintainable and easy to validate.

## 3) Scope and Non-Goals

### In Scope

1. C backend hardening around osqp.CInterface.
2. Pure MATLAB solver completion and validation.
3. Simulink integration using either:
   - MATLAB Level-2 S-Function + TLC/hook workflow, or
   - C S-Function wrapper (direct/native path), or
   - A staged hybrid (recommended).
4. Embedded-oriented performance exploration (assembly-aware and sparse kernels).

### Out of Scope (for initial roadmap horizon)

1. Immediate upstream merge to osqp/osqp-matlab.

## 4) Decision: MATLAB System Object vs C S-Function

This is the key architecture decision for Simulink integration.

### Option A: MATLAB System object

Pros:
1. Fastest development iteration.
2. Easier MATLAB-side debugging.
3. Good for simulation-first usage.

Cons:
1. Code generation constraints can be strict and version-sensitive.
2. Often ends up needing external C integration plumbing anyway.
3. Harder to guarantee deterministic low-overhead behavior for embedded targets.

### Option B: C S-Function wrapper

Pros:
1. Best path for deterministic runtime and embedded workflows.
2. Closer to OSQP C codegen artifacts.
3. Reduces MATLAB interpreter overhead in simulation loops.

Cons:
1. Higher development complexity.
2. Higher maintenance burden for build and toolchain integration.
3. Longer time to first working block.

### Option C: Staged hybrid (recommended)

1. Stage 1: deliver MATLAB Level-2 S-Function path to unblock simulation and codegen experimentation.
2. Stage 2: implement C S-Function for performance-critical and embedded deployment paths.
3. Keep block interface stable so models migrate without rewiring.

## 5) Phased Roadmap

### Phase 0: Baseline and cleanup

Objective: stabilize foundation before feature expansion.

Deliverables:
1. Remove duplicate constants/settings paths and centralize configuration flow.
2. Confirm backend-specific unit test runners are green.

Acceptance criteria:
1. Unit suites pass for C and MATLAB backends in isolated runs.
2. No known duplicated constant definitions in active paths.
3. README and docs agree on backend capabilities and limitations.

### Phase 1: C interface hardening

Objective: make osqp.CInterface predictable and robust.

Deliverables:
1. Complete settings translation/validation path from MATLAB options to C defaults.
2. Expand tests for setup, solve, update, warm/cold start, and status mapping.
3. Add regression tests for prior configuration bugs.
4. Add end-to-end CInterface test runs against `osqp_problems_collection`.

Acceptance criteria:
1. C backend tests pass across supported OS/compiler matrix.
2. Error messages for invalid options are clear and actionable.
3. Backward-incompatible settings changes are documented.
4. CInterface end-to-end runs complete on the problem collection with pass/fail reporting.

### Phase 2: Pure MATLAB backend completion

Objective: close functional and numerical gaps in pure MATLAB solver.

Deliverables:
1. Validate and benchmark linear solver modes:
   - qdldl_c default when available,
   - MATLAB LDL path,
   - pure MATLAB QDLDL fallback.
2. Strengthen parity tests versus C backend on representative QP families.
3. Add performance dashboard scripts for repeatable benchmark runs.
4. Add end-to-end MATLAB solver test runs against `osqp_problems_collection`.

Acceptance criteria:
1. Solver status and solution quality pass parity thresholds on test corpus.
2. Fallback behavior is deterministic when optional binaries are unavailable.
3. Documented performance tradeoffs by solver mode.
4. MATLAB backend end-to-end runs complete on the problem collection with pass/fail reporting.

### Phase 3: Simulink integration MVP

Objective: ship a usable Simulink path with clear boundaries.

Deliverables:
1. Reproduce and modernize the known Level-2 S-Function pattern from historical branches:
   - block files,
   - TLC integration,
   - RTW hook integration.
2. Provide at least one maintained example model showing:
   - simulation path,
   - Simulink code generation path.
3. Define and freeze block I/O contract (dimensions, updates, warm start behavior).

Acceptance criteria:
1. Example model simulates without manual path hacks.
2. Code generation completes for at least one supported target/toolchain profile.
3. Known limitations and required hooks are documented.

### Phase 4: Simulink production path (C S-Function)

Objective: enable lower-overhead deployment path for real-time use.

Deliverables:
1. Implement C S-Function wrapper aligned with stable block interface from Phase 3.
2. Minimize duplication by sharing logic with osqp_mex bridge where safe.
3. Add integration tests for repeated solve/update in closed-loop simulation.

Acceptance criteria:
1. C S-Function benchmark shows measurable overhead reduction vs Level-2 path.
2. Code generation and build procedure is reproducible on supported environments.
3. Fallback to Phase 3 path remains available for unsupported targets.

### Phase 5: Embedded optimization track (ongoing)

Objective: de-risk embedded deployment on embedded targets.

Deliverables:
1. Baseline generated C solver on representative embedded build chain.
2. Profile hotspots with emphasis on sparse matrix-vector and linear solve kernels.
3. Evaluate SparseBLAS/GraphBLAS-style acceleration options where integration cost is justified.

Acceptance criteria:
1. End-to-end benchmark harness available for target profile.
2. Optimization changes gated by reproducible before/after measurements.
3. Memory and timing budgets tracked per benchmark scenario.

## 6) Workstreams (Single Maintainer)

1. Core API and correctness: CInterface and settings mapping.
2. Numerical backend quality: MATLAB solver + linsys variants.
3. Simulink integration: block, TLC, coder hooks, examples.
4. Embedded performance: profiling, toolchain integration, kernel optimization.
5. Quality engineering: CI matrix, regression suites, compatibility docs.

## 7) Risks and Mitigations

1. MATLAB/Simulink version fragility:
   - Mitigation: define explicit support matrix and test only supported versions first.
2. Toolchain and path brittleness in Simulink code generation:
   - Mitigation: automated environment checks and setup scripts; document hook file requirements.
3. Divergence from upstream osqp-matlab:
   - Mitigation: isolate fork-specific features behind clear module boundaries.
4. Premature optimization on embedded targets:
   - Mitigation: profile first, optimize only high-impact hotspots.

## 8) Upstream Strategy Gate

Revisit merge strategy after Phase 3.

Decision criteria:
1. How invasive are API and architecture differences?
2. How much maintenance burden is reduced by upstreaming?
3. Is the Simulink path acceptable to upstream maintainers?

Outcomes:
1. Upstream selected components (incremental PRs), or
2. Maintain long-lived fork with periodic upstream sync policy.

## 9) Near-Term Execution Plan

1. Finalize Phase 0 cleanup list and close remaining settings/constants inconsistencies.
2. Lock CInterface acceptance tests and regression tests, including end-to-end runs on `osqp_problems_collection`.
3. Lock MATLAB solver acceptance tests and parity checks, including end-to-end runs on `osqp_problems_collection`.
4. Prototype Simulink MVP using historical Level-2 S-Function + TLC/hook pattern.
5. Produce one end-to-end example with simulation and codegen walkthrough.
6. Hold architecture review and confirm Stage 2 C S-Function scope.

## 10) Definition of Done for Roadmap v1

This roadmap draft is considered accepted when:

1. Each phase has clear entry/exit criteria and lightweight review triggers.
2. End-to-end validation exists for both CInterface and MATLAB backends using `osqp_problems_collection`.
3. Support matrix (MATLAB/OS/toolchain) is explicitly listed.
4. Simulink architecture path (staged hybrid or alternative) is approved.
5. The next two phases have clear, prioritized task lists sized for one maintainer.
