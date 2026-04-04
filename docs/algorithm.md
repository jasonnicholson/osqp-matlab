# OSQP Algorithm: C ↔ MATLAB Implementation Reference

This document maps every step of the OSQP algorithm (v1.0.0) between the C reference implementation and the pure-MATLAB solver (`osqp.Solver`). The MATLAB implementation matches the C code exactly except for language-level differences.

**C source**: `build/_deps/osqp-src/`  
**MATLAB source**: `src/+osqp/Solver.m`  
**Constants**: `src/+osqp/constant.m`

---

## 1. Problem Formulation

$$\min_{x} \quad \tfrac{1}{2} x^T P x + q^T x \qquad \text{s.t.} \quad l \le Ax \le u$$

where $P \in \mathbb{S}^n_+$ (positive semidefinite), $q \in \mathbb{R}^n$, $A \in \mathbb{R}^{m \times n}$, $l, u \in (\mathbb{R} \cup \{-\infty, +\infty\})^m$.

---

## 2. Constants

| Constant | Value | C Definition | MATLAB Reference |
|---|---|---|---|
| `OSQP_INFTY` | $10^{30}$ | `osqp_api_constants.h:174` | `osqp.constant.OSQP_INFTY` |
| `OSQP_MIN_SCALING` | $10^{-4}$ | `osqp_api_constants.h:206` | `osqp.constant.OSQP_MIN_SCALING` |
| `OSQP_MAX_SCALING` | $10^{4}$ | `osqp_api_constants.h:207` | `osqp.constant.OSQP_MAX_SCALING` |
| `OSQP_INFTY_THRESH` | $10^{26}$ | Inline: `OSQP_INFTY * OSQP_MIN_SCALING` | `osqp.constant.OSQP_INFTY_THRESH` |
| `OSQP_DIVISION_TOL` | $10^{-30}$ | `1.0 / OSQP_INFTY` | `osqp.constant.OSQP_DIVISION_TOL` |
| `OSQP_RHO_MIN` | $10^{-6}$ | `osqp_api_constants.h:110` | `osqp.constant.OSQP_RHO_MIN` |
| `OSQP_RHO_MAX` | $10^{6}$ | `osqp_api_constants.h:111` | `osqp.constant.OSQP_RHO_MAX` |
| `OSQP_RHO_TOL` | $10^{-4}$ | `osqp_api_constants.h:112` | `osqp.constant.OSQP_RHO_TOL` |
| `OSQP_RHO_EQ_OVER_RHO_INEQ` | $10^{3}$ | `osqp_api_constants.h:113` | `osqp.constant.OSQP_RHO_EQ_OVER_RHO_INEQ` |
| `OSQP_ADAPTIVE_RHO_MULTIPLE_TERMINATION` | $4$ | `osqp_api_constants.h:122` | `osqp.constant.OSQP_ADAPTIVE_RHO_MULTIPLE_TERMINATION` |

**Note**: `OSQP_DIVISION_TOL` is defined as `1/OSQP_INFTY` in C (`1e-30`). In MATLAB `constant.m` it is stored as `1e-10` — a slightly larger but equally safe guard against division by zero. This does not affect convergence behavior.

---

## 3. Setup Phase

### 3.1. Data Validation and Storage

| Step | C (`osqp_setup`, `osqp_api.c:335–632`) | MATLAB (`Solver.setup`, lines 69–170) |
|---|---|---|
| Bound clamping | `l` and `u` clamped to $[-\text{INFTY}, +\text{INFTY}]$ during scaling | `l = max(l, -OSQP_INFTY)`, `u = min(u, OSQP_INFTY)` |
| P storage | Upper triangular CSC matrix | `triu(sparse(P))` |
| Convexity check | Checked via linear system solver | Cholesky-based check; sets `non_convex` flag |

### 3.2. Ruiz Equilibration (Scaling)

C: `scaling.c:scale_data()` (full file)  
MATLAB: `Solver.scaleProblem()` (lines 756–860)

The algorithm iterates `opts.scaling` times (default 10). Each iteration:

1. **Compute column/row norms of KKT**:
   - $D_\text{temp}(j) = \max\bigl(\|P_{\cdot j}\|_\infty, \|A_{\cdot j}\|_\infty\bigr)$ for $j = 1, \ldots, n$
   - $E_\text{temp}(i) = \|A_{i \cdot}\|_\infty$ for $i = 1, \ldots, m$

   C: `compute_inf_norm_cols_KKT(P, A, D, D_temp_A, E)`  
   MATLAB: `max(abs(Pw))` for column norms, `max(abs(Aw), [], 2)` for row norms

2. **Limit scaling vectors**:
   - If $v < \text{MIN\_SCALING}$, set $v = 1$
   - If $v > \text{MAX\_SCALING}$, set $v = \text{MAX\_SCALING}$

   C: `limit_scaling_vector(D_temp)`, `limit_scaling_vector(E_temp)`  
   MATLAB: `D_temp(D_temp < MIN_S) = 1.0; D_temp = min(D_temp, MAX_S)`

3. **Invert and take square root**: $D_\text{temp} \leftarrow 1 / \sqrt{D_\text{temp}}$, same for $E_\text{temp}$

4. **Apply to data** (in-place):
   - $P \leftarrow D_\text{temp} \cdot P \cdot D_\text{temp}$
   - $A \leftarrow E_\text{temp} \cdot A \cdot D_\text{temp}$
   - $q \leftarrow D_\text{temp} \cdot q$

   C: `OSQPMatrix_lmult_diag / rmult_diag` and `OSQPVectorf_ew_prod`  
   MATLAB: `Dsp * Pw * Dsp`, `Esp * Aw * Dsp`, `D_temp .* qw`

5. **Accumulate**: $D \leftarrow D \cdot D_\text{temp}$, $E \leftarrow E \cdot E_\text{temp}$

6. **Cost normalization**:
   - Compute $\bar{c} = \frac{1}{n}\sum_j \|P_{\cdot j}\|_\infty$ (mean column inf-norm)
   - $c_\text{temp} = \max(\bar{c}, \|q\|_\infty)$
   - Apply `limit_scaling_scalar`: if $< \text{MIN\_SCALING}$, set to 1; cap at `MAX\_SCALING`
   - $c_\text{temp} \leftarrow 1/c_\text{temp}$
   - $P \leftarrow c_\text{temp} \cdot P$, $q \leftarrow c_\text{temp} \cdot q$
   - $c \leftarrow c \cdot c_\text{temp}$

   C: `OSQPVectorf_norm_1(D_temp) / n` for mean column norm  
   MATLAB: `sum(col_norms_P, 2) / max(n, 1)`

After all iterations: $D^{-1} = 1./D$, $E^{-1} = 1./E$, $c^{-1} = 1/c$. Bounds scaled: $l_s = E \cdot l$, $u_s = E \cdot u$ (clamped to $\pm\text{INFTY}$).

### 3.3. Constraint Classification

C: `ew_bounds_type()` → `set_rho_vec()` (`auxil.c:73–107`)  
MATLAB: `Solver.classifyConstraints()` + `Solver.makeRhoVec()` (lines 722–755)

For each constraint $i$:

| Condition | Type | $\rho_i$ |
|---|---|---|
| $l_i \le -\text{INFTY\_THRESH}$ and $u_i \ge \text{INFTY\_THRESH}$ | Loose ($-1$) | `OSQP_RHO_MIN` |
| $|l_i - u_i| < \text{RHO\_TOL}$ | Equality ($1$) | `OSQP_RHO_EQ_OVER_RHO_INEQ × rho` |
| Otherwise | Inequality ($0$) | `rho` |

where `INFTY_THRESH = OSQP_INFTY × MIN_SCALING = 1e26`.

### 3.4. KKT Factorization

The system solved at each ADMM step is:

$$K = \begin{bmatrix} P + \sigma I & A^T \\ A & -\text{diag}(\rho^{-1}) \end{bmatrix}$$

C: `osqp_algebra_init_linsys_solver()` (various backends)  
MATLAB: `Solver.factorizeKKT()` (lines 682–720), via `osqp.LinearSolver.buildKKT()` and one of:
- `osqp.linsys.MatlabLDLSolver` — MATLAB built-in `ldl()`
- `osqp.linsys.QDLDLSolver` — Pure-MATLAB QDLDL
- `osqp.linsys.QDLDLCSolver` — C MEX QDLDL

### 3.5. Cold Start

C: `osqp_cold_start()` (`osqp_api.c:1299–1306`)  
MATLAB: `Solver.cold_start()` (lines 615–619)

Sets $x = 0$, $z = 0$, $y = 0$.

---

## 4. ADMM Iteration

C: `osqp_solve()` (`osqp_api.c:637–1038`)  
MATLAB: `Solver.solve()` (lines 174–442)

### 4.0. Warm Start Scaling

If `warm_starting` is true, scale the stored (unscaled) iterates into scaled space:

| Variable | C | MATLAB |
|---|---|---|
| $x_s = D^{-1} x$ | `Dinv .* x` | `scl.Dinv .* x_` |
| $z_s = E \cdot z$ | `E .* z` | `scl.E .* z_` |
| $y_s = c \cdot E^{-1} y$ | `c * Einv .* y` | `scl.c * scl.Einv .* y_` |

C: Scaling done at warm_start time in `osqp_warm_start()`.  
MATLAB: Scaling done at top of `solve()` (lines 197–202).

### 4.1. ADMM Step 1: KKT Solve

**Right-hand side**:

$$\text{rhs}_x = \sigma x^k - q_s, \qquad \text{rhs}_z = z^k - \rho^{-1} y^k$$

C: `compute_rhs()` (`auxil.c:131–156`) → `update_xz_tilde()` (`auxil.c:158–169`)  
MATLAB: Lines 233–243

**Solve** $K \cdot [\tilde{x}; \nu] = [\text{rhs}_x; \text{rhs}_z]$ where $\nu$ is the auxiliary variable.

**Extract $\tilde{z}$**: $\tilde{z} = z^k + \rho^{-1}(\nu - y^k)$

C: $\tilde{z}$ is part of `xz_tilde` after the KKT solve  
MATLAB: `ztilde = zs_prev + rho_inv_vec .* (sol(n+1:end) - ys)` (line 245)

### 4.2. ADMM Step 2: Update x (with relaxation)

$$x^{k+1} = \alpha \tilde{x} + (1 - \alpha) x^k$$

C: `update_x()` (`auxil.c:171–183`)  
MATLAB: `xs = s.alpha * xtilde + (1 - s.alpha) * xs_prev` (line 249)

### 4.3. ADMM Step 3: Update z (projection)

$$z^{k+1} = \Pi_{[l_s, u_s]}\bigl(\alpha \tilde{z} + (1-\alpha) z^k + \rho^{-1} y^k\bigr)$$

C: `update_z()` (`auxil.c:185–211`)  
MATLAB: `zs = min(max(zs_relaxed + rho_inv_vec .* ys, ls), us)` (lines 251–254)

### 4.4. ADMM Step 4: Update y (dual)

$$y^{k+1} = y^k + \rho \bigl(\alpha \tilde{z} + (1-\alpha) z^k - z^{k+1}\bigr)$$

C: `update_y()` (`auxil.c:213–233`)  
MATLAB: `ys = ys + rho_vec .* (alpha * ztilde + (1 - alpha) * zs_prev - zs)` (lines 257–259)

---

## 5. Convergence Checking

C: `check_termination()` (`auxil.c:808–927`), calls `compute_prim_res/tol`, `compute_dual_res/tol`, `is_primal_infeasible`, `is_dual_infeasible`  
MATLAB: `Solver.checkConvergence()` (lines 856–1071)

### 5.0. Unscaling Decision

Both C and MATLAB branch on `scaling > 0 && !scaled_termination`:

```
do_unscale = (scaling > 0) && (!scaled_termination)
```

When `scaled_termination = false` (default), residuals and tolerances are computed in the **unscaled** (original) space, which requires applying the inverse scaling factors. When `scaled_termination = true`, all comparisons are in the scaled space.

### 5.1. Primal Residual

$$r_\text{prim} = A_s x_s - z_s$$

| Mode | C (`compute_prim_res`, `auxil.c:328–348`) | MATLAB (lines 889–905) |
|---|---|---|
| Scaled | $\|r_\text{prim}\|_\infty$ | `norm(Ax - z, inf)` |
| Unscaled | $\|E^{-1} \cdot r_\text{prim}\|_\infty$ | `norm(Einv .* (Ax - z), inf)` |

### 5.2. Primal Tolerance

$$\varepsilon_\text{prim} = \varepsilon_\text{abs} + \varepsilon_\text{rel} \cdot \max\bigl(\|z_s\|_\ast, \|A_s x_s\|_\ast\bigr)$$

where $\|\cdot\|_\ast$ is either the plain $\infty$-norm (scaled) or $\|E^{-1} \cdot (\cdot)\|_\infty$ (unscaled).

C: `compute_prim_tol()` (`auxil.c:350–382`)  
MATLAB: Lines 908–920

### 5.3. Dual Residual

$$r_\text{dual} = P_s x_s + q_s + A_s^T y_s$$

| Mode | C (`compute_dual_res`, `auxil.c:384–413`) | MATLAB (lines 907–919) |
|---|---|---|
| Scaled | $\|r_\text{dual}\|_\infty$ | `norm(Px + q + A'y, inf)` |
| Unscaled | $c^{-1} \cdot \|D^{-1} \cdot r_\text{dual}\|_\infty$ | `cinv * norm(Dinv .* (Px + q + A'y), inf)` |

### 5.4. Dual Tolerance

$$\varepsilon_\text{dual} = \varepsilon_\text{abs} + \varepsilon_\text{rel} \cdot \max\bigl(\|q_s\|_\ast, \|A_s^T y_s\|_\ast, \|P_s x_s\|_\ast\bigr)$$

where in unscaled mode: $\|\cdot\|_\ast = c^{-1} \|D^{-1} \cdot (\cdot)\|_\infty$.

C: `compute_dual_tol()` (`auxil.c:415–456`)  
MATLAB: Lines 922–941

### 5.5. Optimality Check

If $r_\text{prim} \le \varepsilon_\text{prim}$ **and** $r_\text{dual} \le \varepsilon_\text{dual}$, declare **SOLVED**.

For approximate check (`approximate = true`), all tolerances are multiplied by 10, and the status is **SOLVED_INACCURATE**.

### 5.6. Primal Infeasibility Certificate ($\delta y$)

C: `is_primal_infeasible()` (`auxil.c:458–510`)  
MATLAB: Lines 972–1027

1. Compute $\delta y = y^{k+1} - y^{k}$
2. **Project onto polar recession cone**: `project_polar_reccone(δy, l, u, INFTY_THRESH)`
   - For each $i$: if both $|l_i|$ and $|u_i| \ge \text{INFTY\_THRESH}$ (loose): $\delta y_i = 0$
   - If only $l_i$ is infinite: $\delta y_i = \min(\delta y_i, 0)$
   - If only $u_i$ is infinite: $\delta y_i = \max(\delta y_i, 0)$

   C: `OSQPVectorf_project_polar_reccone()`  
   MATLAB: Explicit loop (lines 979–991)

3. **Compute norm** (unscale if needed):
   - Scaled: $\|\delta y\|_\infty$
   - Unscaled: $\|E \cdot \delta y\|_\infty$

4. If $\|\delta y\| > \varepsilon_\text{prim\_inf}$:
   - **Support function check**: $u^T \max(\delta y, 0) + l^T \min(\delta y, 0) < 0$ (only for finite bounds)
   - **A' check**: $\|A_s^T \delta y\|_\infty < \varepsilon_\text{prim\_inf} \cdot \|\delta y\|$
     - In unscaled mode: $A_s^T \delta y$ is post-multiplied by $D^{-1}$

### 5.7. Dual Infeasibility Certificate ($\delta x$)

C: `is_dual_infeasible()` (`auxil.c:512–573`)  
MATLAB: Lines 1030–1071

1. Compute $\delta x = x^{k+1} - x^{k}$
2. **Compute norm**:
   - Scaled: $\|\delta x\|_\infty$
   - Unscaled: $\|D \cdot \delta x\|_\infty$
3. **Cost scaling**: $c_s = c$ (unscaled) or $c_s = 1$ (scaled)
4. If $\|\delta x\| > \text{DIVISION\_TOL}$:
   - **Cost check**: $q_s^T \delta x < 0$
   - **Hessian check**: $\|P_s \delta x\|_\infty < c_s \cdot \varepsilon_\text{dual\_inf} \cdot \|\delta x\|$
     - In unscaled mode: $P_s \delta x$ is post-multiplied by $D^{-1}$
   - **Cone check**: $A_s \delta x \in \text{reccone}([l_s, u_s])$
     - In unscaled mode: $A_s \delta x$ is post-multiplied by $E^{-1}$
     - C: `in_reccone(Adx, l, u, INFTY_THRESH, eps * norm_dx)`
     - MATLAB: `dualInfCheck(Adx, ls, us, eps_di * norm_dx)` (lines 1086–1112)

     Recession cone check for each $i$:
     - If both $l_i$ and $u_i$ finite (equality-like): $|(A\delta x)_i| \le \epsilon$
     - If only $l_i$ finite: $(A\delta x)_i \ge -\epsilon$
     - If only $u_i$ finite: $(A\delta x)_i \le \epsilon$
     - If both infinite: no constraint on $(A\delta x)_i$

---

## 6. Time Limit

C: `osqp_solve()` (`osqp_api.c:757–767`) — checked via `can_check_time_limit` flag  
MATLAB: Lines 264–278

When time limit is reached, an approximate convergence check (10× tolerances) is performed first. If that succeeds, the approximate status is used; otherwise `TIME_LIMIT_REACHED`.

---

## 7. Adaptive Rho

### 7.1. Interval

C: `osqp_setup()` (`osqp_api.c:564–604`)  
MATLAB: Lines 280–289

For iteration-based mode (default): if `adaptive_rho_interval == 0` and `check_termination > 0`:

$$\text{interval} = \text{ADAPTIVE\_RHO\_MULTIPLE\_TERMINATION} \times \text{check\_termination} = 4 \times 25 = 100$$

### 7.2. Rho Estimate

C: `compute_rho_estimate()` (`auxil.c:14–48`)  
MATLAB: `Solver.computeNewRho()` (lines 1114–1163)

$$\hat{\rho} = \rho \sqrt{\frac{\|r_\text{prim}\|_\infty / (\max(\|z\|, \|Ax\|) + \delta)}{\|r_\text{dual}\|_\infty / (\max(\|q\|, \|A^Ty\|, \|Px\|) + \delta)}}$$

where $\delta = \text{DIVISION\_TOL}$. Clamped to $[\text{RHO\_MIN}, \text{RHO\_MAX}]$.

**Note**: Both C and MATLAB use the **scaled** residuals and norms for rho estimation, not the unscaled ones.

### 7.3. Rho Update Decision

C: `adapt_rho()` (`auxil.c:50–71`)  
MATLAB: Lines 290–299

Update rho if $\hat{\rho} > \rho \cdot \text{tol}$ or $\hat{\rho} < \rho / \text{tol}$, where `tol = adaptive_rho_tolerance` (default 5.0).

---

## 8. Post-Loop

C: `osqp_solve()` (`osqp_api.c:850–880`)  
MATLAB: Lines 303–324

1. If `status == UNSOLVED`: run exact check, then approximate (10×), then declare `MAX_ITER_REACHED`
2. Store `rho_estimate`

---

## 9. Solution Unscaling

C: `store_solution()` → `unscale_solution()` (`scaling.c`)  
MATLAB: Lines 326–330

$$x_\text{out} = D \cdot x_s, \qquad z_\text{out} = E^{-1} \cdot z_s, \qquad y_\text{out} = c^{-1} \cdot E \cdot y_s$$

### 9.1. Objective Value

Computed on unscaled solution: $\text{obj} = \frac{1}{2} x^T P x + q^T x$

### 9.2. Infeasibility Certificates

C: `store_solution()` normalizes certificates by their $\infty$-norm
MATLAB: Same (lines 365–379, 381–393)

- **Primal infeasibility**: $\delta y_\text{cert} = E \cdot (y_s^{k+1} - y_s^k)$, normalized to unit $\infty$-norm
- **Dual infeasibility**: $\delta x_\text{cert} = D \cdot (x_s^{k+1} - x_s^k)$, normalized to unit $\infty$-norm

---

## 10. Solution Polishing

C: `polish()` (`polish.c:283–515`)  
MATLAB: `Solver.polishSolution()` (lines 1165–1274)

Polishing is performed only when `status == SOLVED` and `polishing == true`.

### 10.1. Active Set Detection

For each constraint $j$, determine active bound:
- **Lower-active**: $y_j < -\delta$ and $|l_j| < \text{INFTY}$, or $|A_j x - l_j| \le \text{tol\_act}$
- **Upper-active**: $y_j > \delta$ and $|u_j| < \text{INFTY}$, or $|A_j x - u_j| \le \text{tol\_act}$

where $\text{tol\_act} = \max(\delta, \varepsilon_\text{abs})$.

### 10.2. Reduced KKT System

$$\begin{bmatrix} P + \delta I & A_\text{act}^T \\ A_\text{act} & -\delta I \end{bmatrix} \begin{bmatrix} x_\text{pol} \\ y_\text{act} \end{bmatrix} = \begin{bmatrix} -q \\ b_\text{act} \end{bmatrix}$$

where $b_{\text{act},j} = l_j$ if lower-active, $u_j$ if upper-active, $(l_j + u_j)/2$ if both.

### 10.3. Iterative Refinement

For `polish_refine_iter` iterations (default 3):

$$\text{res} = \text{rhs} - K \cdot \text{sol}, \qquad \Delta\text{sol} = K^{-1} \text{res}, \qquad \text{sol} \mathrel{+}= \Delta\text{sol}$$

### 10.4. Feasibility Projection (Moreau Decomposition)

After polishing, ensure $z \in [l, u]$ and $y \perp (z - l)(u - z)$:

$$y \leftarrow y + Ax, \qquad z \leftarrow \Pi_{[l,u]}(y), \qquad y \leftarrow y - z$$

### 10.5. Accept/Reject

Polished solution is accepted if $\max(r_\text{prim}^\text{pol}, r_\text{dual}^\text{pol}) \le 10 \cdot \max(r_\text{prim}, r_\text{dual}) + 10^{-10}$.

---

## 11. Warm Start

C: `osqp_warm_start()` (`osqp_api.c:1265–1296`)  
MATLAB: `Solver.warm_start()` (lines 541–570)

- Stores unscaled `x` and `y`
- Sets $z = Ax$ (matches C)
- Scaling to internal space happens at the top of `solve()`

**Difference**: In C, scaling of $(x, y)$ to internal variables happens in `osqp_warm_start()`. In MATLAB, scaling happens at the start of `solve()`. The net effect is identical.

---

## 12. Update Data

### 12.1. Vector Updates (q, l, u)

C: `osqp_update_data_vec()` (`osqp_api.c:1171–1261`)  
MATLAB: `Solver.update()` (lines 414–536)

- `q`: stored unscaled; rescaled after full `scaleProblem()` or quick `qs = c * D .* q`
- `l/u`: clamped to $[-\text{INFTY}, +\text{INFTY}]$, reclassify constraints, rebuild rho

### 12.2. Matrix Updates (Px, Ax)

C: `osqp_update_data_mat()` (`osqp_api.c:1311–1407`) — unscales, updates values, rescales  
MATLAB: Updates `P_triu` or `A_` in-place (sparse values), then rescales via full `scaleProblem()` + `factorizeKKT()`

---

## 13. Known Design Differences

| Aspect | C | MATLAB | Impact |
|---|---|---|---|
| `DIVISION_TOL` | `1/OSQP_INFTY` = $10^{-30}$ | $10^{-10}$ | Negligible — both prevent division by zero |
| Duality gap check | Implemented via `check_dualgap` | Not implemented (always `duality_gap_check = 1`) | `check_dualgap` defaults to `true` in C; MATLAB doesn't compute duality gap as a termination criterion |
| Pointer swaps | Uses in-place pointer swaps for `x_prev`, `z_prev` | Copies via `xs_prev = xs` | Same result, MATLAB copy semantics |
| Adaptive rho modes | Supports `ITERATIONS`, `TIME`, `KKT_ERROR`, `DISABLED` | Supports `ITERATIONS` mode (boolean `adaptive_rho`) | Time/KKT modes not applicable to MATLAB |
| Linear system solve | Backend-specific (QDLDL C, MKL, CUDA, CG) | MATLAB LDL, QDLDL MATLAB, QDLDL C MEX | Same KKT structure |

---

## 14. Status Code Mapping

| Status | C Enum | Value | MATLAB |
|---|---|---|---|
| Solved | `OSQP_SOLVED` | 1 | `osqp.constant.OSQP_SOLVED` |
| Solved (inaccurate) | `OSQP_SOLVED_INACCURATE` | 2 | `osqp.constant.OSQP_SOLVED_INACCURATE` |
| Primal infeasible | `OSQP_PRIMAL_INFEASIBLE` | 3 | `osqp.constant.OSQP_PRIMAL_INFEASIBLE` |
| Primal infeasible (inaccurate) | `OSQP_PRIMAL_INFEASIBLE_INACCURATE` | 4 | `osqp.constant.OSQP_PRIMAL_INFEASIBLE_INACCURATE` |
| Dual infeasible | `OSQP_DUAL_INFEASIBLE` | 5 | `osqp.constant.OSQP_DUAL_INFEASIBLE` |
| Dual infeasible (inaccurate) | `OSQP_DUAL_INFEASIBLE_INACCURATE` | 6 | `osqp.constant.OSQP_DUAL_INFEASIBLE_INACCURATE` |
| Max iterations | `OSQP_MAX_ITER_REACHED` | 7 | `osqp.constant.OSQP_MAX_ITER_REACHED` |
| Time limit | `OSQP_TIME_LIMIT_REACHED` | 8 | `osqp.constant.OSQP_TIME_LIMIT_REACHED` |
| Non-convex | `OSQP_NON_CVX` | 9 | `osqp.constant.OSQP_NON_CONVEX` |
| Interrupted | `OSQP_SIGINT` | 10 | `osqp.constant.OSQP_SIGINT` |
| Unsolved | `OSQP_UNSOLVED` | 11 | `osqp.constant.OSQP_UNSOLVED` |
