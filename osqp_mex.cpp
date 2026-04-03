/*
 * osqp_mex.cpp — MATLAB MEX interface for OSQP v1.0.0
 *
 * Commands:
 *   'new'                         → create handle
 *   'delete'                      → destroy handle
 *   'setup'                       → osqp_setup
 *   'solve'                       → osqp_solve
 *   'update'                      → osqp_update_data_vec / osqp_update_data_mat
 *   'warm_start'                  → osqp_warm_start
 *   'cold_start'                  → osqp_cold_start
 *   'update_settings'             → osqp_update_settings
 *   'update_rho'                  → osqp_update_rho
 *   'default_settings'            → osqp_set_default_settings
 *   'current_settings'            → read solver->settings
 *   'get_dimensions'              → osqp_get_dimensions
 *   'version'                     → osqp_version
 *   'constant'                    → return compile-time constants
 *   'capabilities'                → osqp_capabilities
 *   'codegen'                     → osqp_codegen
 *   'adjoint_derivative_compute'  → osqp_adjoint_derivative_compute
 *   'adjoint_derivative_get_mat'  → osqp_adjoint_derivative_get_mat
 *   'adjoint_derivative_get_vec'  → osqp_adjoint_derivative_get_vec
 */

#include "mex.h"
#include "matrix.h"
#include "osqp_mex.hpp"
#include "osqp.h"

#include <cstring>
#include <cstdlib>
#include <cmath>

/* ------------------------------------------------------------------ */
/* Field-name arrays for MATLAB structs                                */
/* ------------------------------------------------------------------ */

static const char* OSQP_INFO_FIELDS[] = {
    "iter",           // OSQPInt
    "status",         // char[]
    "status_val",     // OSQPInt
    "status_polish",  // OSQPInt
    "obj_val",        // OSQPFloat
    "dual_obj_val",   // OSQPFloat   (new in v1.0.0)
    "prim_res",       // OSQPFloat
    "dual_res",       // OSQPFloat
    "duality_gap",    // OSQPFloat   (new)
    "setup_time",     // OSQPFloat
    "solve_time",     // OSQPFloat
    "update_time",    // OSQPFloat
    "polish_time",    // OSQPFloat
    "run_time",       // OSQPFloat
    "rho_updates",    // OSQPInt
    "rho_estimate",   // OSQPFloat
    "primdual_int",   // OSQPFloat   (new)
    "rel_kkt_error"   // OSQPFloat   (new)
};

static const char* OSQP_SETTINGS_FIELDS[] = {
    "device",
    "linsys_solver",
    "allocate_solution",
    "verbose",
    "profiler_level",
    "warm_starting",
    "scaling",
    "polishing",
    "rho",
    "rho_is_vec",
    "sigma",
    "alpha",
    "cg_max_iter",
    "cg_tol_reduction",
    "cg_tol_fraction",
    "cg_precond",
    "adaptive_rho",
    "adaptive_rho_interval",
    "adaptive_rho_fraction",
    "adaptive_rho_tolerance",
    "max_iter",
    "eps_abs",
    "eps_rel",
    "eps_prim_inf",
    "eps_dual_inf",
    "scaled_termination",
    "check_termination",
    "check_dualgap",
    "time_limit",
    "delta",
    "polish_refine_iter"
};

/* ------------------------------------------------------------------ */
/* Wrapper class: holds the OSQPSolver* handle                         */
/* ------------------------------------------------------------------ */
class OsqpData {
public:
    OsqpData() : solver(nullptr) {}
    OSQPSolver* solver;
};

/* ------------------------------------------------------------------ */
/* Forward declarations — utility helpers                              */
/* ------------------------------------------------------------------ */
static void       castToDoubleArr(const OSQPFloat* arr, double* out, OSQPInt len);
static void       setToNaN(double* out, OSQPInt len);
static OSQPFloat* copyToCfloatVector(const double* src, OSQPInt numel);
static OSQPInt*   copyToCintVector(const mwIndex* src, OSQPInt numel);
static OSQPInt*   copyDoubleToCintVector(const double* src, OSQPInt numel);

static mxArray*   copyInfoToMxStruct(const OSQPInfo* info);
static mxArray*   copySettingsToMxStruct(const OSQPSettings* s);
static void       copyMxStructToSettings(const mxArray* mx, OSQPSettings* s);

/* ================================================================== */
/* mexFunction — main entry point                                      */
/* ================================================================== */
void mexFunction(int nlhs, mxArray* plhs[],
                 int nrhs, const mxArray* prhs[])
{
    OsqpData* osqpData = nullptr;
    OSQPInt   exitflag  = 0;
    char      cmd[64];
    char      stat_string[64];

    if (nrhs < 1 || mxGetString(prhs[0], cmd, sizeof(cmd)))
        mexErrMsgTxt("First input should be a command string less than 64 characters long.");

    /* ---- new ---- */
    if (!strcmp("new", cmd)) {
        if (nlhs != 1)
            mexErrMsgTxt("new: one output expected.");
        plhs[0] = convertPtr2Mat<OsqpData>(new OsqpData);
        return;
    }

    /* Second argument: instance handle or 'static' */
    if (nrhs < 2)
        mexErrMsgTxt("Second input should be a class instance handle or the string 'static'.");

    if (mxGetString(prhs[1], stat_string, sizeof(stat_string))) {
        osqpData = convertMat2Ptr<OsqpData>(prhs[1]);
    } else {
        if (strcmp("static", stat_string))
            mexErrMsgTxt("Second argument for static functions must be the string 'static'.");
    }

    /* ---- delete ---- */
    if (!strcmp("delete", cmd)) {
        if (osqpData && osqpData->solver) {
            osqp_cleanup(osqpData->solver);
            osqpData->solver = nullptr;
        }
        destroyObject<OsqpData>(prhs[1]);
        return;
    }

    /* ---- version ---- */
    if (!strcmp("version", cmd)) {
        plhs[0] = mxCreateString(osqp_version());
        return;
    }

    /* ---- capabilities ---- */
    if (!strcmp("capabilities", cmd)) {
        plhs[0] = mxCreateDoubleScalar((double)osqp_capabilities());
        return;
    }

    /* ---- default_settings ---- */
    if (!strcmp("default_settings", cmd)) {
        OSQPSettings* defaults = (OSQPSettings*)mxCalloc(1, sizeof(OSQPSettings));
        osqp_set_default_settings(defaults);
        plhs[0] = copySettingsToMxStruct(defaults);
        mxFree(defaults);
        return;
    }

    /* ---- current_settings ---- */
    if (!strcmp("current_settings", cmd)) {
        if (!osqpData || !osqpData->solver)
            mexErrMsgTxt("Solver is uninitialized.");
        plhs[0] = copySettingsToMxStruct(osqpData->solver->settings);
        return;
    }

    /* ---- constant ---- */
    if (!strcmp("constant", cmd)) {
        char constant[64];
        mxGetString(prhs[2], constant, sizeof(constant));

        if (!strcmp("OSQP_INFTY", constant))               { plhs[0] = mxCreateDoubleScalar(OSQP_INFTY); return; }
        if (!strcmp("OSQP_NAN", constant))                  { plhs[0] = mxCreateDoubleScalar(mxGetNaN()); return; }
        if (!strcmp("OSQP_SOLVED", constant))               { plhs[0] = mxCreateDoubleScalar(OSQP_SOLVED); return; }
        if (!strcmp("OSQP_SOLVED_INACCURATE", constant))    { plhs[0] = mxCreateDoubleScalar(OSQP_SOLVED_INACCURATE); return; }
        if (!strcmp("OSQP_UNSOLVED", constant))             { plhs[0] = mxCreateDoubleScalar(OSQP_UNSOLVED); return; }
        if (!strcmp("OSQP_PRIMAL_INFEASIBLE", constant))    { plhs[0] = mxCreateDoubleScalar(OSQP_PRIMAL_INFEASIBLE); return; }
        if (!strcmp("OSQP_PRIMAL_INFEASIBLE_INACCURATE", constant)) { plhs[0] = mxCreateDoubleScalar(OSQP_PRIMAL_INFEASIBLE_INACCURATE); return; }
        if (!strcmp("OSQP_DUAL_INFEASIBLE", constant))      { plhs[0] = mxCreateDoubleScalar(OSQP_DUAL_INFEASIBLE); return; }
        if (!strcmp("OSQP_DUAL_INFEASIBLE_INACCURATE", constant)) { plhs[0] = mxCreateDoubleScalar(OSQP_DUAL_INFEASIBLE_INACCURATE); return; }
        if (!strcmp("OSQP_MAX_ITER_REACHED", constant))     { plhs[0] = mxCreateDoubleScalar(OSQP_MAX_ITER_REACHED); return; }
        if (!strcmp("OSQP_TIME_LIMIT_REACHED", constant))   { plhs[0] = mxCreateDoubleScalar(OSQP_TIME_LIMIT_REACHED); return; }
        if (!strcmp("OSQP_NON_CVX", constant))              { plhs[0] = mxCreateDoubleScalar(OSQP_NON_CVX); return; }
        if (!strcmp("OSQP_SIGINT", constant))               { plhs[0] = mxCreateDoubleScalar(OSQP_SIGINT); return; }

        /* Solver types */
        if (!strcmp("OSQP_DIRECT_SOLVER", constant))        { plhs[0] = mxCreateDoubleScalar(OSQP_DIRECT_SOLVER); return; }
        if (!strcmp("OSQP_INDIRECT_SOLVER", constant))      { plhs[0] = mxCreateDoubleScalar(OSQP_INDIRECT_SOLVER); return; }

        /* Preconditioners */
        if (!strcmp("OSQP_NO_PRECONDITIONER", constant))    { plhs[0] = mxCreateDoubleScalar(OSQP_NO_PRECONDITIONER); return; }
        if (!strcmp("OSQP_DIAGONAL_PRECONDITIONER", constant)) { plhs[0] = mxCreateDoubleScalar(OSQP_DIAGONAL_PRECONDITIONER); return; }

        /* Capabilities */
        if (!strcmp("OSQP_CAPABILITY_DIRECT_SOLVER", constant))   { plhs[0] = mxCreateDoubleScalar(OSQP_CAPABILITY_DIRECT_SOLVER); return; }
        if (!strcmp("OSQP_CAPABILITY_INDIRECT_SOLVER", constant)) { plhs[0] = mxCreateDoubleScalar(OSQP_CAPABILITY_INDIRECT_SOLVER); return; }
        if (!strcmp("OSQP_CAPABILITY_CODEGEN", constant))         { plhs[0] = mxCreateDoubleScalar(OSQP_CAPABILITY_CODEGEN); return; }
        if (!strcmp("OSQP_CAPABILITY_UPDATE_MATRICES", constant)) { plhs[0] = mxCreateDoubleScalar(OSQP_CAPABILITY_UPDATE_MATRICES); return; }
        if (!strcmp("OSQP_CAPABILITY_DERIVATIVES", constant))     { plhs[0] = mxCreateDoubleScalar(OSQP_CAPABILITY_DERIVATIVES); return; }

        /* Polish status */
        if (!strcmp("OSQP_POLISH_SUCCESS", constant))              { plhs[0] = mxCreateDoubleScalar(OSQP_POLISH_SUCCESS); return; }
        if (!strcmp("OSQP_POLISH_FAILED", constant))               { plhs[0] = mxCreateDoubleScalar(OSQP_POLISH_FAILED); return; }
        if (!strcmp("OSQP_POLISH_NOT_PERFORMED", constant))        { plhs[0] = mxCreateDoubleScalar(OSQP_POLISH_NOT_PERFORMED); return; }
        if (!strcmp("OSQP_POLISH_LINSYS_ERROR", constant))         { plhs[0] = mxCreateDoubleScalar(OSQP_POLISH_LINSYS_ERROR); return; }
        if (!strcmp("OSQP_POLISH_NO_ACTIVE_SET_FOUND", constant))  { plhs[0] = mxCreateDoubleScalar(OSQP_POLISH_NO_ACTIVE_SET_FOUND); return; }

        mexErrMsgTxt("Constant not recognized.");
        return;
    }

    /* ---- setup ---- */
    if (!strcmp("setup", cmd)) {
        if (osqpData->solver)
            mexErrMsgTxt("Solver is already initialized.");

        /* args: cmd, handle, n, m, P, q, A, l, u, settings */
        OSQPInt n = (OSQPInt)mxGetScalar(prhs[2]);
        OSQPInt m = (OSQPInt)mxGetScalar(prhs[3]);

        const mxArray* mxP = prhs[4];
        const mxArray* mxq = prhs[5];
        const mxArray* mxA = prhs[6];
        const mxArray* mxl = prhs[7];
        const mxArray* mxu = prhs[8];

        /* Vectors */
        OSQPFloat* q_vec = copyToCfloatVector(mxGetPr(mxq), n);
        OSQPFloat* l_vec = copyToCfloatVector(mxGetPr(mxl), m);
        OSQPFloat* u_vec = copyToCfloatVector(mxGetPr(mxu), m);

        /* Matrix P (upper triangular, CSC) */
        OSQPInt*   Pp = copyToCintVector(mxGetJc(mxP), n + 1);
        OSQPInt    Pnnz = Pp[n];
        OSQPInt*   Pi = copyToCintVector(mxGetIr(mxP), Pnnz);
        OSQPFloat* Px = copyToCfloatVector(mxGetPr(mxP), Pnnz);

        OSQPCscMatrix P_csc;
        OSQPCscMatrix_set_data(&P_csc, n, n, Pnnz, Px, Pi, Pp);

        /* Matrix A (CSC) */
        OSQPInt*   Ap = copyToCintVector(mxGetJc(mxA), n + 1);
        OSQPInt    Annz = Ap[n];
        OSQPInt*   Ai = copyToCintVector(mxGetIr(mxA), Annz);
        OSQPFloat* Ax = copyToCfloatVector(mxGetPr(mxA), Annz);

        OSQPCscMatrix A_csc;
        OSQPCscMatrix_set_data(&A_csc, m, n, Annz, Ax, Ai, Ap);

        /* Settings */
        OSQPSettings* settings = (OSQPSettings*)mxCalloc(1, sizeof(OSQPSettings));
        osqp_set_default_settings(settings);
        const mxArray* mxSettings = prhs[9];
        if (!mxIsEmpty(mxSettings))
            copyMxStructToSettings(mxSettings, settings);

        /* Setup solver */
        exitflag = osqp_setup(&(osqpData->solver),
                              &P_csc, q_vec, &A_csc, l_vec, u_vec,
                              m, n, settings);

        /* Cleanup temporaries */
        mxFree(settings);
        free(q_vec); free(l_vec); free(u_vec);
        free(Pp); free(Pi); free(Px);
        free(Ap); free(Ai); free(Ax);

        if (exitflag)
            mexErrMsgTxt("Invalid problem setup.");
        return;
    }

    /* ---- get_dimensions ---- */
    if (!strcmp("get_dimensions", cmd)) {
        if (!osqpData || !osqpData->solver)
            mexErrMsgTxt("Solver is uninitialized.");
        OSQPInt n_out, m_out;
        osqp_get_dimensions(osqpData->solver, &m_out, &n_out);
        plhs[0] = mxCreateDoubleScalar((double)n_out);
        plhs[1] = mxCreateDoubleScalar((double)m_out);
        return;
    }

    /* ---- solve ---- */
    if (!strcmp("solve", cmd)) {
        if (nlhs != 5 || nrhs != 2)
            mexErrMsgTxt("solve: expects 0 extra inputs and 5 outputs.");
        if (!osqpData || !osqpData->solver)
            mexErrMsgTxt("Solver is uninitialized.");

        OSQPSolver* solver = osqpData->solver;
        osqp_solve(solver);

        OSQPInt n, m;
        osqp_get_dimensions(solver, &m, &n);

        /* Allocate outputs */
        plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL); // x
        plhs[1] = mxCreateDoubleMatrix(m, 1, mxREAL); // y
        plhs[2] = mxCreateDoubleMatrix(m, 1, mxREAL); // prim_inf_cert
        plhs[3] = mxCreateDoubleMatrix(n, 1, mxREAL); // dual_inf_cert

        OSQPInt sv = solver->info->status_val;

        if (sv != OSQP_PRIMAL_INFEASIBLE &&
            sv != OSQP_PRIMAL_INFEASIBLE_INACCURATE &&
            sv != OSQP_DUAL_INFEASIBLE &&
            sv != OSQP_DUAL_INFEASIBLE_INACCURATE) {
            /* Normal solve — copy primal/dual solution */
            castToDoubleArr(solver->solution->x, mxGetPr(plhs[0]), n);
            castToDoubleArr(solver->solution->y, mxGetPr(plhs[1]), m);
            setToNaN(mxGetPr(plhs[2]), m);
            setToNaN(mxGetPr(plhs[3]), n);
        } else if (sv == OSQP_PRIMAL_INFEASIBLE ||
                   sv == OSQP_PRIMAL_INFEASIBLE_INACCURATE) {
            setToNaN(mxGetPr(plhs[0]), n);
            setToNaN(mxGetPr(plhs[1]), m);
            castToDoubleArr(solver->solution->prim_inf_cert, mxGetPr(plhs[2]), m);
            setToNaN(mxGetPr(plhs[3]), n);
            solver->info->obj_val = mxGetInf();
        } else {
            /* Dual infeasible */
            setToNaN(mxGetPr(plhs[0]), n);
            setToNaN(mxGetPr(plhs[1]), m);
            setToNaN(mxGetPr(plhs[2]), m);
            castToDoubleArr(solver->solution->dual_inf_cert, mxGetPr(plhs[3]), n);
            solver->info->obj_val = -mxGetInf();
        }

        if (sv == OSQP_NON_CVX)
            solver->info->obj_val = mxGetNaN();

        plhs[4] = copyInfoToMxStruct(solver->info);
        return;
    }

    /* ---- update (data vectors and/or matrices) ---- */
    if (!strcmp("update", cmd)) {
        if (!osqpData || !osqpData->solver)
            mexErrMsgTxt("Solver is uninitialized.");

        /* args: cmd, handle, q, l, u, Px, Px_idx, Px_n, Ax, Ax_idx, Ax_n */
        const mxArray* mxq      = prhs[2];
        const mxArray* mxl      = prhs[3];
        const mxArray* mxu      = prhs[4];
        const mxArray* mxPx     = prhs[5];
        const mxArray* mxPx_idx = prhs[6];
        OSQPInt Px_n = (OSQPInt)mxGetScalar(prhs[7]);
        const mxArray* mxAx     = prhs[8];
        const mxArray* mxAx_idx = prhs[9];
        OSQPInt Ax_n = (OSQPInt)mxGetScalar(prhs[10]);

        OSQPInt n_dim, m_dim;
        osqp_get_dimensions(osqpData->solver, &m_dim, &n_dim);

        /* Vector updates */
        OSQPFloat* q_new = nullptr;
        OSQPFloat* l_new = nullptr;
        OSQPFloat* u_new = nullptr;
        if (!mxIsEmpty(mxq)) q_new = copyToCfloatVector(mxGetPr(mxq), n_dim);
        if (!mxIsEmpty(mxl)) l_new = copyToCfloatVector(mxGetPr(mxl), m_dim);
        if (!mxIsEmpty(mxu)) u_new = copyToCfloatVector(mxGetPr(mxu), m_dim);

        if (q_new || l_new || u_new) {
            exitflag = osqp_update_data_vec(osqpData->solver, q_new, l_new, u_new);
            if (q_new) free(q_new);
            if (l_new) free(l_new);
            if (u_new) free(u_new);
            if (exitflag)
                mexErrMsgTxt("Vector data update error.");
        }

        /* Matrix updates */
        OSQPFloat* Px_vec     = nullptr;
        OSQPInt*   Px_idx_vec = nullptr;
        OSQPFloat* Ax_vec     = nullptr;
        OSQPInt*   Ax_idx_vec = nullptr;

        if (!mxIsEmpty(mxPx)) Px_vec = copyToCfloatVector(mxGetPr(mxPx), Px_n);
        if (!mxIsEmpty(mxAx)) Ax_vec = copyToCfloatVector(mxGetPr(mxAx), Ax_n);
        if (!mxIsEmpty(mxPx_idx)) Px_idx_vec = copyDoubleToCintVector(mxGetPr(mxPx_idx), Px_n);
        if (!mxIsEmpty(mxAx_idx)) Ax_idx_vec = copyDoubleToCintVector(mxGetPr(mxAx_idx), Ax_n);

        if (Px_vec || Ax_vec) {
            exitflag = osqp_update_data_mat(osqpData->solver,
                                            Px_vec, Px_idx_vec, Px_n,
                                            Ax_vec, Ax_idx_vec, Ax_n);
            if (Px_vec)     free(Px_vec);
            if (Ax_vec)     free(Ax_vec);
            if (Px_idx_vec) free(Px_idx_vec);
            if (Ax_idx_vec) free(Ax_idx_vec);
            if (exitflag)
                mexErrMsgTxt("Matrix data update error.");
        }
        return;
    }

    /* ---- warm_start ---- */
    if (!strcmp("warm_start", cmd)) {
        if (!osqpData || !osqpData->solver)
            mexErrMsgTxt("Solver is uninitialized.");

        const mxArray* mxX = prhs[2];
        const mxArray* mxY = prhs[3];

        OSQPInt n_dim, m_dim;
        osqp_get_dimensions(osqpData->solver, &m_dim, &n_dim);

        OSQPFloat* x_vec = nullptr;
        OSQPFloat* y_vec = nullptr;
        if (!mxIsEmpty(mxX)) x_vec = copyToCfloatVector(mxGetPr(mxX), n_dim);
        if (!mxIsEmpty(mxY)) y_vec = copyToCfloatVector(mxGetPr(mxY), m_dim);

        osqp_warm_start(osqpData->solver, x_vec, y_vec);

        if (x_vec) free(x_vec);
        if (y_vec) free(y_vec);
        return;
    }

    /* ---- cold_start ---- */
    if (!strcmp("cold_start", cmd)) {
        if (!osqpData || !osqpData->solver)
            mexErrMsgTxt("Solver is uninitialized.");
        osqp_cold_start(osqpData->solver);
        return;
    }

    /* ---- update_settings ---- */
    if (!strcmp("update_settings", cmd)) {
        if (!osqpData || !osqpData->solver)
            mexErrMsgTxt("Solver is uninitialized.");

        /* Build a full settings struct from current + overrides */
        OSQPSettings new_settings;
        memcpy(&new_settings, osqpData->solver->settings, sizeof(OSQPSettings));
        copyMxStructToSettings(prhs[2], &new_settings);

        exitflag = osqp_update_settings(osqpData->solver, &new_settings);
        if (exitflag)
            mexErrMsgTxt("Settings update error.");
        return;
    }

    /* ---- update_rho ---- */
    if (!strcmp("update_rho", cmd)) {
        if (!osqpData || !osqpData->solver)
            mexErrMsgTxt("Solver is uninitialized.");
        OSQPFloat rho_new = (OSQPFloat)mxGetScalar(prhs[2]);
        exitflag = osqp_update_rho(osqpData->solver, rho_new);
        if (exitflag)
            mexErrMsgTxt("Rho update error.");
        return;
    }

    /* ---- codegen ---- */
    if (!strcmp("codegen", cmd)) {
        if (!osqpData || !osqpData->solver)
            mexErrMsgTxt("Solver is uninitialized.");

        /* args: cmd, handle, output_dir, prefix, defines_struct */
        char output_dir[512];
        char prefix[128];
        mxGetString(prhs[2], output_dir, sizeof(output_dir));
        mxGetString(prhs[3], prefix, sizeof(prefix));

        const mxArray* mxDefines = prhs[4];
        OSQPCodegenDefines* defines = OSQPCodegenDefines_new();
        osqp_set_default_codegen_defines(defines);

        if (!mxIsEmpty(mxDefines)) {
            mxArray* f;
            if ((f = mxGetField(mxDefines, 0, "embedded_mode")))
                defines->embedded_mode = (OSQPInt)mxGetScalar(f);
            if ((f = mxGetField(mxDefines, 0, "float_type")))
                defines->float_type = (OSQPInt)mxGetScalar(f);
            if ((f = mxGetField(mxDefines, 0, "printing_enable")))
                defines->printing_enable = (OSQPInt)mxGetScalar(f);
            if ((f = mxGetField(mxDefines, 0, "profiling_enable")))
                defines->profiling_enable = (OSQPInt)mxGetScalar(f);
            if ((f = mxGetField(mxDefines, 0, "interrupt_enable")))
                defines->interrupt_enable = (OSQPInt)mxGetScalar(f);
            if ((f = mxGetField(mxDefines, 0, "derivatives_enable")))
                defines->derivatives_enable = (OSQPInt)mxGetScalar(f);
        }

        exitflag = osqp_codegen(osqpData->solver, output_dir, prefix, defines);
        OSQPCodegenDefines_free(defines);

        if (exitflag)
            mexErrMsgTxt("Code generation failed.");
        return;
    }

    /* ---- adjoint_derivative_compute ---- */
    if (!strcmp("adjoint_derivative_compute", cmd)) {
        if (!osqpData || !osqpData->solver)
            mexErrMsgTxt("Solver is uninitialized.");

        OSQPInt n_dim, m_dim;
        osqp_get_dimensions(osqpData->solver, &m_dim, &n_dim);

        OSQPFloat* dx = copyToCfloatVector(mxGetPr(prhs[2]), n_dim);
        OSQPFloat* dy = copyToCfloatVector(mxGetPr(prhs[3]), m_dim);

        exitflag = osqp_adjoint_derivative_compute(osqpData->solver, dx, dy);

        free(dx);
        free(dy);

        if (exitflag)
            mexErrMsgTxt("Adjoint derivative computation failed.");
        return;
    }

    /* ---- adjoint_derivative_get_mat ---- */
    if (!strcmp("adjoint_derivative_get_mat", cmd)) {
        if (!osqpData || !osqpData->solver)
            mexErrMsgTxt("Solver is uninitialized.");

        OSQPInt n_dim, m_dim;
        osqp_get_dimensions(osqpData->solver, &m_dim, &n_dim);

        /* Allocate output CSC matrices */
        OSQPCscMatrix* dP = OSQPCscMatrix_zeros(n_dim, n_dim);
        OSQPCscMatrix* dA = OSQPCscMatrix_zeros(m_dim, n_dim);

        exitflag = osqp_adjoint_derivative_get_mat(osqpData->solver, dP, dA);

        if (exitflag) {
            OSQPCscMatrix_free(dP);
            OSQPCscMatrix_free(dA);
            mexErrMsgTxt("Adjoint derivative matrix retrieval failed.");
        }

        /* Convert dP to MATLAB sparse */
        OSQPInt dP_nnz = dP->p[n_dim];
        mxArray* mxdP = mxCreateSparse(n_dim, n_dim, dP_nnz > 0 ? dP_nnz : 1, mxREAL);
        for (OSQPInt j = 0; j <= n_dim; j++) ((mwIndex*)mxGetJc(mxdP))[j] = (mwIndex)dP->p[j];
        for (OSQPInt k = 0; k < dP_nnz; k++) {
            ((mwIndex*)mxGetIr(mxdP))[k] = (mwIndex)dP->i[k];
            mxGetPr(mxdP)[k] = (double)dP->x[k];
        }
        plhs[0] = mxdP;

        /* Convert dA to MATLAB sparse */
        OSQPInt dA_nnz = dA->p[n_dim];
        mxArray* mxdA = mxCreateSparse(m_dim, n_dim, dA_nnz > 0 ? dA_nnz : 1, mxREAL);
        for (OSQPInt j = 0; j <= n_dim; j++) ((mwIndex*)mxGetJc(mxdA))[j] = (mwIndex)dA->p[j];
        for (OSQPInt k = 0; k < dA_nnz; k++) {
            ((mwIndex*)mxGetIr(mxdA))[k] = (mwIndex)dA->i[k];
            mxGetPr(mxdA)[k] = (double)dA->x[k];
        }
        plhs[1] = mxdA;

        OSQPCscMatrix_free(dP);
        OSQPCscMatrix_free(dA);
        return;
    }

    /* ---- adjoint_derivative_get_vec ---- */
    if (!strcmp("adjoint_derivative_get_vec", cmd)) {
        if (!osqpData || !osqpData->solver)
            mexErrMsgTxt("Solver is uninitialized.");

        OSQPInt n_dim, m_dim;
        osqp_get_dimensions(osqpData->solver, &m_dim, &n_dim);

        OSQPFloat* dq = (OSQPFloat*)calloc(n_dim, sizeof(OSQPFloat));
        OSQPFloat* dl = (OSQPFloat*)calloc(m_dim, sizeof(OSQPFloat));
        OSQPFloat* du = (OSQPFloat*)calloc(m_dim, sizeof(OSQPFloat));

        exitflag = osqp_adjoint_derivative_get_vec(osqpData->solver, dq, dl, du);

        if (exitflag) {
            free(dq); free(dl); free(du);
            mexErrMsgTxt("Adjoint derivative vector retrieval failed.");
        }

        plhs[0] = mxCreateDoubleMatrix(n_dim, 1, mxREAL);
        plhs[1] = mxCreateDoubleMatrix(m_dim, 1, mxREAL);
        plhs[2] = mxCreateDoubleMatrix(m_dim, 1, mxREAL);
        castToDoubleArr(dq, mxGetPr(plhs[0]), n_dim);
        castToDoubleArr(dl, mxGetPr(plhs[1]), m_dim);
        castToDoubleArr(du, mxGetPr(plhs[2]), m_dim);

        free(dq); free(dl); free(du);
        return;
    }

    mexErrMsgTxt("Command not recognized.");
}

/* ================================================================== */
/* Utility functions                                                   */
/* ================================================================== */

static OSQPFloat* copyToCfloatVector(const double* src, OSQPInt numel) {
    OSQPFloat* out = (OSQPFloat*)malloc(numel * sizeof(OSQPFloat));
    for (OSQPInt i = 0; i < numel; i++)
        out[i] = (OSQPFloat)src[i];
    return out;
}

static OSQPInt* copyToCintVector(const mwIndex* src, OSQPInt numel) {
    OSQPInt* out = (OSQPInt*)malloc(numel * sizeof(OSQPInt));
    for (OSQPInt i = 0; i < numel; i++)
        out[i] = (OSQPInt)src[i];
    return out;
}

static OSQPInt* copyDoubleToCintVector(const double* src, OSQPInt numel) {
    OSQPInt* out = (OSQPInt*)malloc(numel * sizeof(OSQPInt));
    for (OSQPInt i = 0; i < numel; i++)
        out[i] = (OSQPInt)src[i];
    return out;
}

static void castToDoubleArr(const OSQPFloat* arr, double* out, OSQPInt len) {
    for (OSQPInt i = 0; i < len; i++)
        out[i] = (double)arr[i];
}

static void setToNaN(double* out, OSQPInt len) {
    for (OSQPInt i = 0; i < len; i++)
        out[i] = mxGetNaN();
}

/* ------------------------------------------------------------------ */
/* Info struct → MATLAB struct                                         */
/* ------------------------------------------------------------------ */
static mxArray* copyInfoToMxStruct(const OSQPInfo* info) {
    int nfields = sizeof(OSQP_INFO_FIELDS) / sizeof(OSQP_INFO_FIELDS[0]);
    mxArray* mx = mxCreateStructMatrix(1, 1, nfields, OSQP_INFO_FIELDS);

    mxSetField(mx, 0, "iter",           mxCreateDoubleScalar(info->iter));
    mxSetField(mx, 0, "status",         mxCreateString(info->status));
    mxSetField(mx, 0, "status_val",     mxCreateDoubleScalar(info->status_val));
    mxSetField(mx, 0, "status_polish",  mxCreateDoubleScalar(info->status_polish));
    mxSetField(mx, 0, "obj_val",        mxCreateDoubleScalar(info->obj_val));
    mxSetField(mx, 0, "dual_obj_val",   mxCreateDoubleScalar(info->dual_obj_val));
    mxSetField(mx, 0, "prim_res",       mxCreateDoubleScalar(info->prim_res));
    mxSetField(mx, 0, "dual_res",       mxCreateDoubleScalar(info->dual_res));
    mxSetField(mx, 0, "duality_gap",    mxCreateDoubleScalar(info->duality_gap));
    mxSetField(mx, 0, "setup_time",     mxCreateDoubleScalar(info->setup_time));
    mxSetField(mx, 0, "solve_time",     mxCreateDoubleScalar(info->solve_time));
    mxSetField(mx, 0, "update_time",    mxCreateDoubleScalar(info->update_time));
    mxSetField(mx, 0, "polish_time",    mxCreateDoubleScalar(info->polish_time));
    mxSetField(mx, 0, "run_time",       mxCreateDoubleScalar(info->run_time));
    mxSetField(mx, 0, "rho_updates",    mxCreateDoubleScalar(info->rho_updates));
    mxSetField(mx, 0, "rho_estimate",   mxCreateDoubleScalar(info->rho_estimate));
    mxSetField(mx, 0, "primdual_int",   mxCreateDoubleScalar(info->primdual_int));
    mxSetField(mx, 0, "rel_kkt_error",  mxCreateDoubleScalar(info->rel_kkt_error));

    return mx;
}

/* ------------------------------------------------------------------ */
/* Settings → MATLAB struct                                            */
/* ------------------------------------------------------------------ */
static mxArray* copySettingsToMxStruct(const OSQPSettings* s) {
    int nfields = sizeof(OSQP_SETTINGS_FIELDS) / sizeof(OSQP_SETTINGS_FIELDS[0]);
    mxArray* mx = mxCreateStructMatrix(1, 1, nfields, OSQP_SETTINGS_FIELDS);

    mxSetField(mx, 0, "device",                mxCreateDoubleScalar(s->device));
    mxSetField(mx, 0, "linsys_solver",         mxCreateDoubleScalar(s->linsys_solver));
    mxSetField(mx, 0, "allocate_solution",     mxCreateDoubleScalar(s->allocate_solution));
    mxSetField(mx, 0, "verbose",               mxCreateDoubleScalar(s->verbose));
    mxSetField(mx, 0, "profiler_level",        mxCreateDoubleScalar(s->profiler_level));
    mxSetField(mx, 0, "warm_starting",         mxCreateDoubleScalar(s->warm_starting));
    mxSetField(mx, 0, "scaling",               mxCreateDoubleScalar(s->scaling));
    mxSetField(mx, 0, "polishing",             mxCreateDoubleScalar(s->polishing));
    mxSetField(mx, 0, "rho",                   mxCreateDoubleScalar(s->rho));
    mxSetField(mx, 0, "rho_is_vec",            mxCreateDoubleScalar(s->rho_is_vec));
    mxSetField(mx, 0, "sigma",                 mxCreateDoubleScalar(s->sigma));
    mxSetField(mx, 0, "alpha",                 mxCreateDoubleScalar(s->alpha));
    mxSetField(mx, 0, "cg_max_iter",           mxCreateDoubleScalar(s->cg_max_iter));
    mxSetField(mx, 0, "cg_tol_reduction",      mxCreateDoubleScalar(s->cg_tol_reduction));
    mxSetField(mx, 0, "cg_tol_fraction",       mxCreateDoubleScalar(s->cg_tol_fraction));
    mxSetField(mx, 0, "cg_precond",            mxCreateDoubleScalar(s->cg_precond));
    mxSetField(mx, 0, "adaptive_rho",          mxCreateDoubleScalar(s->adaptive_rho));
    mxSetField(mx, 0, "adaptive_rho_interval", mxCreateDoubleScalar(s->adaptive_rho_interval));
    mxSetField(mx, 0, "adaptive_rho_fraction", mxCreateDoubleScalar(s->adaptive_rho_fraction));
    mxSetField(mx, 0, "adaptive_rho_tolerance",mxCreateDoubleScalar(s->adaptive_rho_tolerance));
    mxSetField(mx, 0, "max_iter",              mxCreateDoubleScalar(s->max_iter));
    mxSetField(mx, 0, "eps_abs",               mxCreateDoubleScalar(s->eps_abs));
    mxSetField(mx, 0, "eps_rel",               mxCreateDoubleScalar(s->eps_rel));
    mxSetField(mx, 0, "eps_prim_inf",          mxCreateDoubleScalar(s->eps_prim_inf));
    mxSetField(mx, 0, "eps_dual_inf",          mxCreateDoubleScalar(s->eps_dual_inf));
    mxSetField(mx, 0, "scaled_termination",    mxCreateDoubleScalar(s->scaled_termination));
    mxSetField(mx, 0, "check_termination",     mxCreateDoubleScalar(s->check_termination));
    mxSetField(mx, 0, "check_dualgap",         mxCreateDoubleScalar(s->check_dualgap));
    mxSetField(mx, 0, "time_limit",            mxCreateDoubleScalar(s->time_limit));
    mxSetField(mx, 0, "delta",                 mxCreateDoubleScalar(s->delta));
    mxSetField(mx, 0, "polish_refine_iter",    mxCreateDoubleScalar(s->polish_refine_iter));

    return mx;
}

/* ------------------------------------------------------------------ */
/* MATLAB struct → Settings (field-by-field copy, skip missing fields)  */
/* ------------------------------------------------------------------ */
static void copyMxStructToSettings(const mxArray* mx, OSQPSettings* s) {
    mxArray* f;

    if ((f = mxGetField(mx, 0, "device")))                s->device                = (OSQPInt)mxGetScalar(f);
    if ((f = mxGetField(mx, 0, "linsys_solver")))         s->linsys_solver         = (enum osqp_linsys_solver_type)(OSQPInt)mxGetScalar(f);
    if ((f = mxGetField(mx, 0, "allocate_solution")))     s->allocate_solution     = (OSQPInt)mxGetScalar(f);
    if ((f = mxGetField(mx, 0, "verbose")))               s->verbose               = (OSQPInt)mxGetScalar(f);
    if ((f = mxGetField(mx, 0, "profiler_level")))        s->profiler_level        = (OSQPInt)mxGetScalar(f);
    if ((f = mxGetField(mx, 0, "warm_starting")))         s->warm_starting         = (OSQPInt)mxGetScalar(f);
    if ((f = mxGetField(mx, 0, "scaling")))               s->scaling               = (OSQPInt)mxGetScalar(f);
    if ((f = mxGetField(mx, 0, "polishing")))             s->polishing             = (OSQPInt)mxGetScalar(f);
    if ((f = mxGetField(mx, 0, "rho")))                   s->rho                   = (OSQPFloat)mxGetScalar(f);
    if ((f = mxGetField(mx, 0, "rho_is_vec")))            s->rho_is_vec            = (OSQPInt)mxGetScalar(f);
    if ((f = mxGetField(mx, 0, "sigma")))                 s->sigma                 = (OSQPFloat)mxGetScalar(f);
    if ((f = mxGetField(mx, 0, "alpha")))                 s->alpha                 = (OSQPFloat)mxGetScalar(f);
    if ((f = mxGetField(mx, 0, "cg_max_iter")))           s->cg_max_iter           = (OSQPInt)mxGetScalar(f);
    if ((f = mxGetField(mx, 0, "cg_tol_reduction")))      s->cg_tol_reduction      = (OSQPInt)mxGetScalar(f);
    if ((f = mxGetField(mx, 0, "cg_tol_fraction")))       s->cg_tol_fraction       = (OSQPFloat)mxGetScalar(f);
    if ((f = mxGetField(mx, 0, "cg_precond")))            s->cg_precond            = (osqp_precond_type)(OSQPInt)mxGetScalar(f);
    if ((f = mxGetField(mx, 0, "adaptive_rho")))          s->adaptive_rho          = (OSQPInt)mxGetScalar(f);
    if ((f = mxGetField(mx, 0, "adaptive_rho_interval"))) s->adaptive_rho_interval = (OSQPInt)mxGetScalar(f);
    if ((f = mxGetField(mx, 0, "adaptive_rho_fraction"))) s->adaptive_rho_fraction = (OSQPFloat)mxGetScalar(f);
    if ((f = mxGetField(mx, 0, "adaptive_rho_tolerance")))s->adaptive_rho_tolerance= (OSQPFloat)mxGetScalar(f);
    if ((f = mxGetField(mx, 0, "max_iter")))              s->max_iter              = (OSQPInt)mxGetScalar(f);
    if ((f = mxGetField(mx, 0, "eps_abs")))               s->eps_abs               = (OSQPFloat)mxGetScalar(f);
    if ((f = mxGetField(mx, 0, "eps_rel")))               s->eps_rel               = (OSQPFloat)mxGetScalar(f);
    if ((f = mxGetField(mx, 0, "eps_prim_inf")))          s->eps_prim_inf          = (OSQPFloat)mxGetScalar(f);
    if ((f = mxGetField(mx, 0, "eps_dual_inf")))          s->eps_dual_inf          = (OSQPFloat)mxGetScalar(f);
    if ((f = mxGetField(mx, 0, "scaled_termination")))    s->scaled_termination    = (OSQPInt)mxGetScalar(f);
    if ((f = mxGetField(mx, 0, "check_termination")))     s->check_termination     = (OSQPInt)mxGetScalar(f);
    if ((f = mxGetField(mx, 0, "check_dualgap")))         s->check_dualgap         = (OSQPInt)mxGetScalar(f);
    if ((f = mxGetField(mx, 0, "time_limit")))            s->time_limit            = (OSQPFloat)mxGetScalar(f);
    if ((f = mxGetField(mx, 0, "delta")))                 s->delta                 = (OSQPFloat)mxGetScalar(f);
    if ((f = mxGetField(mx, 0, "polish_refine_iter")))    s->polish_refine_iter    = (OSQPInt)mxGetScalar(f);
}
