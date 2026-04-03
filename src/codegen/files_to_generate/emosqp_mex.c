/*
 * emosqp_mex.c — MATLAB MEX interface for OSQP embedded (codegen) solver
 *
 * Commands:
 *   'solve'                → solve the embedded problem
 *   'update_data_vec'      → update q, l, u vectors
 *   'update_data_mat'      → update P and/or A matrix elements (embedded mode 2)
 *   'warm_start'           → warm start x, y
 *   'cold_start'           → reset iterates to zero
 *   'update_settings'      → update solver settings
 */

#include <string.h>
#include "mex.h"
#include "osqp.h"
#include "workspace.h"


/* ------------------------------------------------------------------ */
/* Utility functions                                                   */
/* ------------------------------------------------------------------ */

static OSQPFloat* copyToCfloatVector(double* vecData, OSQPInt numel) {
    OSQPFloat* out = (OSQPFloat*)mxMalloc(numel * sizeof(OSQPFloat));
    for (OSQPInt i = 0; i < numel; i++)
        out[i] = (OSQPFloat)vecData[i];
    return out;
}

static void castToDoubleArr(OSQPFloat* arr, double* arr_out, OSQPInt len) {
    for (OSQPInt i = 0; i < len; i++)
        arr_out[i] = (double)arr[i];
}

static void setToNaN(double* arr_out, OSQPInt len) {
    for (OSQPInt i = 0; i < len; i++)
        arr_out[i] = mxGetNaN();
}

#if OSQP_EMBEDDED_MODE == 2
static OSQPInt* copyDoubleToCintVector(double* vecData, OSQPInt numel) {
    OSQPInt* out = (OSQPInt*)mxMalloc(numel * sizeof(OSQPInt));
    for (OSQPInt i = 0; i < numel; i++)
        out[i] = (OSQPInt)vecData[i];
    return out;
}
#endif


/* ------------------------------------------------------------------ */
/* mexFunction                                                         */
/* ------------------------------------------------------------------ */
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    char cmd[64];
    if (nrhs < 1 || mxGetString(prhs[0], cmd, sizeof(cmd)))
        mexErrMsgTxt("First input should be a command string less than 64 characters long.");

    /* ---- solve ---- */
    if (!strcmp("solve", cmd)) {
        if (nlhs != 5 || nrhs != 1)
            mexErrMsgTxt("solve: expects 0 inputs and 5 outputs.");

        OSQPInt n, m;
        osqp_get_dimensions(&solver, &m, &n);

        osqp_solve(&solver);

        /* Allocate outputs */
        plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);  /* x */
        plhs[1] = mxCreateDoubleMatrix(m, 1, mxREAL);  /* y */
        plhs[2] = mxCreateDoubleScalar(solver.info->status_val);
        plhs[3] = mxCreateDoubleScalar(solver.info->iter);
        plhs[4] = mxCreateDoubleScalar(solver.info->run_time);

        OSQPInt sv = solver.info->status_val;
        if (sv != OSQP_PRIMAL_INFEASIBLE &&
            sv != OSQP_PRIMAL_INFEASIBLE_INACCURATE &&
            sv != OSQP_DUAL_INFEASIBLE &&
            sv != OSQP_DUAL_INFEASIBLE_INACCURATE) {
            castToDoubleArr(solver.solution->x, mxGetPr(plhs[0]), n);
            castToDoubleArr(solver.solution->y, mxGetPr(plhs[1]), m);
        } else {
            setToNaN(mxGetPr(plhs[0]), n);
            setToNaN(mxGetPr(plhs[1]), m);
        }
        return;
    }

    /* ---- update_data_vec ---- */
    if (!strcmp("update_data_vec", cmd)) {
        const mxArray* mxq = prhs[1];
        const mxArray* mxl = prhs[2];
        const mxArray* mxu = prhs[3];

        OSQPInt n, m;
        osqp_get_dimensions(&solver, &m, &n);

        OSQPFloat* q_new = NULL;
        OSQPFloat* l_new = NULL;
        OSQPFloat* u_new = NULL;
        if (!mxIsEmpty(mxq)) q_new = copyToCfloatVector(mxGetPr(mxq), n);
        if (!mxIsEmpty(mxl)) l_new = copyToCfloatVector(mxGetPr(mxl), m);
        if (!mxIsEmpty(mxu)) u_new = copyToCfloatVector(mxGetPr(mxu), m);

        osqp_update_data_vec(&solver, q_new, l_new, u_new);

        if (q_new) mxFree(q_new);
        if (l_new) mxFree(l_new);
        if (u_new) mxFree(u_new);
        return;
    }

#if OSQP_EMBEDDED_MODE == 2
    /* ---- update_data_mat ---- */
    if (!strcmp("update_data_mat", cmd)) {
        const mxArray* mxPx     = prhs[1];
        const mxArray* mxPx_idx = prhs[2];
        OSQPInt        Px_n     = (OSQPInt)mxGetScalar(prhs[3]);
        const mxArray* mxAx     = prhs[4];
        const mxArray* mxAx_idx = prhs[5];
        OSQPInt        Ax_n     = (OSQPInt)mxGetScalar(prhs[6]);

        OSQPFloat* Px     = NULL;
        OSQPInt*   Px_idx = NULL;
        OSQPFloat* Ax     = NULL;
        OSQPInt*   Ax_idx = NULL;

        if (!mxIsEmpty(mxPx))     Px     = copyToCfloatVector(mxGetPr(mxPx), Px_n);
        if (!mxIsEmpty(mxPx_idx)) Px_idx = copyDoubleToCintVector(mxGetPr(mxPx_idx), Px_n);
        if (!mxIsEmpty(mxAx))     Ax     = copyToCfloatVector(mxGetPr(mxAx), Ax_n);
        if (!mxIsEmpty(mxAx_idx)) Ax_idx = copyDoubleToCintVector(mxGetPr(mxAx_idx), Ax_n);

        osqp_update_data_mat(&solver, Px, Px_idx, Px_n, Ax, Ax_idx, Ax_n);

        if (Px)     mxFree(Px);
        if (Px_idx) mxFree(Px_idx);
        if (Ax)     mxFree(Ax);
        if (Ax_idx) mxFree(Ax_idx);
        return;
    }
#endif

    /* ---- warm_start ---- */
    if (!strcmp("warm_start", cmd)) {
        const mxArray* mxX = prhs[1];
        const mxArray* mxY = prhs[2];

        OSQPInt n, m;
        osqp_get_dimensions(&solver, &m, &n);

        OSQPFloat* x_vec = NULL;
        OSQPFloat* y_vec = NULL;
        if (!mxIsEmpty(mxX)) x_vec = copyToCfloatVector(mxGetPr(mxX), n);
        if (!mxIsEmpty(mxY)) y_vec = copyToCfloatVector(mxGetPr(mxY), m);

        osqp_warm_start(&solver, x_vec, y_vec);

        if (x_vec) mxFree(x_vec);
        if (y_vec) mxFree(y_vec);
        return;
    }

    /* ---- cold_start ---- */
    if (!strcmp("cold_start", cmd)) {
        osqp_cold_start(&solver);
        return;
    }

    /* ---- update_settings ---- */
    if (!strcmp("update_settings", cmd)) {
        /* Accept a settings struct and update specific fields */
        const mxArray* mxSettings = prhs[1];
        OSQPSettings new_settings;
        memcpy(&new_settings, solver.settings, sizeof(OSQPSettings));

        mxArray* f;
        if ((f = mxGetField(mxSettings, 0, "max_iter")))
            new_settings.max_iter = (OSQPInt)mxGetScalar(f);
        if ((f = mxGetField(mxSettings, 0, "eps_abs")))
            new_settings.eps_abs = (OSQPFloat)mxGetScalar(f);
        if ((f = mxGetField(mxSettings, 0, "eps_rel")))
            new_settings.eps_rel = (OSQPFloat)mxGetScalar(f);
        if ((f = mxGetField(mxSettings, 0, "eps_prim_inf")))
            new_settings.eps_prim_inf = (OSQPFloat)mxGetScalar(f);
        if ((f = mxGetField(mxSettings, 0, "eps_dual_inf")))
            new_settings.eps_dual_inf = (OSQPFloat)mxGetScalar(f);
        if ((f = mxGetField(mxSettings, 0, "rho")))
            new_settings.rho = (OSQPFloat)mxGetScalar(f);
        if ((f = mxGetField(mxSettings, 0, "alpha")))
            new_settings.alpha = (OSQPFloat)mxGetScalar(f);
        if ((f = mxGetField(mxSettings, 0, "warm_starting")))
            new_settings.warm_starting = (OSQPInt)mxGetScalar(f);
        if ((f = mxGetField(mxSettings, 0, "polishing")))
            new_settings.polishing = (OSQPInt)mxGetScalar(f);
        if ((f = mxGetField(mxSettings, 0, "verbose")))
            new_settings.verbose = (OSQPInt)mxGetScalar(f);
        if ((f = mxGetField(mxSettings, 0, "scaled_termination")))
            new_settings.scaled_termination = (OSQPInt)mxGetScalar(f);
        if ((f = mxGetField(mxSettings, 0, "check_termination")))
            new_settings.check_termination = (OSQPInt)mxGetScalar(f);
        if ((f = mxGetField(mxSettings, 0, "time_limit")))
            new_settings.time_limit = (OSQPFloat)mxGetScalar(f);

        osqp_update_settings(&solver, &new_settings);
        return;
    }

    mexErrMsgTxt("Command not recognized.");
}
