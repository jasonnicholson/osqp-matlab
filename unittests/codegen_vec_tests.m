classdef codegen_vec_tests < matlab.unittest.TestCase
    % CODEGEN_VEC_TESTS Test code generation (embedded mode 1) for OSQP v1.0.0

    properties
        P
        q
        A
        u
        l
        solver
        tol
        codegen_dir
    end

    methods(TestMethodSetup)
        function setup_problem(testCase)
            testCase.P = sparse([11 0; 0 0]);
            testCase.q = [3; 4];
            testCase.A = sparse([-1 0; 0 -1; -1 -3; 2 5; 3 4]);
            testCase.u = [0; 0; -15; 100; 80];
            testCase.l = -1e30 * ones(5, 1);
            testCase.tol = 1e-4;
            testCase.codegen_dir = tempname;

            testCase.solver = osqp;
            testCase.solver.setup(testCase.P, testCase.q, ...
                testCase.A, testCase.l, testCase.u, ...
                'verbose', false, 'polishing', true, ...
                'eps_abs', 1e-8, 'eps_rel', 1e-8);
        end
    end

    methods(TestMethodTeardown)
        function teardown(testCase)
            delete(testCase.solver);
            if exist(testCase.codegen_dir, 'dir')
                rmdir(testCase.codegen_dir, 's');
            end
        end
    end

    methods(Test)
        function test_codegen_vec(testCase)
            testCase.solver.codegen(testCase.codegen_dir, ...
                'parameters', 'vectors', 'force_rewrite', true);

            % Verify generated files exist (v1.0.0 layout)
            testCase.verifyTrue(exist(fullfile(testCase.codegen_dir, ...
                'Makefile'), 'file') > 0);
            testCase.verifyTrue(exist(fullfile(testCase.codegen_dir, ...
                'inc', 'public', 'osqp.h'), 'file') > 0);
            testCase.verifyTrue(exist(fullfile(testCase.codegen_dir, ...
                'src', 'osqp_api.c'), 'file') > 0);
        end
    end
end
