classdef codegen_mat_tests < matlab.unittest.TestCase
    % CODEGEN_MAT_TESTS Test code generation (embedded mode 2) for OSQP v1.0.0

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
            testCase.P = sparse([11 0; 0 0.1]);
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
        function test_codegen_mat(testCase)
            testCase.solver.codegen(testCase.codegen_dir, ...
                'parameters', 'matrices', 'force_rewrite', true);

            % Verify generated files exist (v1.0.0 layout)
            testCase.verifyTrue(exist(fullfile(testCase.codegen_dir, ...
                'Makefile'), 'file') > 0);
            testCase.verifyTrue(exist(fullfile(testCase.codegen_dir, ...
                'inc', 'public', 'osqp.h'), 'file') > 0);
            testCase.verifyTrue(exist(fullfile(testCase.codegen_dir, ...
                'src', 'osqp_api.c'), 'file') > 0);

            % Verify workspace/configure/emosqp files land in target dir
            ws = dir(fullfile(testCase.codegen_dir, '*workspace*'));
            testCase.verifyGreaterThan(numel(ws), 0, ...
                'Workspace files missing from codegen output');
            cf = dir(fullfile(testCase.codegen_dir, '*configure*'));
            testCase.verifyGreaterThan(numel(cf), 0, ...
                'Configure file missing from codegen output');
            em = dir(fullfile(testCase.codegen_dir, '*emosqp*'));
            testCase.verifyGreaterThan(numel(em), 0, ...
                'Emosqp file missing from codegen output');
        end

        function test_codegen_mat_with_prefix(testCase)
            testCase.solver.codegen(testCase.codegen_dir, ...
                'parameters', 'matrices', 'force_rewrite', true, ...
                'prefix', 'myqp_');

            % Verify prefixed workspace files exist inside target dir
            ws = dir(fullfile(testCase.codegen_dir, 'myqp_*'));
            testCase.verifyGreaterThan(numel(ws), 0, ...
                'Prefixed files missing from codegen output');

            % Verify no stray files outside target dir
            parent = fileparts(testCase.codegen_dir);
            dname = regexp(testCase.codegen_dir, '[^/\\]+$', 'match', 'once');
            stray = dir(fullfile(parent, [dname '*']));
            stray = stray(~ismember({stray.name}, {dname}));
            testCase.verifyEmpty(stray, ...
                'Stray codegen files found outside target directory');
        end
    end
end
