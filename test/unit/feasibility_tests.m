classdef feasibility_tests < matlab.unittest.TestCase
    % FEASIBILITY_TESTS Test feasibility problems for OSQP v1.0.0

    properties
        P
        q
        A
        u
        l
        m
        n
        solver
        tol
    end

    methods(TestMethodSetup)
        function setup_problem(testCase)
            rng(4);
            testCase.n = 30;
            testCase.m = 30;
            testCase.P = sparse(testCase.n, testCase.n);
            testCase.q = zeros(testCase.n, 1);
            testCase.A = sprandn(testCase.m, testCase.n, 0.8);
            b = randn(testCase.m, 1);
            testCase.l = b;
            testCase.u = b;
            testCase.tol = 1e-4;

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
        end
    end

    methods(Test)
        function test_feasibility_problem(testCase)
            [x_ref, y_ref, obj_ref] = load_high_accuracy('test_feasibility_problem');
            res = testCase.solver.solve();
            testCase.verifyEqual(res.x, x_ref, 'AbsTol', testCase.tol);
            testCase.verifyEqual(res.y, y_ref, 'AbsTol', testCase.tol);
            testCase.verifyEqual(res.info.obj_val, obj_ref, 'AbsTol', testCase.tol);
        end
    end
end
