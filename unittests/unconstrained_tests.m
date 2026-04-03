classdef unconstrained_tests < matlab.unittest.TestCase
    % UNCONSTRAINED_TESTS Test unconstrained QP for OSQP v1.0.0

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
            testCase.m = 0;
            Pt = sprandn(testCase.n, testCase.n, 0.6);
            testCase.P = Pt * Pt' + speye(testCase.n);
            testCase.q = randn(testCase.n, 1);
            testCase.A = sparse(0, testCase.n);
            testCase.l = zeros(0, 1);
            testCase.u = zeros(0, 1);
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
        function test_unconstrained_problem(testCase)
            [x_ref, y_ref, obj_ref] = load_high_accuracy('test_unconstrained_problem');
            res = testCase.solver.solve();
            testCase.verifyEqual(res.x, x_ref, 'AbsTol', testCase.tol);
            testCase.verifyEqual(res.info.obj_val, obj_ref, 'AbsTol', testCase.tol);
        end
    end
end
