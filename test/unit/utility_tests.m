classdef utility_tests < matlab.unittest.TestCase
    % UTILITY_TESTS  Tests for osqp package-level utility functions.

    methods (Test)
        function test_default_settings(testCase)
            opts = osqp.default_settings();
            testCase.verifyClass(opts, 'osqp.Options');
            testCase.verifyEqual(opts.rho, 0.1);
            testCase.verifyEqual(opts.max_iter, 4000);
        end

        function test_version(testCase)
            v = osqp.version();
            testCase.verifyTrue(ischar(v) || isstring(v));
            testCase.verifyTrue(contains(v, '1.0.0'));
        end

        function test_capabilities(testCase)
            cap = osqp.capabilities();
            testCase.verifyTrue(isnumeric(cap));
            testCase.verifyGreaterThanOrEqual(cap, 0);
        end

        function test_constant_dot_access(testCase)
            testCase.verifyEqual(osqp.constant.OSQP_INFTY, 1e30);
            testCase.verifyEqual(osqp.constant.OSQP_SOLVED, 1);
            testCase.verifyEqual(osqp.constant.OSQP_PRIMAL_INFEASIBLE, 3);
            testCase.verifyEqual(osqp.constant.OSQP_DUAL_INFEASIBLE, 5);
            testCase.verifyEqual(osqp.constant.OSQP_MAX_ITER_REACHED, 7);
            testCase.verifyEqual(osqp.constant.OSQP_TIME_LIMIT_REACHED, 8);
            testCase.verifyEqual(osqp.constant.OSQP_NON_CONVEX, 9);
            testCase.verifyEqual(osqp.constant.OSQP_UNSOLVED, 11);
        end

        function test_constant_polish_values(testCase)
            testCase.verifyEqual(osqp.constant.OSQP_POLISH_SUCCESS, 1);
            testCase.verifyEqual(osqp.constant.OSQP_POLISH_FAILED, -1);
            testCase.verifyEqual(osqp.constant.OSQP_POLISH_NOT_PERFORMED, 0);
            testCase.verifyEqual(osqp.constant.OSQP_POLISH_NO_ACTIVE_SET_FOUND, 2);
            testCase.verifyEqual(osqp.constant.OSQP_POLISH_LINSYS_ERROR, -2);
        end

        function test_constant_lookup(testCase)
            testCase.verifyEqual(osqp.constant.lookup('OSQP_SOLVED'), 1);
            testCase.verifyEqual(osqp.constant.lookup('OSQP_INFTY'), 1e30);
        end

        function test_constant_lookup_unknown(testCase)
            testCase.verifyError( ...
                @() osqp.constant.lookup('NOT_A_CONSTANT'), '');
        end

        function test_constant_function_call_errors(testCase)
            testCase.verifyError( ...
                @() osqp.constant('OSQP_SOLVED'), ...
                'OSQP:constant');
        end

        function test_options_fromStruct_invalid_value_skipped(testCase)
            % The catch block in fromStruct should silently skip bad values
            s.rho = -999;  % Invalid: must be positive
            s.max_iter = 500;  % Valid
            opts = osqp.Options.fromStruct(s);
            % max_iter should be set, rho should keep default
            testCase.verifyEqual(opts.max_iter, 500);
            testCase.verifyEqual(opts.rho, 0.1);  % default, -999 was skipped
        end
    end
end
