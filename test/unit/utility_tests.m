classdef utility_tests < matlab.unittest.TestCase
    % UTILITY_TESTS  Tests for osqp package-level utility functions.

    methods (Test)
        function test_default_settings(testCase)
            opts = osqp.default_settings();
            testCase.verifyClass(opts, 'osqp.SolverOptions');
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

        function test_options_fromStruct_invalid_value_skipped(testCase)
            % The catch block in fromStruct should silently skip bad values
            s.rho = -999;  % Invalid: must be positive
            s.max_iter = 500;  % Valid
            opts = osqp.SolverOptions.fromStruct(s);
            % max_iter should be set, rho should keep default
            testCase.verifyEqual(opts.max_iter, 500);
            testCase.verifyEqual(opts.rho, 0.1);  % default, -999 was skipped
        end
    end
end
