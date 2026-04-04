classdef options_tests < matlab.unittest.TestCase
    % OPTIONS_TESTS  Unit tests for osqp.Options

    methods (Test)
        function testDefaultConstruction(testCase)
            opts = osqp.Options();
            testCase.verifyEqual(opts.rho, 0.1);
            testCase.verifyEqual(opts.sigma, 1e-6);
            testCase.verifyEqual(opts.alpha, 1.6);
            testCase.verifyEqual(opts.max_iter, 4000);
            testCase.verifyEqual(opts.eps_abs, 1e-3);
            testCase.verifyEqual(opts.eps_rel, 1e-3);
            testCase.verifyEqual(opts.eps_prim_inf, 1e-4);
            testCase.verifyEqual(opts.eps_dual_inf, 1e-4);
            testCase.verifyEqual(opts.scaling, 10);
            testCase.verifyEqual(opts.adaptive_rho, true);
            testCase.verifyEqual(opts.adaptive_rho_interval, 0);
            testCase.verifyEqual(opts.adaptive_rho_tolerance, 5);
            testCase.verifyEqual(opts.adaptive_rho_fraction, 0.4);
            testCase.verifyEqual(opts.polishing, false);
            testCase.verifyEqual(opts.delta, 1e-6);
            testCase.verifyEqual(opts.polish_refine_iter, 3);
            testCase.verifyEqual(opts.verbose, true);
            testCase.verifyEqual(opts.warm_starting, true);
            testCase.verifyEqual(opts.check_termination, 25);
            testCase.verifyEqual(opts.time_limit, 0);
            testCase.verifyEqual(opts.linear_solver, 'qdldl_c');
        end

        function testPropertyModification(testCase)
            opts = osqp.Options();
            opts.max_iter = 1000;
            testCase.verifyEqual(opts.max_iter, 1000);
            opts.rho = 0.5;
            testCase.verifyEqual(opts.rho, 0.5);
            opts.linear_solver = 'qdldl';
            testCase.verifyEqual(opts.linear_solver, 'qdldl');
            opts.verbose = false;
            testCase.verifyEqual(opts.verbose, false);
        end

        function testToStruct(testCase)
            opts = osqp.Options();
            opts.max_iter = 500;
            s = opts.toStruct();
            testCase.verifyTrue(isstruct(s));
            testCase.verifyEqual(s.max_iter, 500);
            testCase.verifyEqual(s.rho, 0.1);
            testCase.verifyEqual(s.linear_solver, 'qdldl_c');
        end

        function testFromStruct(testCase)
            s.rho = 0.5;
            s.max_iter = 2000;
            s.linear_solver = 'qdldl';
            opts = osqp.Options.fromStruct(s);
            testCase.verifyEqual(opts.rho, 0.5);
            testCase.verifyEqual(opts.max_iter, 2000);
            testCase.verifyEqual(opts.linear_solver, 'qdldl');
            % Unset fields keep defaults
            testCase.verifyEqual(opts.sigma, 1e-6);
        end

        function testToStructRoundTrip(testCase)
            opts = osqp.Options();
            opts.rho = 0.3;
            opts.alpha = 1.2;
            opts.max_iter = 100;
            s = opts.toStruct();
            opts2 = osqp.Options.fromStruct(s);
            testCase.verifyEqual(opts2.rho, opts.rho);
            testCase.verifyEqual(opts2.alpha, opts.alpha);
            testCase.verifyEqual(opts2.max_iter, opts.max_iter);
        end

        function testValidationRhoPositive(testCase)
            opts = osqp.Options();
            testCase.verifyError(@() setfield_helper(opts, 'rho', -1), ...
                'MATLAB:validators:mustBePositive');
        end

        function testValidationAlphaRange(testCase)
            opts = osqp.Options();
            testCase.verifyError(@() setfield_helper(opts, 'alpha', 3), ...
                'MATLAB:validators:mustBeInRange');
            testCase.verifyError(@() setfield_helper(opts, 'alpha', -0.5), ...
                'MATLAB:validators:mustBeInRange');
        end

        function testValidationLinearSolver(testCase)
            opts = osqp.Options();
            testCase.verifyError(@() setfield_helper(opts, 'linear_solver', 'bad_solver'), ...
                'MATLAB:validators:mustBeMember');
        end

        function testValidationMaxIterInteger(testCase)
            opts = osqp.Options();
            testCase.verifyError(@() setfield_helper(opts, 'max_iter', 3.5), ...
                'MATLAB:validators:mustBeInteger');
        end

        function testValidationEpsNonnegative(testCase)
            opts = osqp.Options();
            testCase.verifyError(@() setfield_helper(opts, 'eps_abs', -0.01), ...
                'MATLAB:validators:mustBeNonnegative');
        end

        function testUpdatableClassification(testCase)
            opts = osqp.Options();
            testCase.verifyTrue(opts.isUpdatable('max_iter'));
            testCase.verifyTrue(opts.isUpdatable('rho'));
            testCase.verifyTrue(opts.isUpdatable('alpha'));
            testCase.verifyFalse(opts.isUpdatable('sigma'));
            testCase.verifyFalse(opts.isUpdatable('scaling'));
            testCase.verifyFalse(opts.isUpdatable('linear_solver'));
        end

        function testSetupOnlyClassification(testCase)
            opts = osqp.Options();
            testCase.verifyTrue(opts.isSetupOnly('sigma'));
            testCase.verifyTrue(opts.isSetupOnly('scaling'));
            testCase.verifyTrue(opts.isSetupOnly('linear_solver'));
            testCase.verifyFalse(opts.isSetupOnly('max_iter'));
            testCase.verifyFalse(opts.isSetupOnly('rho'));
        end

        function testUpdatableAndSetupOnlyAreMutuallyExclusive(testCase)
            % Every non-constant property should be in exactly one list
            opts = osqp.Options();
            mc = metaclass(opts);
            for k = 1:numel(mc.PropertyList)
                p = mc.PropertyList(k);
                if p.Constant || p.Dependent
                    continue
                end
                u = opts.isUpdatable(p.Name);
                s = opts.isSetupOnly(p.Name);
                testCase.verifyTrue(xor(u, s), ...
                    sprintf('Property ''%s'' must be in exactly one list', p.Name));
            end
        end
    end
end


function setfield_helper(obj, name, val)
    obj.(name) = val;
end
