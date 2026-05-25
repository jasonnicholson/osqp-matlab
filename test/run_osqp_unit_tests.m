function result = run_osqp_unit_tests(generateCoverage, scope)
  % RUN_OSQP_UNIT_TESTS  Run OSQP unit tests with optional backend scope and coverage.
  %
  %   result = run_osqp_unit_tests()                    % all unit tests, with coverage
  %   result = run_osqp_unit_tests(true)                % all unit tests, with coverage
  %   result = run_osqp_unit_tests(false, 'cinterface') % CInterface-focused tests only
  %   result = run_osqp_unit_tests(false, 'matlab')     % pure-MATLAB solver tests only
  %
  %   The coverage report (when enabled) is written to coverage_report/ and
  %   covers osqp.m plus everything under codegen/ and utils/.


  arguments
    generateCoverage (1,1) logical = true;
    scope (1,:) char {mustBeMember(scope, {'all', 'cinterface', 'matlab'})} = 'all';
  end

  import matlab.unittest.TestSuite
  import matlab.unittest.TestRunner
  import matlab.unittest.plugins.CodeCoveragePlugin
  import matlab.unittest.plugins.codecoverage.CoverageReport

  % Get the directory containing this file (= osqp-matlab root)
  osqp_root = fileparts(fileparts(mfilename('fullpath')));

  srcFolder = fullfile(osqp_root, "src");

  % Ensure osqp is on the path
  assert(~isempty(which('osqp')), "osqp.m not found on the MATLAB path. Run setupOSQPdevelopmentPath first.");

  unittest_dir = fullfile(osqp_root, 'test', 'unit');

  cinterfaceTests = { ...
    'basic_tests.m', ...
    'codegen_mat_tests.m', ...
    'codegen_vec_tests.m', ...
    'derivative_tests.m', ...
    'dual_infeasibility_tests.m', ...
    'feasibility_tests.m', ...
    'non_cvx_tests.m', ...
    'polishing_tests.m', ...
    'primal_infeasibility_tests.m', ...
    'unconstrained_tests.m', ...
    'update_matrices_tests.m', ...
    'utility_tests.m', ...
    'warm_start_tests.m'};

  matlabTests = { ...
    'linsys_tests.m', ...
    'options_tests.m', ...
    'qdldl_solver_tests.m', ...
    'solver_cparity_tests.m', ...
    'solver_coverage_tests.m', ...
    'solver_tests.m'};

  switch scope
    case 'cinterface'
      selected = cinterfaceTests;
    case 'matlab'
      selected = matlabTests;
    otherwise
      selected = [cinterfaceTests, matlabTests];
  end

  suite = matlab.unittest.Test.empty;
  for k = 1:numel(selected)
    suite = [suite, TestSuite.fromFile(fullfile(unittest_dir, selected{k}))]; %#ok<AGROW>
  end

  if generateCoverage

    % Collect source files to measure coverage on
    sourceFiles = fullfile(srcFolder, "osqp.m");

    % Add +osqp package files
    pkgFiles = dir(fullfile(srcFolder, '+osqp', '*.m'));
    for k = 1:numel(pkgFiles)
      sourceFiles(end+1) = string(fullfile(pkgFiles(k).folder, pkgFiles(k).name)); %#ok<AGROW>
    end

    % Add +osqp/+linsys package files
    linsysFiles = dir(fullfile(srcFolder, '+osqp', '+linsys', '*.m'));
    for k = 1:numel(linsysFiles)
      sourceFiles(end+1) = string(fullfile(linsysFiles(k).folder, linsysFiles(k).name)); %#ok<AGROW>
    end

    % Add codegen .m files
    codegenFiles = dir(fullfile(srcFolder, 'codegen', '*.m'));
    for k = 1:numel(codegenFiles)
      sourceFiles(end+1) = string(fullfile(codegenFiles(k).folder, codegenFiles(k).name)); %#ok<AGROW>
    end

    reportDir = fullfile(osqp_root, 'coverage_report');
    reportFormat = CoverageReport(reportDir);
    coveragePlugin = CodeCoveragePlugin.forFile(sourceFiles, 'Producing', reportFormat);

    runner = TestRunner.withTextOutput;
    runner.addPlugin(coveragePlugin);
    result = runner.run(suite);

    fprintf('\nHTML coverage report: %s\n', fullfile(reportDir, 'index.html'));
  else
    result = run(suite);
  end

  % Print summary
  fprintf('\n=== Test Summary ===\n');
  fprintf('  Total:   %d\n', numel(result));
  fprintf('  Passed:  %d\n', sum([result.Passed]));
  fprintf('  Failed:  %d\n', sum([result.Failed]));
  fprintf('  Errors:  %d\n', sum([result.Incomplete]));

  if nargout == 0
    clear result
  end

end
