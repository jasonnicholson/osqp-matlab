function result = run_osqp_unit_tests(generateCoverage)
  % RUN_OSQP_UNIT_TESTS  Run OSQP MATLAB unit tests with optional HTML coverage report.
  %
  %   result = run_osqp_tests()          — run all tests, no coverage
  %   result = run_osqp_tests(true)      — run all tests + HTML coverage report
  %
  %   The coverage report (when enabled) is written to coverage_report/ and
  %   covers osqp.m plus everything under codegen/ and utils/.


  arguments
    generateCoverage (1,1) logical = true;
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

  % Build test suite from the unittests folder
  suite = TestSuite.fromFolder(unittest_dir);

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
