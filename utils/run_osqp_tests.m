function result = run_osqp_tests(generateCoverage)
% RUN_OSQP_TESTS  Run OSQP MATLAB unit tests with optional HTML coverage report.
%
%   result = run_osqp_tests()          — run all tests, no coverage
%   result = run_osqp_tests(true)      — run all tests + HTML coverage report
%
%   The coverage report (when enabled) is written to coverage_report/ and
%   covers osqp.m plus everything under codegen/ and utils/.

import matlab.unittest.TestSuite
import matlab.unittest.TestRunner
import matlab.unittest.plugins.CodeCoveragePlugin
import matlab.unittest.plugins.codecoverage.CoverageReport

if nargin < 1
    generateCoverage = false;
end

% Get the directory containing this file (= osqp-matlab root)
osqp_path = fileparts(mfilename('fullpath'));
if isempty(osqp_path)
    osqp_path = pwd;
end

% Ensure osqp is on the path
addpath(osqp_path);
unittest_dir = fullfile(osqp_path, 'unittests');

% Build test suite from the unittests folder
suite = TestSuite.fromFolder(unittest_dir);

if generateCoverage

    % Collect source files to measure coverage on
    sourceFiles = {fullfile(osqp_path, 'osqp.m')};

    % Add codegen .m files
    codegenFiles = dir(fullfile(osqp_path, 'codegen', '*.m'));
    for k = 1:numel(codegenFiles)
        sourceFiles{end+1} = fullfile(codegenFiles(k).folder, codegenFiles(k).name); %#ok<AGROW>
    end

    % Add utils .m files
    utilFiles = dir(fullfile(osqp_path, 'utils', '*.m'));
    for k = 1:numel(utilFiles)
        sourceFiles{end+1} = fullfile(utilFiles(k).folder, utilFiles(k).name); %#ok<AGROW>
    end

    reportDir = fullfile(osqp_path, 'coverage_report');
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
