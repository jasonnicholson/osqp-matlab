function result = run_osqp_matlab_unit_tests(generateCoverage)
% RUN_OSQP_MATLAB_UNIT_TESTS  Run only pure-MATLAB OSQP unit tests.

  arguments
    generateCoverage (1,1) logical = false;
  end

  result = run_osqp_unit_tests(generateCoverage, 'matlab');

  if nargout == 0
    clear result
  end
end
