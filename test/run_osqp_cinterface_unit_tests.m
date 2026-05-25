function result = run_osqp_cinterface_unit_tests(generateCoverage)
% RUN_OSQP_CINTERFACE_UNIT_TESTS  Run only CInterface-focused OSQP unit tests.

  arguments
    generateCoverage (1,1) logical = false;
  end

  result = run_osqp_unit_tests(generateCoverage, 'cinterface');

  if nargout == 0
    clear result
  end
end
