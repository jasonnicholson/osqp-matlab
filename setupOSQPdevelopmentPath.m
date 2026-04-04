function setupOSQPdevelopmentPath(shouldSetup)
  %SETUPOSQPDEVELOPMENTPATH Setup OSQP development paths
  %   Adds or removes the OSQP development source and utility folders to/from
  %   the MATLAB path depending on the input flag.
  arguments (Input)
    shouldSetup (1,1) logical = true;
  end

  repoRoot = fileparts(mfilename('fullpath'));

  p = convertStringsToChars([...
    fullfile(repoRoot, "src"), ...
    fullfile(repoRoot, "utils"), ...
    fullfile(repoRoot, "test"), ...
    fullfile(repoRoot, "test", "unit"), ...
    fullfile(repoRoot, "test", "integration"), ...
    fullfile(repoRoot, "examples"), ...
  ]);

  if shouldSetup
    addpath(p{:});
  else
    rmpath(p{:});
  end
end