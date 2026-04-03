function run_osqp_examples()
% Runs all OSQP MATLAB example scripts located in the 'examples' folder.

  repoRoot = fileparts(mfilename('fullpath'));
  startingDirectory = pwd;
  
  examples = dir(fullfile(repoRoot, 'examples', '*.m'));
  examples = string({examples.name});
  nExamples = numel(examples);

  for i = 1:nExamples
    fprintf('\n==== Running %d of %d: %s ====\n', i, nExamples, examples(i));
    try
      cd(fullfile(repoRoot, 'examples'));
      run(examples(i));
      cd(startingDirectory);
      fprintf('==== %s: OK ====\n', examples(i));
    catch e
      cd(startingDirectory);
      fprintf('==== %s: FAILED ====\n%s\n', examples(i), e.message);
    end
  end
end