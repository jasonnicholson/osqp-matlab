function setupOSQPdevelopmentPath(shouldSetup)
  %SETUPOSQPDEVELOPMENTPATH Summary of this function goes here
  %   Detailed explanation goes here
  arguments (Input)
    shouldSetup (1,1) logical = true;
  end

  thisFolder = fileparts(mfilename('fullpath'));

  p = convertStringsToChars([...
    fullfile(thisFolder, "src"), ...
    fullfile(thisFolder, "utils"), ...
  ]);

  if shouldSetup
    addpath(p{:});
  else
    rmpath(p{:});
  end
end