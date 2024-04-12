function make_osqp(command, verboseFlag, debugFlag,cmakeExtraOptions, mexExtraOptions)
  % Matlab MEX makefile for OSQP.
  %
  %   make_osqp
  %   make_osqp(command)
  %   make_osqp(command, verboseFlag)
  %   make_osqp(command, verboseFlag, debugFlag)
  %
  %% Inputs
  % - command - string scalar or array. Must be one or more of the following:
  %     - "all" - compile osqp and osqp_mex. Call with this command alone.
  %     - "codegen" - copy source files for codegen.
  %     - "osqp" - builds the OSQP solver using CMake
  %     - "osqp_mex" - builds the OSQP mex interface and links it to the OSQP library
  %     - "clean" - delete all object files (.o and .obj). Call with this command alone.
  %     - "purge" - same as clean, and also delete the mex files. Call with this command alone.
  %       DEFAULT: "all"
  % - verboseFlag - logical scalar. Sets extra communication to the user.
  %       DEFAULT: false
  % - debugFlag - logical scalar. Enables debug symbols in the mex function and shared libraries.
  %       DEFAULT: false
  %
  %% Description
  %    MAKE_OSQP is a make file for OSQP solver. It builds OSQP and its components from source.
  %

  %% Input checking
  arguments
    command (:,1) string {checkCommandString(command)} = "all";
    verboseFlag (1,1) logical = false;
    debugFlag (1,1) logical = false;
    cmakeExtraOptions (1,1) string = "";
    mexExtraOptions (1,1) string = "";
  end

  %% Try to unlock any pre-existing version of osqp_mex

  % this prevents compile errors if a user builds, runs osqp
  % and then tries to recompile
  if(mislocked('osqp_mex'))
    munlock('osqp_mex');
  end

  %% Basic compile commands

  % Get make and mex commands
  make_cmd = "cmake --build .";
  if ~verboseFlag
    mexoptflags = "-silent";
  else
    mexoptflags = "-v";
  end

  % Add arguments to cmake and mex compiler
  cmake_args = "-DMATLAB=ON";
  mexoptflags = [mexoptflags; "-DMATLAB"];

  % Add specific generators for windows linux or mac
  if (ispc)
    cmake_args = sprintf("%s %s", cmake_args, '-G "MinGW Makefiles"');
  else
    cmake_args = sprintf("%s %s", cmake_args, '-G "Unix Makefiles"');
  end

  % Pass Matlab root to cmake
  Matlab_ROOT = strrep(matlabroot, '\', '/');
  cmake_args = sprintf('%s %s%s%s', cmake_args, '-DMatlab_ROOT_DIR="', Matlab_ROOT, '"');

  % Configure debug symbols
  if debugFlag
    cmake_args = sprintf("%s -DCMAKE_BUILD_TYPE=Debug",cmake_args);
  end

  % Add extra options to cmake
  cmake_args = sprintf("%s %s", cmake_args, cmakeExtraOptions);

  % Add parameters options to mex and cmake
  % CTRLC
  mex_libs = "-lut";
  if (ispc)
      % not sure why this is necessary but it is.
    folderContainingLibut_lib = fullfile(matlabroot,"extern","lib","win64","mingw64");
    mex_libs = [mex_libs; "-L"+folderContainingLibut_lib];
  end
  % Shared library loading
  if (isunix && ~ismac)
    mex_libs = [mex_libs; "-ldl"];
  end

  % Add large arrays support if computer is 64 bit and a pre-2018 version Release R2018a corresponds to Matlab version 9.4
  if (contains(computer, '64') && verLessThan('matlab', '9.4'))
    mexoptflags = [mexoptflags; '-largeArrayDims'];
  end

  %Force Matlab to respect usage of mxGetPr in releases after 2018a, Note -R2017b signifies the mex API.
  if ~verLessThan('matlab', '9.4')
    mexoptflags = [mexoptflags; '-R2017b'];
  end

  % Set optimizer flag
  if (~ispc)
    mexoptflags = [mexoptflags; "COPTIMFLAGS='-O3'"];
  end

  % configure mex debug symbols
  if debugFlag
    mexoptflags = [mexoptflags; "-g"];
  end

  % Add extra options to mex
  mexoptflags = [mexoptflags; mexExtraOptions];

  % Set library extension
  lib_ext = '.a';
  lib_name = sprintf('libosqp%s', lib_ext);

  % Set osqp directory and osqp_build directory
  makefile_path = fileparts(mfilename("fullpath"));
  osqp_dir = fullfile(makefile_path, 'osqp_sources');
  osqp_build_dir = fullfile(osqp_dir, 'build');
  qdldl_dir = fullfile(osqp_dir, 'lin_sys', 'direct', 'qdldl');
  cg_sources_dir = fullfile(makefile_path, 'codegen', 'sources');

  % Include directory
  inc_dir = [
    fullfile("-I" + osqp_dir, "include"); ...
    "-I" + qdldl_dir; ...
    fullfile("-I" + qdldl_dir, "qdldl_sources", "include")];


  %% OSQP Solver
  if any(strcmpi(command,'osqp')) || any(strcmpi(command,'all')) || (any(strcmpi(command,'osqp_mex')) && ~exist(fullfile(makefile_path, lib_name), 'file'))
    fprintf('Compiling OSQP solver...\n');

    % Create build directory and go inside
    if exist(osqp_build_dir, 'dir')
      rmdir(osqp_build_dir, 's');
    end
    mkdir(osqp_build_dir);

    % Restore the directory if there is a failure
    oldDirectory = cd(osqp_build_dir);
    try
      % Extend path for CMake mac (via Homebrew)
      PATH = getenv('PATH');
      if ((ismac) && (~contains(PATH, '/usr/local/bin')))
        setenv('PATH', [PATH ':/usr/local/bin']);
      end

      % Compile static library with CMake
      cmakeConfigure = sprintf('%s %s ..', 'cmake', cmake_args);
      if verboseFlag
        fprintf("\n%s\n", cmakeConfigure);
        [status, output] = system(cmakeConfigure,"-echo");
      else
        [status, output] = system(cmakeConfigure);
      end
      if(status)
        fprintf('\n');
        disp(output);
        error('Error configuring CMake environment');
      end

      cmakeBuild = sprintf('%s %s', make_cmd, '--target osqpstatic');
      if verboseFlag
        fprintf("\n%s\n", cmakeBuild);
        [status, output] = system(cmakeBuild,"-echo");
      else
        [status, output] = system(cmakeBuild);
      end
      if (status)
        fprintf('\n');
        disp(output);
        error('Error compiling OSQP');
      end

      % Copy static library to current folder
      lib_origin = fullfile(osqp_build_dir, 'out', lib_name);
      lib_destination = fullfile(makefile_path, lib_name);
      copyfile(lib_origin, lib_destination);

      cd(oldDirectory);
    catch e
      cd(oldDirectory);
      rethrow(e);
    end

    fprintf('Done\n\n\n');
  end

  %% osqpmex
  if( any(strcmpi(command,'osqp_mex')) || any(strcmpi(command,'all')) )
    % Compile interface
    fprintf('Compiling and linking osqpmex...\n');

    % Change directory back to matlab interface
    oldDirectory = cd(makefile_path);
    try

      % Compile
      osqp_mex_cpp = "osqp_mex.cpp";
      combined = cellfun(@(x) x, vertcat(mexoptflags, inc_dir, mex_libs, lib_name),"UniformOutput", false); % make into a cell so we can dump in variable arguments.
      mex(osqp_mex_cpp, combined{:});
      cd(oldDirectory);

      fprintf('Done\n\n\n');
    catch e
      cd(oldDirectory);
      rethrow(e);
    end
  end


  %% codegen
  if( any(strcmpi(command,'codegen')) || any(strcmpi(command,'all')) )
    fprintf('Copying source files for codegen...\n');

    % Copy C files
    cg_src_dir = fullfile(cg_sources_dir, 'src');
    if ~exist(cg_src_dir, 'dir')
      mkdir(cg_src_dir);
    end
    cdirs  = {fullfile(osqp_dir, 'src'),...
      fullfile(qdldl_dir),...
      fullfile(qdldl_dir, 'qdldl_sources', 'src')};
    for j = 1:length(cdirs)
      cfiles = dir(fullfile(cdirs{j},'*.c'));
      for i = 1 : length(cfiles)
        if ~any(strcmp(cfiles(i).name, {'cs.c', 'ctrlc.c', 'lin_sys.c', 'polish.c'}))
          copyfile(fullfile(cdirs{j}, cfiles(i).name), ...
            fullfile(cg_src_dir, cfiles(i).name));
        end
      end
    end

    % Copy H files
    cg_include_dir = fullfile(cg_sources_dir, 'include');
    if ~exist(cg_include_dir, 'dir')
      mkdir(cg_include_dir);
    end
    hdirs  = {fullfile(osqp_dir, 'include'),...
      fullfile(qdldl_dir),...
      fullfile(qdldl_dir, 'qdldl_sources', 'include')};
    for j = 1:length(hdirs)
      hfiles = dir(fullfile(hdirs{j},'*.h'));
      for i = 1 : length(hfiles)
        if ~any(strcmp(hfiles(i).name, {'qdldl_types.h', 'osqp_configure.h', ...
            'cs.h', 'ctrlc.h', 'lin_sys.h', 'polish.h'}))
          copyfile(fullfile(hdirs{j}, hfiles(i).name), ...
            fullfile(cg_include_dir, hfiles(i).name));
        end
      end
    end

    % Copy configure files
    cg_configure_dir = fullfile(cg_sources_dir, 'configure');
    if ~exist(cg_configure_dir, 'dir')
      mkdir(cg_configure_dir);
    end
    configure_dirs  = {fullfile(osqp_dir, 'configure'),...
      fullfile(qdldl_dir, 'qdldl_sources', 'configure')};
    for j = 1:length(configure_dirs)
      configure_files = dir(fullfile(configure_dirs{j},'*.h.in'));
      for i = 1 : length(configure_files)
        copyfile(fullfile(configure_dirs{j}, configure_files(i).name), ...
          fullfile(cg_configure_dir, configure_files(i).name));
      end
    end

    % Copy cmake files
    copyfile(fullfile(osqp_dir, 'src', 'CMakeLists.txt'), ...
      fullfile(cg_src_dir, 'CMakeLists.txt'));
    copyfile(fullfile(osqp_dir, 'include', 'CMakeLists.txt'), ...
      fullfile(cg_include_dir, 'CMakeLists.txt'));

    fprintf('Done\n\n\n');

  end


  %% clean
  if( any(strcmpi(command,'clean')) || any(strcmpi(command,'purge')) )
    fprintf('Cleaning mex files and library...\n');

    % Delete mex file
    mexfiles = dir(fullfile(makefile_path, '*.', mexext));
    for i = 1 : length(mexfiles)
      delete(mexfiles(i).name);
    end

    % Delete static library
    lib_full_path = fullfile(makefile_path, lib_name);
    if( exist(lib_full_path,'file') )
      delete(lib_full_path);
    end

    fprintf('Done\n\n\n');
  end


  %% purge
  if( any(strcmpi(command,'purge')) )
    fprintf('Cleaning OSQP build and codegen directories...\n');

    % Delete OSQP build directory
    if exist(osqp_build_dir, 'dir')
      rmdir(osqp_build_dir, 's');
    end

    % Delete codegen files
    if exist(cg_sources_dir, 'dir')
      rmdir(cg_sources_dir, 's');
    end

    fprintf('Done\n\n\n');
  end
end

%%
function checkCommandString(command)
  if any(contains(command, ["clean";"purge"],"IgnoreCase",true))
    if ~isscalar(command)
      error("When ""purge"" or ""clean"" is specified, no other commands can be specified.");
    end
  elseif any(contains(command,"all","IgnoreCase",true))
    if ~isscalar(command)
      error("When ""all"" is specified, no other commands can be specified.");
    end
  else
    if ~all(contains(command, ["osqp";"osqp_mex";"codegen"],"IgnoreCase",true))
      error("Command must be one or more of the following when ""all"", ""purge"", ""clean"" is not" + newline + ...
        "specified: ""osqp"", ""osqp_mex"", ""codegen"".\nYou entered: %s", strjoin(command, ', '));
    end
  end

end