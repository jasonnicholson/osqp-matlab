function make_osqp(varargin)
% MAKE_OSQP  Build the OSQP v1.0.0 MATLAB MEX interface.
%
%   MAKE_OSQP              Build everything (osqp library + mex) with default
%                          algebra backend (builtin).
%
%   MAKE_OSQP('algebra','builtin')   Build with builtin (QDLDL) backend.
%   MAKE_OSQP('algebra','mkl')       Build with Intel MKL backend.
%   MAKE_OSQP('algebra','cuda')      Build with CUDA backend.
%
%   MAKE_OSQP('-verbose')            Print full compiler output.
%
%   Targets (last argument):
%     'all'       — build osqp library + mex (default)
%     'osqp'      — build the OSQP C library only via CMake
%     'osqp_mex'  — compile the MEX interface only (library must exist)
%     'clean'     — remove build artefacts
%     'purge'     — clean + remove the CMake build directory

% ---------------------------------------------------------------
% Parse arguments
% ---------------------------------------------------------------
verbose = false;
algebra = 'builtin';
what    = 'all';

i = 1;
while i <= nargin
    arg = varargin{i};
    if strcmpi(arg, '-verbose')
        verbose = true;
    elseif strcmpi(arg, 'algebra')
        i = i + 1;
        algebra = lower(varargin{i});
        assert(ismember(algebra, {'builtin','mkl','cuda'}), ...
            'Algebra backend must be builtin, mkl, or cuda.');
    elseif any(strcmpi(arg, {'all','osqp','osqp_mex','clean','purge'}))
        what = lower(arg);
    else
        fprintf('Ignoring unrecognised argument "%s".\n', arg);
    end
    i = i + 1;
end

% ---------------------------------------------------------------
% Unlock any existing MEX
% ---------------------------------------------------------------
if mislocked('osqp_mex')
    munlock('osqp_mex');
end

% ---------------------------------------------------------------
% Paths
% ---------------------------------------------------------------
[makefile_path, ~, ~] = fileparts(which('make_osqp.m'));
osqp_root = fileparts(makefile_path);  % go up from src/ to project root
build_dir = fullfile(osqp_root, 'build');

% ---------------------------------------------------------------
% Compiler / linker helpers
% ---------------------------------------------------------------
mex_cmd = 'mex -O -silent';
mex_libs = '';

% Ctrl-C support
if ispc
    ut = fullfile(matlabroot, 'extern', 'lib', computer('arch'), ...
                  'mingw64', 'libut.lib');
    mex_libs = sprintf('%s "%s"', mex_libs, ut);
else
    mex_libs = sprintf('%s %s', mex_libs, '-lut');
end
% Shared library loading (Linux)
if isunix && ~ismac
    mex_libs = sprintf('%s %s', mex_libs, '-ldl');
end

% MEX flags
mexoptflags = '-DMATLAB';
if ~isempty(strfind(computer, '64')) && verLessThan('matlab', '9.4') %#ok<STREMP>
    mexoptflags = sprintf('%s %s', mexoptflags, '-largeArrayDims');
end
if ~verLessThan('matlab', '9.4')
    mexoptflags = sprintf('%s %s', mexoptflags, '-R2017b');
end
if ~ispc
    mexoptflags = sprintf('%s %s', mexoptflags, 'COPTIMFLAGS=''-O3''');
end

% ---------------------------------------------------------------
% CMake generator
% ---------------------------------------------------------------
if ispc
    cmake_generator = '-G "MinGW Makefiles"';
else
    cmake_generator = '-G "Unix Makefiles"';
end

% Matlab root for CMake
Matlab_ROOT = strrep(matlabroot, '\', '/');

% ---------------------------------------------------------------
% Build OSQP library via CMake + FetchContent
% ---------------------------------------------------------------
if any(strcmpi(what, {'osqp','all'}))
    fprintf('Building OSQP v1.0.0 (%s backend)...', algebra);

    if exist(build_dir, 'dir')
        rmdir(build_dir, 's');
    end
    mkdir(build_dir);

    % Extend PATH for Homebrew on macOS
    PATH = getenv('PATH');
    if ismac && isempty(strfind(PATH, '/usr/local/bin')) %#ok<STREMP>
        setenv('PATH', [PATH ':/usr/local/bin']);
    end

    cmake_args = sprintf([ ...
        '%s ' ...                               % generator
        '-DOSQP_ALGEBRA_BACKEND=%s ' ...
        '-DOSQP_BUILD_SHARED_LIB=OFF ' ...
        '-DOSQP_BUILD_STATIC_LIB=ON ' ...
        '-DOSQP_BUILD_DEMO_EXE=OFF ' ...
        '-DOSQP_BUILD_UNITTESTS=OFF ' ...
        '-DOSQP_ENABLE_PRINTING=ON ' ...
        '-DOSQP_ENABLE_PROFILING=ON ' ...
        '-DOSQP_ENABLE_INTERRUPT=ON ' ...
        '-DMatlab_ROOT_DIR="%s"'], ...
        cmake_generator, algebra, Matlab_ROOT);

    % Configure
    orig_dir = pwd;
    cd(build_dir);
    [status, output] = system(sprintf('cmake %s "%s"', cmake_args, osqp_root));
    if status
        cd(orig_dir);
        fprintf('\n'); disp(output);
        error('CMake configuration failed.');
    elseif verbose
        fprintf('\n'); disp(output);
    end

    % Build
    [status, output] = system('cmake --build . --target osqpstatic --config Release');
    if status
        cd(orig_dir);
        fprintf('\n'); disp(output);
        error('OSQP build failed.');
    elseif verbose
        fprintf('\n'); disp(output);
    end

    % Build codegen source staging target
    [status, output] = system('cmake --build . --target copy_codegen_files --config Release');
    if status && verbose
        fprintf('\n(codegen staging skipped)\n');
    elseif verbose
        fprintf('\n'); disp(output);
    end

    cd(orig_dir);
    fprintf('\t[done]\n');
end

% ---------------------------------------------------------------
% Compile MEX interface
% ---------------------------------------------------------------
if any(strcmpi(what, {'osqp_mex','all'}))
    fprintf('Compiling osqp_mex...');

    % Discover include dirs from the CMake build
    % The FetchContent checkout lands under build/_deps/osqp-src
    osqp_src_dir = fullfile(build_dir, '_deps', 'osqp-src');
    osqp_bin_dir = fullfile(build_dir, '_deps', 'osqp-build');

    inc_dirs = { ...
        fullfile(osqp_src_dir, 'include', 'public'), ...
        fullfile(osqp_src_dir, 'include', 'private'), ...
        fullfile(osqp_bin_dir, 'include', 'public') ...
    };

    % Also find the algebra-specific include dir
    algebra_inc = fullfile(osqp_src_dir, 'algebra', 'builtin', 'include');
    if exist(algebra_inc, 'dir')
        inc_dirs{end+1} = algebra_inc;
    end

    inc_flags = '';
    for k = 1:length(inc_dirs)
        inc_flags = sprintf('%s -I"%s"', inc_flags, inc_dirs{k});
    end

    % Find the static library
    lib_name = findlib(build_dir);

    % Build MEX
    cmd = sprintf('%s %s %s -outdir "%s" "%s" "%s" %s', ...
        mex_cmd, mexoptflags, inc_flags, ...
        makefile_path, ...
        fullfile(makefile_path, 'osqp_mex.cpp'), ...
        lib_name, mex_libs);

    eval(cmd);
    fprintf('\t\t\t\t[done]\n');
end

% ---------------------------------------------------------------
% Clean / Purge
% ---------------------------------------------------------------
if any(strcmpi(what, {'clean','purge'}))
    fprintf('Cleaning...');
    mexfiles = dir(fullfile(makefile_path, ['*.', mexext]));
    for k = 1:length(mexfiles)
        delete(fullfile(mexfiles(k).folder, mexfiles(k).name));
    end
    fprintf('\t\t\t\t\t[done]\n');
end

if strcmpi(what, 'purge')
    fprintf('Purging build directory...');
    if exist(build_dir, 'dir')
        rmdir(build_dir, 's');
    end
    fprintf('\t\t\t[done]\n');
end

end

% ---------------------------------------------------------------
% Helper: recursively find libosqp.a / osqp.lib in the build tree
% ---------------------------------------------------------------
function libpath = findlib(build_dir)
    if ispc
        patterns = {'osqpstatic.lib', 'osqp.lib'};
    else
        patterns = {'libosqpstatic.a', 'libosqp.a'};
    end
    for p = 1:numel(patterns)
        f = dir(fullfile(build_dir, '**', patterns{p}));
        if ~isempty(f)
            libpath = fullfile(f(1).folder, f(1).name);
            return;
        end
    end
    error('Cannot find OSQP static library under %s. Run make_osqp(''osqp'') first.', ...
          build_dir);
end
