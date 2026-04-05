function make_osqp(options)
% MAKE_OSQP  Build the OSQP MATLAB MEX interface.
%
%   make_osqp("Name", "value")
%   make_osqp(Name=value)
%
% Name Value Pairs:
%   algebra ((1,1) string, default="builtin"): Linear algebra backend for KKT factorization. Options:
%     - "builtin"  - QDLDL backend.
%     - "mkl"      - Intel MKL pardiso. May require extra modification and thus requires expertise.
%     - "cuda"     - CUDA backend.
%
%   verbose ((1,1) logical, default=false): Print full compiler output.
%
%   target ((1,1) string, default="all"): Build target. Options:
%     - "all"      - Build OSQP library + MEX interface (default).
%     - "osqp"     - Build the OSQP C library only via CMake.
%     - "osqp_mex" - Compile the MEX interface only (library must already exist).
%     - "clean"    - Remove build artefacts.
%     - "purge"    - Clean + remove the CMake build directory.

arguments
    options.algebra (1,1) string {mustBeMember(options.algebra, ["builtin","mkl","cuda"])} = "builtin"
    options.verbose (1,1) logical = false
    options.target  (1,1) string {mustBeMember(options.target,  ["all","osqp","osqp_mex","clean","purge"])} = "all"
end

algebra = options.algebra;
verbose = options.verbose;
what    = options.target;

% ---------------------------------------------------------------
% Unlock any existing MEX
% ---------------------------------------------------------------
if mislocked("osqp_mex")
    munlock("osqp_mex");
end

% ---------------------------------------------------------------
% Paths
% ---------------------------------------------------------------
[makefile_path, ~, ~] = fileparts(which("make_osqp.m"));
osqp_root = fileparts(makefile_path);  % go up from src/ to project root
build_dir = fullfile(osqp_root, "build");

% ---------------------------------------------------------------
% Compiler / linker helpers
% ---------------------------------------------------------------
% MEX compile flags
mex_flags = ["-O", "-silent", "-DMATLAB"];
if contains(computer, "64") && verLessThan("matlab", "9.4")
    mex_flags(end+1) = "-largeArrayDims";
end
if ~verLessThan("matlab", "9.4")
    mex_flags(end+1) = "-R2017b";
end
if ~ispc
    mex_flags(end+1) = "COPTIMFLAGS='-O3'";
end

% Link libraries
mex_link = strings(1,0);
if ispc
    ut = fullfile(matlabroot, "extern", "lib", computer("arch"), "mingw64", "libut.lib");
    mex_link(end+1) = string(ut);
else
    mex_link(end+1) = "-lut";
end
if isunix && ~ismac
    mex_link(end+1) = "-ldl";
end

% ---------------------------------------------------------------
% CMake generator
% ---------------------------------------------------------------
if ispc
    cmake_generator = "-G ""MinGW Makefiles""";
else
    cmake_generator = "-G ""Unix Makefiles""";
end

% Matlab root for CMake
Matlab_ROOT = replace(matlabroot, "\", "/");

% ---------------------------------------------------------------
% Build OSQP library via CMake + FetchContent
% ---------------------------------------------------------------
if any(what == ["osqp", "all"])
    fprintf("Building OSQP (%s backend)...", algebra);

    if exist(build_dir, "dir")
        rmdir(build_dir, "s");
    end
    mkdir(build_dir);

    % Extend PATH for Homebrew on macOS
    PATH = string(getenv("PATH"));
    if ismac && ~contains(PATH, "/usr/local/bin")
        setenv("PATH", PATH + ":/usr/local/bin");
    end

    % On Linux, remove MATLAB's bundled libraries from LD_LIBRARY_PATH
    % to avoid version conflicts with cmake and system tools
    if isunix && ~ismac
        LD_LIBRARY_PATH = string(getenv("LD_LIBRARY_PATH"));
        libraryPaths = strsplit(LD_LIBRARY_PATH, ":");
        matlab_sys_lib = string(fullfile(matlabroot, "sys/os/glnxa64"));
        indexMatch = matches(libraryPaths, matlab_sys_lib);
        if any(indexMatch)
            % Remove MATLAB's sys/os/glnxa64 path from LD_LIBRARY_PATH
            libraryPaths(indexMatch) = [];
            LD_LIBRARY_PATH_NEW = strjoin(libraryPaths, ":");
            setenv("LD_LIBRARY_PATH", LD_LIBRARY_PATH_NEW);

            % Reset to original LD_LIBRARY_PATH
            reset_LD_LIBRARY_PATH = onCleanup(@() setenv("LD_LIBRARY_PATH", LD_LIBRARY_PATH));
        end
    end

    cmake_args = sprintf( ...
        "%s "                               + ... % generator
        "-DOSQP_ALGEBRA_BACKEND=%s "        + ...
        "-DOSQP_BUILD_SHARED_LIB=OFF "      + ...
        "-DOSQP_BUILD_STATIC_LIB=ON "       + ...
        "-DOSQP_BUILD_DEMO_EXE=OFF "        + ...
        "-DOSQP_BUILD_UNITTESTS=OFF "       + ...
        "-DOSQP_ENABLE_PRINTING=ON "        + ...
        "-DOSQP_ENABLE_PROFILING=ON "       + ...
        "-DOSQP_ENABLE_INTERRUPT=ON "       + ...
        "-DMatlab_ROOT_DIR=""%s""", ...
        cmake_generator, algebra, Matlab_ROOT);

    % Configure
    orig_dir = pwd;
    cleanup_dir = onCleanup(@() cd(orig_dir));
    cd(build_dir);
    [status, output] = system(sprintf("cmake %s ""%s""", cmake_args, makefile_path));
    if status
        fprintf("\n"); disp(output);
        error("CMake configuration failed.");
    elseif verbose
        fprintf("\n"); disp(output);
    end

    % Build
    [status, output] = system("cmake --build . --target osqpstatic --config Release");
    if status
        fprintf("\n"); disp(output);
        error("OSQP build failed.");
    elseif verbose
        fprintf("\n"); disp(output);
    end

    % Build codegen source staging target
    [status, output] = system("cmake --build . --target copy_codegen_files --config Release");
    if status && verbose
        fprintf("\n(codegen staging skipped)\n");
    elseif verbose
        fprintf("\n"); disp(output);
    end
    fprintf("\t[done]\n");
end

% ---------------------------------------------------------------
% Compile MEX interface
% ---------------------------------------------------------------
if any(what == ["osqp_mex", "all"])
    fprintf("Compiling osqp_mex...");

    % Discover include dirs from the CMake build
    % The FetchContent checkout lands under build/_deps/osqp-src
    osqp_src_dir = fullfile(build_dir, "_deps", "osqp-src");
    osqp_bin_dir = fullfile(build_dir, "_deps", "osqp-build");

    inc_dirs = [ ...
        string(fullfile(osqp_src_dir, "include", "public")), ...
        string(fullfile(osqp_src_dir, "include", "private")), ...
        string(fullfile(osqp_bin_dir, "include", "public")) ...
    ];

    % Also find the algebra-specific include dir
    algebra_inc = fullfile(osqp_src_dir, "algebra", "builtin", "include");
    if exist(algebra_inc, "dir")
        inc_dirs(end+1) = string(algebra_inc);
    end

    mex_inc = "-I" + inc_dirs;

    % Find the static library
    lib_name = findlib(build_dir);

    % Compile MEX
    mex_flag_args = cellstr(mex_flags);
    mex_inc_args = cellstr(mex_inc);
    mex_link_args = cellstr(mex_link);
    mex(mex_flag_args{:}, mex_inc_args{:}, "-outdir", makefile_path, ...
        fullfile(makefile_path, "osqp_mex.cpp"), lib_name, mex_link_args{:});
    fprintf("\t\t\t\t[done]\n");
end

% ---------------------------------------------------------------
% Clean / Purge
% ---------------------------------------------------------------
if any(what == ["clean", "purge"])
    fprintf("Cleaning...");
    mexfiles = dir(fullfile(makefile_path, "*." + mexext));
    for k = 1:length(mexfiles)
        delete(fullfile(mexfiles(k).folder, mexfiles(k).name));
    end
    fprintf("\t\t\t\t\t[done]\n");
end

if what == "purge"
    fprintf("Purging build directory...");
    if exist(build_dir, "dir")
        rmdir(build_dir, "s");
    end
    fprintf("\t\t\t[done]\n");
end

end

% ---------------------------------------------------------------
% Helper: recursively find libosqp.a / osqp.lib in the build tree
% ---------------------------------------------------------------
function libpath = findlib(build_dir)
    if ispc
        patterns = ["osqpstatic.lib", "osqp.lib"];
    else
        patterns = ["libosqpstatic.a", "libosqp.a"];
    end
    for p = 1:numel(patterns)
        f = dir(fullfile(build_dir, "**", patterns(p)));
        if ~isempty(f)
            libpath = fullfile(f(1).folder, f(1).name);
            return;
        end
    end
    error("Cannot find OSQP static library under %s. Run make_osqp('osqp') first.", ...
          build_dir);
end
