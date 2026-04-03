function install_osqp
    % Install the OSQP solver MATLAB interface
    %
    % This function clones osqp-matlab, builds the MEX interface using
    % CMake FetchContent (automatically downloads OSQP v1.0.0), and adds
    % the installation directory to the MATLAB path.

    fprintf('Cloning osqp-matlab...');
    if exist('osqp-matlab', 'dir')
        rmdir('osqp-matlab', 's');
    end
    system('git clone https://github.com/osqp/osqp-matlab.git');
    fprintf('\t\t\t[done]\n');

    fprintf('Building MEX interface...\n');
    cd('osqp-matlab');
    make_osqp;
    fprintf('Building MEX interface...\t\t[done]\n');

    fprintf('Updating path...');
    addpath(pwd);
    savepath;
    fprintf('\t\t\t\t[done]\n');

    cd('..');
    fprintf('OSQP v%s is successfully installed!\n', osqp.version());
end
