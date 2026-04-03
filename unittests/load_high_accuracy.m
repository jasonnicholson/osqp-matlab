function [x_val, y_val, obj] = load_high_accuracy(test_name)
% LOAD_HIGH_ACCURACY Load high-accuracy reference solution for a test.
%
%   [x_val, y_val, obj] = load_high_accuracy('test_basic_QP')

    sol_dir = fullfile(fileparts(mfilename('fullpath')), 'solutions');
    data = load(fullfile(sol_dir, [test_name, '.mat']));
    x_val = data.x_val(:);
    y_val = data.y_val(:);
    obj = double(data.obj);
end
