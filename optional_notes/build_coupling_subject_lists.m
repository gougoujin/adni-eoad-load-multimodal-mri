% Repository usage summary
% Description: Build subject lists for main and sensitivity coupling analyses.
% Usage: Optional cohort bookkeeping utility.
% Outputs: Useful for documenting flagged low-coverage subjects.
% Note: Update input paths, toolboxes, and filenames for your local environment.

% Public repository version: update file paths, toolboxes, and local settings before running.
% This script/function was lightly sanitized for sharing and may require project-specific inputs.

clear; clc;

inFile = '';

T = readtable(inFile, 'FileType', 'text', 'Delimiter', '\t');

% 手动低覆盖标记
flag_list = {
    'group1/sub005'
    'group1/sub039'
    'group1/sub040'
    'group2/sub053'
    'group2/sub054'
    'group3/sub003'
    'group3/sub063'
    'group3/sub090'
};

tags = strcat(string(T.group), "/", string(T.sub));

is_nan = isnan(T.mean_coupling) | T.valid_roi <= 0;
is_flag = ismember(cellstr(tags), flag_list);
is_low_valid = T.valid_roi < 50;

main_keep = ~is_nan & ~is_low_valid;   % 主分析
sens_flag = ~is_nan & (is_flag | is_low_valid); % 敏感性分析标记
drop = is_nan; % 直接排除

writetable(T(main_keep,:), '', ...
    'Delimiter', '\t', 'FileType', 'text');
writetable(T(sens_flag,:), '', ...
    'Delimiter', '\t', 'FileType', 'text');
writetable(T(drop,:), '', ...
    'Delimiter', '\t', 'FileType', 'text');

disp('done');
