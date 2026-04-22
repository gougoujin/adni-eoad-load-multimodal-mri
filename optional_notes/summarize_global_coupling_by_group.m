% Repository usage summary
% Description: Print simple descriptive summaries for global coupling by group.
% Usage: Optional descriptive utility.
% Outputs: Useful for quick checks, not required for the full pipeline.
% Note: Update input paths, toolboxes, and filenames for your local environment.

% Public repository version: update file paths, toolboxes, and local settings before running.
% This script/function was lightly sanitized for sharing and may require project-specific inputs.

clear; clc;

T = readtable('', 'FileType', 'text', 'Delimiter', '\t');

% 先排除完全无效
T = T(~isnan(T.mean_coupling) & T.valid_roi > 0, :);

groups = unique(T.group);

for i = 1:numel(groups)
    idx = strcmp(T.group, groups{i});
    x = T.mean_coupling(idx);

    fprintf('\n%s\n', groups{i});
    fprintf('n = %d\n', numel(x));
    fprintf('mean = %.4f\n', mean(x));
    fprintf('std = %.4f\n', std(x));
    fprintf('median = %.4f\n', median(x));
    fprintf('min = %.4f\n', min(x));
    fprintf('max = %.4f\n', max(x));
end

% 简单单因素 ANOVA
[p, tbl, stats] = anova1(T.mean_coupling, T.group, 'off');
fprintf('\nANOVA p = %.6g\n', p);

% 两两比较
results = multcompare(stats, 'Display', 'off');
disp(results);
