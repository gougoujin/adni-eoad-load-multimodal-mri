% Repository usage summary
% Description: Batch-build structural connectivity matrices for all subjects from probtrackx outputs.
% Usage: Run after setting the cohort root directory and ROI count.
% Outputs: Calls the subject-level builder across all group/sub directories.
% Note: Update input paths, toolboxes, and filenames for your local environment.

% Public repository version: update file paths, toolboxes, and local settings before running.
% This script/function was lightly sanitized for sharing and may require project-specific inputs.

clear; clc;
% addpath('/path/to/required/toolbox');

root_dir = '';
nROI = 116;

groups = dir(fullfile(root_dir, 'group*'));
groups = groups([groups.isdir]);

for g = 1:numel(groups)
    gdir = fullfile(root_dir, groups(g).name);
    subs = dir(fullfile(gdir, 'sub*'));
    subs = subs([subs.isdir]);

    for s = 1:numel(subs)
        subj_dir = fullfile(gdir, subs(s).name);

        fprintf('\n============================\n');
        fprintf('处理 %s\\%s\n', groups(g).name, subs(s).name);
        fprintf('============================\n');

        try
            build_subject_structural_connectivity_from_probtrackx(subj_dir, nROI);
        catch ME
            warning('失败: %s\\%s -> %s', groups(g).name, subs(s).name, ME.message);
        end
    end
end

fprintf('\n全部完成。\n');
