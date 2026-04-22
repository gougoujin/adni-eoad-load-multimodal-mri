% Repository usage summary
% Description: Match subjects with available structural and functional connectivity files.
% Usage: matched = match_structural_and_functional_subjects(...) or run as a script after setting paths.
% Outputs: Creates matched subject lists and file paths for downstream coupling/SDI analyses.
% Note: Update input paths, toolboxes, and filenames for your local environment.

% Public repository version: update file paths, toolboxes, and local settings before running.
% This script/function was lightly sanitized for sharing and may require project-specific inputs.

clear; clc;

% ---------- 路径 ----------
sc_root = '';

fc_roots = {
    'group1', '';
    'group2', '';
    'group3', '';
};

out_txt = '';
out_mat = '';

matched = {};
only_sc = {};
only_fc = {};

for g = 1:size(fc_roots,1)
    group_name = fc_roots{g,1};
    fc_root = fc_roots{g,2};

    % ----- SC 被试 -----
    sc_group_dir = fullfile(sc_root, group_name);
    sc_subs = dir(fullfile(sc_group_dir, 'sub*'));
    sc_subs = sc_subs([sc_subs.isdir]);

    sc_ids = {};
    sc_map = struct();

    for i = 1:numel(sc_subs)
        sid = sc_subs(i).name;   % sub001
        sc_file = fullfile(sc_group_dir, sid, 'tract', 'SC_counts_116x116.mat');
        if exist(sc_file, 'file')
            sc_ids{end+1,1} = sid; %#ok<SAGROW>
            sc_map.(sid) = sc_file;
        end
    end

    % ----- FC 被试 -----
    fc_txt = dir(fullfile(fc_root, 'ROICorrelation_FisherZ_sub*.txt'));
    fc_mat = dir(fullfile(fc_root, 'ROICorrelation_FisherZ_sub*.mat'));

    fc_ids = {};
    fc_map = struct();

    for i = 1:numel(fc_txt)
        name = fc_txt(i).name;
        tok = regexp(name, 'ROICorrelation_FisherZ_(sub\d+)\.txt$', 'tokens', 'once');
        if ~isempty(tok)
            sid = tok{1};
            fc_ids{end+1,1} = sid; %#ok<SAGROW>
            fc_map.(sid) = fullfile(fc_root, name);
        end
    end

    for i = 1:numel(fc_mat)
        name = fc_mat(i).name;
        tok = regexp(name, 'ROICorrelation_FisherZ_(sub\d+)\.mat$', 'tokens', 'once');
        if ~isempty(tok)
            sid = tok{1};
            if ~isfield(fc_map, sid) % 优先 txt，没 txt 才用 mat
                fc_ids{end+1,1} = sid; %#ok<SAGROW>
                fc_map.(sid) = fullfile(fc_root, name);
            end
        end
    end

    sc_ids = unique(sc_ids, 'stable');
    fc_ids = unique(fc_ids, 'stable');

    common_ids = intersect(sc_ids, fc_ids, 'stable');
    only_sc_ids = setdiff(sc_ids, fc_ids, 'stable');
    only_fc_ids = setdiff(fc_ids, sc_ids, 'stable');

    fprintf('\n%/path/to/project', group_name);
    fprintf('  SC 被试数: %d\n', numel(sc_ids));
    fprintf('  FC 被试数: %d\n', numel(fc_ids));
    fprintf('  共同被试数: %d\n', numel(common_ids));
    fprintf('  仅 SC: %d\n', numel(only_sc_ids));
    fprintf('  仅 FC: %d\n', numel(only_fc_ids));

    for i = 1:numel(common_ids)
        sid = common_ids{i};
        matched(end+1,:) = {group_name, sid, sc_map.(sid), fc_map.(sid)}; %#ok<SAGROW>
    end

    for i = 1:numel(only_sc_ids)
        sid = only_sc_ids{i};
        only_sc(end+1,:) = {group_name, sid, sc_map.(sid)}; %#ok<SAGROW>
    end

    for i = 1:numel(only_fc_ids)
        sid = only_fc_ids{i};
        only_fc(end+1,:) = {group_name, sid, fc_map.(sid)}; %#ok<SAGROW>
    end
end

% ----- 导出 txt -----
fid = fopen(out_txt, 'w');
fprintf(fid, 'group\tsub\tsc_file\tfc_file\n');
for i = 1:size(matched,1)
    fprintf(fid, '%s\t%s\t%s\t%s\n', matched{i,1}, matched{i,2}, matched{i,3}, matched{i,4});
end
fclose(fid);

save(out_mat, 'matched', 'only_sc', 'only_fc');

fprintf('\n输出：\n%s\n%s\n', out_txt, out_mat);
