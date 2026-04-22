% Repository usage summary
% Description: Summarize included/excluded subjects across SC and FC availability.
% Usage: Optional QC/reporting utility.
% Outputs: Useful for transparent sample accounting.
% Note: Update input paths, toolboxes, and filenames for your local environment.

% Public repository version: update file paths, toolboxes, and local settings before running.
% This script/function was lightly sanitized for sharing and may require project-specific inputs.

clear; clc;

%% paths
sc_root = '';

fc_roots = {
    'group1', '';
    'group2', '';
    'group3', '';
};

out_txt = '';
out_mat = '';

all_subjects = {};
with_sc = {};
with_fc = {};
matched = {};
excluded = {};

for g = 1:size(fc_roots,1)
    group_name = fc_roots{g,1};
    fc_root = fc_roots{g,2};

    %% all subjects from PROC/groupX/sub*
    sc_group_dir = fullfile(sc_root, group_name);
    subs = dir(fullfile(sc_group_dir, 'sub*'));
    subs = subs([subs.isdir]);

    all_ids = {};
    sc_ids = {};

    for i = 1:numel(subs)
        sid = subs(i).name;
        tag = [group_name '/' sid];
        all_ids{end+1,1} = sid; %#ok<SAGROW>
        all_subjects(end+1,1) = {tag}; %#ok<SAGROW>

        sc_file = fullfile(sc_group_dir, sid, 'tract', 'SC_counts_116x116.mat');
        if exist(sc_file, 'file')
            sc_ids{end+1,1} = sid; %#ok<SAGROW>
            with_sc(end+1,:) = {group_name, sid, sc_file}; %#ok<SAGROW>
        end
    end

    %% FC ids
    fc_txt = dir(fullfile(fc_root, 'ROICorrelation_FisherZ_sub*.txt'));
    fc_mat = dir(fullfile(fc_root, 'ROICorrelation_FisherZ_sub*.mat'));

    fc_map = struct();
    fc_ids = {};

    for i = 1:numel(fc_txt)
        name = fc_txt(i).name;
        tok = regexp(name, 'ROICorrelation_FisherZ_(sub\d+)\.txt$', 'tokens', 'once');
        if ~isempty(tok)
            sid = tok{1};
            if ~isfield(fc_map, sid)
                fc_map.(sid) = fullfile(fc_root, name);
                fc_ids{end+1,1} = sid; %#ok<SAGROW>
                with_fc(end+1,:) = {group_name, sid, fc_map.(sid)}; %#ok<SAGROW>
            end
        end
    end

    for i = 1:numel(fc_mat)
        name = fc_mat(i).name;
        tok = regexp(name, 'ROICorrelation_FisherZ_(sub\d+)\.mat$', 'tokens', 'once');
        if ~isempty(tok)
            sid = tok{1};
            if ~isfield(fc_map, sid)
                fc_map.(sid) = fullfile(fc_root, name);
                fc_ids{end+1,1} = sid; %#ok<SAGROW>
                with_fc(end+1,:) = {group_name, sid, fc_map.(sid)}; %#ok<SAGROW>
            end
        end
    end

    all_ids = unique(all_ids, 'stable');
    sc_ids = unique(sc_ids, 'stable');
    fc_ids = unique(fc_ids, 'stable');

    common_ids = intersect(sc_ids, fc_ids, 'stable');

    for i = 1:numel(common_ids)
        sid = common_ids{i};
        matched(end+1,:) = {group_name, sid}; %#ok<SAGROW>
    end

    %% excluded relative to total PROC subjects
    for i = 1:numel(all_ids)
        sid = all_ids{i};
        has_sc = ismember(sid, sc_ids);
        has_fc = ismember(sid, fc_ids);

        if ~(has_sc && has_fc)
            if ~has_sc && has_fc
                reason = 'missing_SC';
            elseif has_sc && ~has_fc
                reason = 'missing_FC';
            else
                reason = 'missing_SC_and_FC';
            end

            excluded(end+1,:) = {group_name, sid, reason}; %#ok<SAGROW>
        end
    end

    fprintf('\n%s\n', group_name);
    fprintf('  total in PROC = %d\n', numel(all_ids));
    fprintf('  with SC = %d\n', numel(sc_ids));
    fprintf('  with FC = %d\n', numel(fc_ids));
    fprintf('  matched SC+FC = %d\n', numel(common_ids));
    fprintf('  excluded = %d\n', numel(all_ids) - numel(common_ids));
end

%% write txt report
fid = fopen(out_txt, 'w');

fprintf(fid, '=== MATCH SUMMARY ===\n\n');

groups = unique(excluded(:,1));
for g = 1:numel(groups)
    group_name = groups{g};

    total_n = sum(startsWith(all_subjects, [group_name '/']));
    sc_n = sum(strcmp(with_sc(:,1), group_name));
    fc_n = sum(strcmp(with_fc(:,1), group_name));
    matched_n = sum(strcmp(matched(:,1), group_name));
    excl_n = sum(strcmp(excluded(:,1), group_name));

    fprintf(fid, '%s\n', group_name);
    fprintf(fid, 'total_in_PROC = %d\n', total_n);
    fprintf(fid, 'with_SC = %d\n', sc_n);
    fprintf(fid, 'with_FC = %d\n', fc_n);
    fprintf(fid, 'matched_SC_FC = %d\n', matched_n);
    fprintf(fid, 'excluded = %d\n\n', excl_n);
end

fprintf(fid, '=== EXCLUDED SUBJECTS ===\n');
fprintf(fid, 'group\tsub\treason\n');
for i = 1:size(excluded,1)
    fprintf(fid, '%s\t%s\t%s\n', excluded{i,1}, excluded{i,2}, excluded{i,3});
end

fclose(fid);

save(out_mat, 'all_subjects', 'with_sc', 'with_fc', 'matched', 'excluded');

disp(out_txt);
disp(out_mat);
