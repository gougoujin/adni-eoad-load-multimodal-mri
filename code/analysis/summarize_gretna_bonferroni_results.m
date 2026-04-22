% Repository usage summary
% Description: Summarize Bonferroni-corrected GRETNA nodal results into master/detail tables.
% Usage: Run after GRETNA result folders have been generated.
% Outputs: Creates machine-readable summary tables used by downstream figures.
% Note: Update input paths, toolboxes, and filenames for your local environment.

% Public repository version: update file paths, toolboxes, and local settings before running.
% This script/function was lightly sanitized for sharing and may require project-specific inputs.

function summarize_gretna_bonferroni_results()
clear; clc;

%% ===================== 路径配置 =====================
root_dir   = '';
label_file = '';   % 按实际修改
out_dir    = fullfile(root_dir, '_summary_autoTF');

if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

master_file = fullfile(out_dir, 'gretna_bonferroni_master_summary_autoTF.txt');
detail_file = fullfile(out_dir, 'gretna_bonferroni_detail_summary_autoTF.txt');

labels = read_labels(label_file, 116);

dirs = dir(root_dir);
dirs = dirs([dirs.isdir]);
dirs = dirs(~ismember({dirs.name}, {'.','..','_summary_from_chatgpt', ...
    '_summary_from_chatgpt_v2','_summary_autoTF'}));

master_rows = {};
detail_rows = {};

fprintf('扫描目录: %s\n', root_dir);

for i = 1:numel(dirs)
    dname = dirs(i).name;
    dpath = fullfile(root_dir, dname);

    [metric_name, contrast_name] = parse_folder_name_loose(dname);

    fprintf('\n[%d/%d] %s\n', i, numel(dirs), dname);

    file_sets = detect_all_stat_files(dpath);

    if isempty(file_sets)
        fprintf('  [SKIP] 未检测到 T/F Vector 文件\n');
        continue;
    end

    for j = 1:numel(file_sets)
        stat_type = file_sets(j).stat_type;
        p_file    = file_sets(j).p_file;
        s_file    = file_sets(j).s_file;

        fprintf('  -> 读取 %s\n', stat_type);

        pvec = read_numeric_vector(p_file);
        svec = read_numeric_vector(s_file);

        if isempty(pvec)
            fprintf('     [SKIP] PVector 为空: %s\n', p_file);
            continue;
        end
        if isempty(svec)
            fprintf('     [SKIP] Statistic Vector 为空: %s\n', s_file);
            continue;
        end

        % 长度对齐
        n = min(numel(pvec), numel(svec));
        pvec = pvec(1:n);
        svec = svec(1:n);

        is_nodal = (n == 116);

        % 显著阈值
        sig_idx = find(isfinite(pvec) & pvec < 0.05);

        if isempty(sig_idx)
            n_sig = 0;
            top_hit = '';
            top_p = NaN;
            top_stat = NaN;
        else
            n_sig = numel(sig_idx);
            [~, ix] = min(pvec(sig_idx));
            top_idx = sig_idx(ix);

            if is_nodal
                top_hit = sprintf('%d:%s', top_idx, labels{top_idx});
            else
                top_hit = sprintf('Index%d', top_idx);
            end

            top_p = pvec(top_idx);
            top_stat = svec(top_idx);
        end

        % master
        master_rows(end+1,:) = { ...
            dname, stat_type, metric_name, contrast_name, ...
            n, is_nodal, n_sig, top_hit, top_p, top_stat, ...
            p_file, s_file}; %#ok<AGROW>

        % detail
        for k = 1:numel(sig_idx)
            idx = sig_idx(k);

            if is_nodal
                roi_id = idx;
                roi_label = labels{idx};
            else
                roi_id = idx;
                roi_label = 'GLOBAL_OR_NON116';
            end

            detail_rows(end+1,:) = { ...
                dname, stat_type, metric_name, contrast_name, ...
                roi_id, roi_label, pvec(idx), svec(idx), ...
                p_file, s_file}; %#ok<AGROW>
        end

        fprintf('     指标=%s | 比较=%s | 向量长度=%d | 显著数=%d\n', ...
            metric_name, contrast_name, n, n_sig);
    end
end

write_master(master_file, master_rows);
write_detail(detail_file, detail_rows);

fprintf('\n完成。\n');
fprintf('主汇总: %s\n', master_file);
fprintf('详细表: %s\n', detail_file);
end


%% ===================== 核心函数 =====================

function file_sets = detect_all_stat_files(dpath)
file_sets = struct('stat_type', {}, 'p_file', {}, 's_file', {});

% ---- T 检验 ----
t_p = fullfile(dpath, 'T2_PVector.txt');
t_s = fullfile(dpath, 'T2_TVector.txt');

if exist(t_p, 'file') && exist(t_s, 'file')
    file_sets(end+1).stat_type = 'T'; %#ok<AGROW>
    file_sets(end).p_file = t_p;
    file_sets(end).s_file = t_s;
end

% ---- F 检验 ----
f_p = fullfile(dpath, 'F_PVector.txt');
f_s = fullfile(dpath, 'F_FVector.txt');

if exist(f_p, 'file') && exist(f_s, 'file')
    file_sets(end+1).stat_type = 'F'; %#ok<AGROW>
    file_sets(end).p_file = f_p;
    file_sets(end).s_file = f_s;
end
end


function [metric_name, contrast_name] = parse_folder_name_loose(dname)
% 只解析“指标”和“比较”，不再解析 T/F
% 例子:
% 2BC
% 3BC
% 3BC-1and2
% 3NodalEfficiency-2and3
% 2NetworkEfficiency

metric_name = '';
contrast_name = '';

name_str = string(dname);

% 去掉开头数字前缀 2 / 3
name_str = regexprep(name_str, '^\d+', '');
parts = split(name_str, '-');

metric_raw = char(parts(1));
metric_name = normalize_metric_name(metric_raw);

if numel(parts) == 1
    % 没有写 pair 的情况：
    % 2... 默认按 group1_vs_group2
    % 3... 默认按 group_main_effect
    if startsWith(dname, '2')
        contrast_name = 'group1_vs_group2';
    else
        contrast_name = 'group_main_effect';
    end
else
    contrast_raw = char(parts(2));
    contrast_name = normalize_contrast_name(contrast_raw);
end
end


function name = normalize_metric_name(x)
x = char(string(x));
x_low = lower(strtrim(x));

switch x_low
    case 'bc'
        name = 'BetweennessCentrality';
    case 'dc'
        name = 'DegreeCentrality';
    case 'networkefficiency'
        name = 'NetworkEfficiency';
    case 'nodalefficiency'
        name = 'NodalEfficiency';
    case 'nodallocalefficiency'
        name = 'NodalLocalEfficiency';
    case 'nodalclustcoeff'
        name = 'NodalClustCoeff';
    case 'nodalshortestpath'
        name = 'NodalShortestPath';
    case 'smallworld'
        name = 'SmallWorld';
    case 'richclub'
        name = 'RichClub';
    case 'assortativity'
        name = 'Assortativity';
    case 'hierarchy'
        name = 'Hierarchy';
    case 'synchronization'
        name = 'Synchronization';
    otherwise
        name = x;
end
end


function name = normalize_contrast_name(x)
x = char(string(x));
x_low = lower(strtrim(x));

switch x_low
    case '1and2'
        name = 'group1_vs_group2';
    case '1and3'
        name = 'group1_vs_group3';
    case '2and3'
        name = 'group2_vs_group3';
    otherwise
        name = x;
end
end


function vec = read_numeric_vector(f)
vec = [];

try
    A = importdata(f);
catch
    return;
end

if isnumeric(A)
    vec = A(:);
elseif isstruct(A)
    if isfield(A, 'data') && isnumeric(A.data) && ~isempty(A.data)
        vec = A.data(:);
    elseif isfield(A, 'textdata') && ~isempty(A.textdata)
        tmp = str2double(string(A.textdata(:)));
        vec = tmp(:);
    else
        vec = [];
    end
else
    vec = [];
end

vec = vec(:);
vec = vec(isfinite(vec));
end


function labels = read_labels(label_file, nROI)
labels = cell(nROI,1);
for i = 1:nROI
    labels{i} = sprintf('ROI_%d', i);
end

if ~exist(label_file, 'file')
    warning('找不到 AAL 标签文件: %s', label_file);
    return;
end

fid = fopen(label_file, 'r');
if fid < 0
    warning('无法打开标签文件: %s', label_file);
    return;
end

k = 0;
while ~feof(fid) && k < nROI
    line = strtrim(fgetl(fid));
    if isempty(line)
        continue;
    end
    k = k + 1;
    labels{k} = line;
end
fclose(fid);
end


function write_master(outfile, rows)
fid = fopen(outfile, 'w');
if fid < 0
    error('无法写文件: %s', outfile);
end

fprintf(fid, ['Folder\tStatType\tMetric\tContrast\tVectorLength\tIsNodal116\t' ...
    'N_Significant\tTopHit\tTopHit_P\tTopHit_Stat\tPVectorFile\tStatVectorFile\n']);

for i = 1:size(rows,1)
    fprintf(fid, '%s\t%s\t%s\t%s\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t%s\n', ...
        tochar(rows{i,1}), ...
        tochar(rows{i,2}), ...
        tochar(rows{i,3}), ...
        tochar(rows{i,4}), ...
        rows{i,5}, ...
        rows{i,6}, ...
        rows{i,7}, ...
        tochar(rows{i,8}), ...
        num2str_nan(rows{i,9}), ...
        num2str_nan(rows{i,10}), ...
        tochar(rows{i,11}), ...
        tochar(rows{i,12}));
end

fclose(fid);
end


function write_detail(outfile, rows)
fid = fopen(outfile, 'w');
if fid < 0
    error('无法写文件: %s', outfile);
end

fprintf(fid, 'Folder\tStatType\tMetric\tContrast\tROI\tLabel\tP\tStat\tPVectorFile\tStatVectorFile\n');

for i = 1:size(rows,1)
    fprintf(fid, '%s\t%s\t%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s\n', ...
        tochar(rows{i,1}), ...
        tochar(rows{i,2}), ...
        tochar(rows{i,3}), ...
        tochar(rows{i,4}), ...
        rows{i,5}, ...
        tochar(rows{i,6}), ...
        num2str_nan(rows{i,7}), ...
        num2str_nan(rows{i,8}), ...
        tochar(rows{i,9}), ...
        tochar(rows{i,10}));
end

fclose(fid);
end


function s = num2str_nan(x)
if isnan(x)
    s = 'NaN';
else
    s = sprintf('%.6g', x);
end
end


function s = tochar(x)
s = char(string(x));
end
