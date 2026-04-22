% Repository usage summary
% Description: Run multimodal clinical association models integrating GRETNA and functional metrics.
% Usage: Run after the relevant coupling, GRETNA, and ALFF/fALFF inputs are prepared.
% Outputs: Outputs selected imaging-clinical association summaries.
% Note: Update input paths, toolboxes, and filenames for your local environment.

% Public repository version: update file paths, toolboxes, and local settings before running.
% This script/function was lightly sanitized for sharing and may require project-specific inputs.

function run_multimodal_clinical_association_models()
clear; clc;

%% ===================== 用户配置 =====================
cov_file = '';

aal_label_file = '';
aal_atlas_file = '';   % 改成你的实际 AAL116 图谱

gretna_root = '';

dpabi_group_roots = { ...
    '', ...   % group1 = EOAD
    '', ...  % group2 = LOAD
    '' ...      % group3 = HC
    };

out_dir = '';
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

%% ===================== 预设：只跑最值得看的指标 =====================
% GRETNA 节点指标及对应文件名
gretna_specs = { ...
    'BetweennessCentrality', 'aBc.txt'; ...
    'DegreeCentrality',      'aDc.txt'; ...
    'NodalEfficiency',       'aNe.txt'; ...
    'NodalLocalEfficiency',  'aNLe.txt' ...
    };

% GRETNA 代表性 ROI
gretna_roi_labels = { ...
    'Vermis_9', ...
    'Heschl_L', ...
    'Temporal_Sup_L', ...
    'Cerebelum_4_5_L', ...
    'Frontal_Inf_Oper_L', ...
    'Frontal_Sup_Medial_L' ...
    };

% ALFF 代表性 ROI
alff_roi_labels = { ...
    'Hippocampus_R', ...
    'Hippocampus_L', ...
    'Thalamus_R', ...
    'Caudate_R', ...
    'Pallidum_L', ...
    'Insula_L' ...
    };

% fALFF 代表性 ROI
falff_roi_labels = { ...
    'Cuneus_L', ...
    'Cerebelum_Crus1_L', ...
    'Occipital_Inf_R', ...
    'Temporal_Mid_L', ...
    'Temporal_Pole_Sup_R', ...
    'Precuneus_L' ...
    };

%% ===================== 读 AAL 标签 =====================
aal_labels = read_labels(aal_label_file, 116);

gretna_roi_idx = labels_to_idx(gretna_roi_labels, aal_labels);
alff_roi_idx   = labels_to_idx(alff_roi_labels, aal_labels);
falff_roi_idx  = labels_to_idx(falff_roi_labels, aal_labels);

%% ===================== 读取临床表 =====================
% 第1列=group, 第3列=sub, 第4列=age, 第5列=sex, 第6:11列=量表
raw = readcell(cov_file);
header = raw(1,:);
raw = raw(2:end,:);

n = size(raw,1);

group_num = nan(n,1);
sub_id    = cell(n,1);
age       = nan(n,1);
sex       = nan(n,1);

for i = 1:n
    group_num(i) = to_num(raw{i,1});

    s = raw{i,3};
    if isnumeric(s) && ~isnan(s)
        sub_id{i} = sprintf('sub%03d', round(s));
    else
        tok = regexp(char(string(s)), '(\d+)', 'tokens', 'once');
        if ~isempty(tok)
            sub_id{i} = sprintf('sub%03d', str2double(tok{1}));
        else
            sub_id{i} = '';
        end
    end

    age(i) = to_num(raw{i,4});
    sex(i) = to_num(raw{i,5});
end

scale_cols = 6:11; % F:K
nScale = numel(scale_cols);
scale_names = cell(nScale,1);
scale_mat = nan(n, nScale);

for j = 1:nScale
    scale_names{j} = safe_var_name(header{scale_cols(j)}, sprintf('Scale%d', j));
    for i = 1:n
        scale_mat(i,j) = to_num(raw{i, scale_cols(j)});
    end
end

T = table;
T.group_num = group_num;
T.sub = string(sub_id);
T.age = age;
T.sex = sex;
for j = 1:nScale
    T.(scale_names{j}) = scale_mat(:,j);
end
T.tag = string(T.group_num) + "_" + T.sub;
T = T(isfinite(T.group_num) & strlength(T.sub)>0 & isfinite(T.age) & isfinite(T.sex), :);
T_clinical = T;

%% ===================== 读取 GRETNA（按 aBc/aDc/aNe/aNLe） =====================
fprintf('\n========== 读取 GRETNA ==========\n');

for m = 1:size(gretna_specs,1)
    metric_dir_name = gretna_specs{m,1};
    metric_file_name = gretna_specs{m,2};

    fprintf('GRETNA metric: %s (%s)\n', metric_dir_name, metric_file_name);

    G = table();

    for g = 1:3
        gfile = fullfile(gretna_root, metric_dir_name, sprintf('group%d', g), metric_file_name);
        if ~exist(gfile, 'file')
            fprintf('  [WARN] 文件不存在: %s\n', gfile);
            continue;
        end

        A = readmatrix(gfile);
        A = double(A);

        % 去掉全空行/列
        if isempty(A)
            fprintf('  [WARN] 空文件: %s\n', gfile);
            continue;
        end

        % 这一组的被试列表：从临床表里取，再按 sub 编号排序
        subs_g = T_clinical.sub(T_clinical.group_num == g);
        subs_g = sort_sub_ids(unique(subs_g));
        nsub = numel(subs_g);

        sa = size(A);

        % 自动判断方向
        % 目标：B = Nsub × 116
        if sa(2) == 116
            B = A;
        elseif sa(1) == 116
            B = A';
        else
            fprintf('  [WARN] 无法识别矩阵方向: %s size=[%d %d]\n', gfile, sa(1), sa(2));
            continue;
        end

        % 裁剪到共同最小长度
        n_use = min(size(B,1), nsub);
        B = B(1:n_use, :);
        subs_use = subs_g(1:n_use);

        X = table();
        X.group_num = repmat(g, n_use, 1);
        X.sub = subs_use(:);
        X.tag = string(X.group_num) + "_" + X.sub;

        for r = 1:numel(gretna_roi_idx)
            roi = gretna_roi_idx(r);
            if isnan(roi), continue; end
            vname = safe_var_name([metric_dir_name '_' aal_labels{roi}], sprintf('%s_ROI%d', metric_dir_name, roi));
            X.(vname) = B(:, roi);
        end

        if isempty(G)
            G = X;
        else
            G = [G; X]; %#ok<AGROW>
        end
    end

    if ~isempty(G)
        T = outerjoin(T, G, 'Keys', 'tag', 'MergeKeys', true);
    end
end

%% ===================== 提取 zALFF / zfALFF =====================
fprintf('\n========== 提取 zALFF / zfALFF ==========\n');

atlas = double(read_nifti_auto(aal_atlas_file));

for g = 1:3
    group_root = dpabi_group_roots{g};
    if ~exist(group_root, 'dir')
        fprintf('[WARN] 目录不存在: %s\n', group_root);
        continue;
    end

    fprintf('扫描 group%d: %s\n', g, group_root);

    subs = dir(fullfile(group_root, 'sub*'));
    subs = subs([subs.isdir]);

    rows = {};
    ALFF_vals = [];
    fALFF_vals = [];

    for s = 1:numel(subs)
        sid = subs(s).name;

        alff_file  = fullfile(group_root, 'Results', 'ALFF_FunImgARCW',  sprintf('zALFFMap_%s.nii', sid));
        falff_file = fullfile(group_root, 'Results', 'fALFF_FunImgARCW', sprintf('zfALFFMap_%s.nii', sid));

        if ~exist(alff_file, 'file')
            gz = [alff_file '.gz'];
            if exist(gz, 'file'), alff_file = gz; else, alff_file = ''; end
        end

        if ~exist(falff_file, 'file')
            gz = [falff_file '.gz'];
            if exist(gz, 'file'), falff_file = gz; else, falff_file = ''; end
        end

        if isempty(alff_file) && isempty(falff_file)
            continue;
        end

        arow = nan(1, numel(alff_roi_idx));
        frow = nan(1, numel(falff_roi_idx));

        if ~isempty(alff_file)
            A = double(read_nifti_auto(alff_file));
            if isequal(size(A), size(atlas))
                for r = 1:numel(alff_roi_idx)
                    roi = alff_roi_idx(r);
                    if isnan(roi), continue; end
                    mask = atlas == roi;
                    arow(r) = mean_valid(A(mask));
                end
            end
        end

        if ~isempty(falff_file)
            F = double(read_nifti_auto(falff_file));
            if isequal(size(F), size(atlas))
                for r = 1:numel(falff_roi_idx)
                    roi = falff_roi_idx(r);
                    if isnan(roi), continue; end
                    mask = atlas == roi;
                    frow(r) = mean_valid(F(mask));
                end
            end
        end

        rows(end+1,:) = {g, sid, sprintf('%d_%s', g, sid)}; %#ok<AGROW>
        ALFF_vals(end+1,:) = arow; %#ok<AGROW>
        fALFF_vals(end+1,:) = frow; %#ok<AGROW>
    end

    if isempty(rows), continue; end

    X = table;
    X.group_num = cell2mat(rows(:,1));
    X.sub = string(rows(:,2));
    X.tag = string(rows(:,3));

    for r = 1:numel(alff_roi_idx)
        roi = alff_roi_idx(r);
        if isnan(roi), continue; end
        vname = safe_var_name(['zALFF_' aal_labels{roi}], sprintf('zALFF_ROI%d', roi));
        X.(vname) = ALFF_vals(:,r);
    end

    for r = 1:numel(falff_roi_idx)
        roi = falff_roi_idx(r);
        if isnan(roi), continue; end
        vname = safe_var_name(['zfALFF_' aal_labels{roi}], sprintf('zfALFF_ROI%d', roi));
        X.(vname) = fALFF_vals(:,r);
    end

    T = outerjoin(T, X, 'Keys', 'tag', 'MergeKeys', true);
end

%% ===================== 只保留 AD（EOAD + LOAD） =====================
Tad = T(T.group_num==1 | T.group_num==2, :);

%% ===================== 影像变量清单 =====================
all_vars = Tad.Properties.VariableNames;
img_vars = {};

% GRETNA vars
for m = 1:size(gretna_specs,1)
    metric_dir_name = gretna_specs{m,1};
    for r = 1:numel(gretna_roi_idx)
        roi = gretna_roi_idx(r);
        if isnan(roi), continue; end
        vname = safe_var_name([metric_dir_name '_' aal_labels{roi}], sprintf('%s_ROI%d', metric_dir_name, roi));
        if ismember(vname, all_vars)
            img_vars{end+1} = vname; %#ok<AGROW>
        end
    end
end

% zALFF vars
for r = 1:numel(alff_roi_idx)
    roi = alff_roi_idx(r);
    if isnan(roi), continue; end
    vname = safe_var_name(['zALFF_' aal_labels{roi}], sprintf('zALFF_ROI%d', roi));
    if ismember(vname, all_vars)
        img_vars{end+1} = vname; %#ok<AGROW>
    end
end

% zfALFF vars
for r = 1:numel(falff_roi_idx)
    roi = falff_roi_idx(r);
    if isnan(roi), continue; end
    vname = safe_var_name(['zfALFF_' aal_labels{roi}], sprintf('zfALFF_ROI%d', roi));
    if ismember(vname, all_vars)
        img_vars{end+1} = vname; %#ok<AGROW>
    end
end

fprintf('\n最终用于临床关联的影像变量数: %d\n', numel(img_vars));

%% ===================== 跑临床关联 =====================
assoc_rows = {};
p_all = [];

for i = 1:numel(img_vars)
    img = img_vars{i};
    x = Tad.(img);

    for j = 1:numel(scale_names)
        sc = scale_names{j};
        y = Tad.(sc);

        ok = isfinite(x) & isfinite(y) & isfinite(Tad.age) & isfinite(Tad.sex) & isfinite(Tad.group_num);

        if sum(ok) < 12
            assoc_rows(end+1,:) = {img, sc, sum(ok), NaN, NaN, NaN, NaN, NaN, NaN}; %#ok<AGROW>
            p_all(end+1,1) = NaN; %#ok<AGROW>
            continue;
        end

        D = table(y(ok), x(ok), categorical(Tad.group_num(ok)), Tad.age(ok), Tad.sex(ok), ...
            'VariableNames', {'scale','img','group','age','sex'});

        lm = fitlm(D, 'scale ~ img + group + age + sex');
        [b1, t1, p1] = get_img_coef(lm);

        [rho, p_rho] = corr(x(ok), y(ok), 'type', 'Spearman', 'rows', 'complete');

        ok2 = ok;
        z = nan(size(x));
        z(ok) = zscore(x(ok));
        ok2(abs(z) > 3) = false;

        if sum(ok2) >= 12
            D2 = table(y(ok2), x(ok2), categorical(Tad.group_num(ok2)), Tad.age(ok2), Tad.sex(ok2), ...
                'VariableNames', {'scale','img','group','age','sex'});
            lm2 = fitlm(D2, 'scale ~ img + group + age + sex');
            [~, ~, p2] = get_img_coef(lm2);
        else
            p2 = NaN;
        end

        assoc_rows(end+1,:) = {img, sc, sum(ok), b1, t1, p1, rho, p_rho, p2}; %#ok<AGROW>
        p_all(end+1,1) = p1; %#ok<AGROW>
    end
end

q_all = fill_q(p_all);

%% ===================== 输出 =====================
out_all = fullfile(out_dir, 'gretna_zalff_zfalff_clinical_AD_only.txt');
fid = fopen(out_all, 'w');
fprintf(fid, ['Imaging\tScale\tN\tBeta_raw\tT_raw\tP_raw\tQ_raw\t' ...
    'SpearmanRho\tP_spearman\tP_cleanOutlier\n']);

for k = 1:size(assoc_rows,1)
    fprintf(fid, '%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', ...
        assoc_rows{k,1}, assoc_rows{k,2}, assoc_rows{k,3}, ...
        ns(assoc_rows{k,4}), ns(assoc_rows{k,5}), ns(assoc_rows{k,6}), ns(q_all(k)), ...
        ns(assoc_rows{k,7}), ns(assoc_rows{k,8}), ns(assoc_rows{k,9}));
end
fclose(fid);

out_sig = fullfile(out_dir, 'gretna_zalff_zfalff_clinical_AD_only_FDRsig.txt');
fid = fopen(out_sig, 'w');
fprintf(fid, 'Imaging\tScale\tN\tBeta_raw\tT_raw\tP_raw\tQ_raw\tSpearmanRho\tP_spearman\tP_cleanOutlier\n');

for k = 1:size(assoc_rows,1)
    if isfinite(q_all(k)) && q_all(k) < 0.05
        fprintf(fid, '%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', ...
            assoc_rows{k,1}, assoc_rows{k,2}, assoc_rows{k,3}, ...
            ns(assoc_rows{k,4}), ns(assoc_rows{k,5}), ns(assoc_rows{k,6}), ns(q_all(k)), ...
            ns(assoc_rows{k,7}), ns(assoc_rows{k,8}), ns(assoc_rows{k,9}));
    end
end
fclose(fid);

fprintf('\n完成。\n');
fprintf('总表: %s\n', out_all);
fprintf('FDR显著: %s\n', out_sig);

end

%% ===================== 子函数 =====================

function subs = sort_sub_ids(subs)
subs = cellstr(string(subs));
nums = nan(numel(subs),1);
for i = 1:numel(subs)
    tok = regexp(subs{i}, 'sub(\d+)', 'tokens', 'once');
    if ~isempty(tok)
        nums(i) = str2double(tok{1});
    end
end
[~, idx] = sort(nums);
subs = string(subs(idx));
end

function nii = read_nifti_auto(f)
if endsWith(lower(f), '.gz')
    tmpdir = tempname;
    mkdir(tmpdir);
    gunzip(f, tmpdir);
    [~,name,~] = fileparts(f);
    if endsWith(lower(name), '.nii')
        nii_file = fullfile(tmpdir, name);
    else
        L = dir(fullfile(tmpdir, '*.nii'));
        nii_file = fullfile(L(1).folder, L(1).name);
    end
    nii = niftiread(nii_file);
else
    nii = niftiread(f);
end
end

function idx = labels_to_idx(target_labels, all_labels)
idx = nan(numel(target_labels),1);
for i = 1:numel(target_labels)
    hit = find(strcmpi(strtrim(target_labels{i}), strtrim(all_labels)), 1);
    if ~isempty(hit)
        idx(i) = hit;
    end
end
end

function labels = read_labels(label_file, nROI)
labels = arrayfun(@(x) sprintf('ROI_%d',x), 1:nROI, 'UniformOutput', false)';

if ~exist(label_file,'file')
    return;
end

fid = fopen(label_file,'r');
if fid < 0
    return;
end

k = 0;
while ~feof(fid) && k < nROI
    line = strtrim(fgetl(fid));
    if isempty(line)
        continue;
    end

    % 格式示例：
    % 97 Cerebelum_4_5_L 9031
    parts = strsplit(line);
    if numel(parts) >= 2
        k = k + 1;
        labels{k} = strtrim(parts{2});
    end
end

fclose(fid);
end
function [b, t, p] = get_img_coef(lm)
b = NaN; t = NaN; p = NaN;
coef = lm.Coefficients;
rn = coef.Properties.RowNames;
idx = find(strcmp(rn, 'img') | contains(rn,'img'), 1);
if ~isempty(idx)
    b = coef.Estimate(idx);
    t = coef.tStat(idx);
    p = coef.pValue(idx);
end
end

function q = fill_q(p)
q = nan(size(p));
ok = isfinite(p);
if any(ok)
    q(ok) = mafdr(p(ok), 'BHFDR', true);
end
end

function x = to_num(v)
if isnumeric(v)
    x = v;
elseif ischar(v) || isstring(v)
    x = str2double(string(v));
else
    x = NaN;
end
if isempty(x), x = NaN; end
end

function s = ns(x)
if isnan(x)
    s = 'NaN';
else
    s = sprintf('%.6g', x);
end
end

function n = safe_var_name(x, fallback)
try
    if ismissing(x)
        n = fallback;
        return;
    end
catch
end

if ischar(x)
    s = x;
elseif isstring(x)
    s = char(x);
elseif isnumeric(x)
    s = num2str(x);
else
    s = fallback;
end

s = strtrim(s);
if isempty(s)
    s = fallback;
end

n = matlab.lang.makeValidName(s);
if isempty(n)
    n = fallback;
end
end

function m = mean_valid(x)
x = x(isfinite(x));
if isempty(x)
    m = NaN;
else
    m = mean(x);
end
end
