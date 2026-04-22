% Repository usage summary
% Description: Run sensitivity and clinical-association models focused on coupling results.
% Usage: Run after global and regional coupling outputs are available.
% Outputs: Writes sensitivity tables and AD-only clinical association summaries.
% Note: Update input paths, toolboxes, and filenames for your local environment.

% Public repository version: update file paths, toolboxes, and local settings before running.
% This script/function was lightly sanitized for sharing and may require project-specific inputs.

function run_coupling_sensitivity_and_clinical_models()
clear; clc;

%% ===================== 路径 =====================
cov_file      = '';
mean_mat_file = '';
roi_mat_file  = '';

out_dir = '';
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

%% ===================== 读取协变量与临床量表 =====================
% 第1列=组别, 第3列=sub, 第4列=年龄, 第5列=性别, 第6:11列=量表
raw = readcell(cov_file);
header = raw(1,:);
raw = raw(2:end,:);

n = size(raw,1);

cov_group = nan(n,1);
cov_sub   = cell(n,1);
cov_age   = nan(n,1);
cov_sex   = nan(n,1);

for i = 1:n
    cov_group(i) = to_num(raw{i,1});

    s = raw{i,3};
    if isnumeric(s) && ~isnan(s)
        cov_sub{i} = sprintf('sub%03d', round(s));
    else
        tok = regexp(char(string(s)), '(\d+)', 'tokens', 'once');
        if ~isempty(tok)
            cov_sub{i} = sprintf('sub%03d', str2double(tok{1}));
        else
            cov_sub{i} = '';
        end
    end

    cov_age(i) = to_num(raw{i,4});
    cov_sex(i) = to_num(raw{i,5});
end

scale_cols = 6:11;   % F:K
nScale = numel(scale_cols);

scale_names = cell(nScale,1);
for j = 1:nScale
    scale_names{j} = safe_var_name(header{scale_cols(j)}, sprintf('Scale%d', j));
end

scale_mat = nan(n, nScale);
for j = 1:nScale
    for i = 1:n
        scale_mat(i,j) = to_num(raw{i, scale_cols(j)});
    end
end

C = table;
C.group_num = cov_group;
C.sub       = string(cov_sub);
C.age       = cov_age;
C.sex       = cov_sex;

for j = 1:nScale
    C.(scale_names{j}) = scale_mat(:,j);
end

C.tag = string(C.group_num) + "_" + C.sub;
C = C(isfinite(C.group_num) & strlength(C.sub)>0 & isfinite(C.age) & isfinite(C.sex), :);

%% ===================== 读取 mean coupling =====================
S = load(mean_mat_file);

matched_mean = [];
subj_mean = [];

if isfield(S, 'results')
    Rm = S.results;
    matched_mean = Rm(:,1:2);
    subj_mean = cellfun(@to_num, Rm(:,3));
elseif isfield(S, 'matched') && isfield(S, 'all_roi_coupling')
    matched_mean = S.matched(:,1:2);
    subj_mean = mean(S.all_roi_coupling, 2, 'omitnan');
else
    error('sc_fc_coupling_all.mat 中未找到可用变量');
end

mean_group_num = nan(size(matched_mean,1),1);
mean_sub = strings(size(matched_mean,1),1);

for i = 1:size(matched_mean,1)
    g = string(matched_mean{i,1});
    s = string(matched_mean{i,2});
    tok = regexp(char(g), 'group(\d+)', 'tokens', 'once');
    mean_group_num(i) = str2double(tok{1});
    mean_sub(i) = s;
end

MeanT = table;
MeanT.group_num = mean_group_num;
MeanT.sub = mean_sub;
MeanT.mean_coupling = subj_mean;
MeanT.tag = string(MeanT.group_num) + "_" + MeanT.sub;

%% ===================== 读取 ROI-level coupling =====================
R = load(roi_mat_file);

if ~isfield(R, 'matched') || ~isfield(R, 'roi_coupling_all')
    error('roi_level_sc_fc_coupling_all.mat 中缺少 matched 或 roi_coupling_all');
end

n_match = size(R.matched, 1);
A = R.roi_coupling_all;
sa = size(A);

roi_group_num = nan(n_match,1);
roi_sub = strings(n_match,1);

for i = 1:n_match
    g = string(R.matched{i,1});
    s = string(R.matched{i,2});
    tok = regexp(char(g), 'group(\d+)', 'tokens', 'once');
    roi_group_num(i) = str2double(tok{1});
    roi_sub(i) = s;
end

% 自动判断方向
if numel(sa) ~= 2
    error('roi_coupling_all 不是二维矩阵');
end

if sa(1) == n_match && sa(2) >= 23
    B = A;      % Nsub × ROI
elseif sa(2) == n_match && sa(1) >= 23
    B = A';     % ROI × Nsub -> Nsub × ROI
else
    error('roi_coupling_all 尺寸无法和 matched 对齐：matched=%d, roi_coupling_all=[%d %d]', ...
        n_match, sa(1), sa(2));
end

RoiT = table;
RoiT.group_num = roi_group_num;
RoiT.sub = roi_sub;
RoiT.tag = string(RoiT.group_num) + "_" + RoiT.sub;
RoiT.ROI8  = B(:,8);
RoiT.ROI10 = B(:,10);
RoiT.ROI14 = B(:,14);
RoiT.ROI23 = B(:,23);

%% ===================== 合并总表 =====================
T = outerjoin(C, MeanT, 'Keys','tag', 'MergeKeys', true);
T = outerjoin(T, RoiT(:, {'tag','ROI8','ROI10','ROI14','ROI23'}), 'Keys','tag', 'MergeKeys', true);

keep = isfinite(T.mean_coupling) | isfinite(T.ROI8) | isfinite(T.ROI10) | isfinite(T.ROI14) | isfinite(T.ROI23);
T = T(keep,:);

vnames = T.Properties.VariableNames;

if ismember('group_num_MeanT', vnames)
    T.group_num = T.group_num_MeanT;
elseif ismember('group_num_C', vnames)
    T.group_num = T.group_num_C;
elseif ~ismember('group_num', vnames)
    error('合并后未找到 group_num');
end

if ismember('sub_MeanT', vnames)
    T.sub = T.sub_MeanT;
elseif ismember('sub_C', vnames)
    T.sub = T.sub_C;
elseif ~ismember('sub', vnames)
    error('合并后未找到 sub');
end

keep_cols = {'tag','group_num','sub','age','sex','mean_coupling','ROI8','ROI10','ROI14','ROI23'};
for j = 1:nScale
    if ismember(scale_names{j}, T.Properties.VariableNames)
        keep_cols{end+1} = scale_names{j}; %#ok<AGROW>
    end
end
T = T(:, keep_cols);

%% ===================== 1) 敏感性分析：EOAD vs LOAD =====================
T12 = T(T.group_num==1 | T.group_num==2, :);
img_vars = {'mean_coupling','ROI8','ROI10','ROI14','ROI23'};

sens_file = fullfile(out_dir, 'coupling_sensitivity_EOAD_vs_LOAD.txt');
fid = fopen(sens_file, 'w');
fprintf(fid, ['Metric\tN_raw\tBeta_LOAD_vs_EOAD_raw\tT_raw\tP_raw\t' ...
    'N_clean\tBeta_LOAD_vs_EOAD_clean\tT_clean\tP_clean\n']);

for i = 1:numel(img_vars)
    v = img_vars{i};
    y = T12.(v);

    ok = isfinite(y) & isfinite(T12.age) & isfinite(T12.sex) & isfinite(T12.group_num);
    D = table(y(ok), categorical(T12.group_num(ok)), T12.age(ok), T12.sex(ok), ...
        'VariableNames', {'y','group','age','sex'});

    b_raw = NaN; t_raw = NaN; p_raw = NaN; n_raw = height(D);
    if height(D) >= 10
        lm = fitlm(D, 'y ~ group + age + sex');
        [b_raw, t_raw, p_raw] = get_group2_coef(lm);
    end

    ok2 = ok;
    z = nan(size(y));
    z(ok) = zscore(y(ok));
    ok2(abs(z) > 3) = false;

    D2 = table(y(ok2), categorical(T12.group_num(ok2)), T12.age(ok2), T12.sex(ok2), ...
        'VariableNames', {'y','group','age','sex'});

    b_clean = NaN; t_clean = NaN; p_clean = NaN; n_clean = height(D2);
    if height(D2) >= 10
        lm2 = fitlm(D2, 'y ~ group + age + sex');
        [b_clean, t_clean, p_clean] = get_group2_coef(lm2);
    end

    fprintf(fid, '%s\t%d\t%s\t%s\t%s\t%d\t%s\t%s\t%s\n', ...
        v, n_raw, ns(b_raw), ns(t_raw), ns(p_raw), ...
        n_clean, ns(b_clean), ns(t_clean), ns(p_clean));
end
fclose(fid);

%% ===================== 2) 临床量表关联：仅 AD（EOAD+LOAD） =====================
Tad = T(T.group_num==1 | T.group_num==2, :);

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

assoc_file = fullfile(out_dir, 'clinical_association_AD_only.txt');
fid = fopen(assoc_file, 'w');
fprintf(fid, ['Imaging\tScale\tN\tBeta_raw\tT_raw\tP_raw\tQ_raw\t' ...
    'SpearmanRho\tP_spearman\tP_cleanOutlier\n']);

for k = 1:size(assoc_rows,1)
    fprintf(fid, '%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', ...
        assoc_rows{k,1}, assoc_rows{k,2}, assoc_rows{k,3}, ...
        ns(assoc_rows{k,4}), ns(assoc_rows{k,5}), ns(assoc_rows{k,6}), ns(q_all(k)), ...
        ns(assoc_rows{k,7}), ns(assoc_rows{k,8}), ns(assoc_rows{k,9}));
end
fclose(fid);

sig_file = fullfile(out_dir, 'clinical_association_AD_only_FDRsig.txt');
fid = fopen(sig_file, 'w');
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
fprintf('敏感性分析: %s\n', sens_file);
fprintf('临床关联总表: %s\n', assoc_file);
fprintf('临床关联FDR显著: %s\n', sig_file);

end

%% ===================== 子函数 =====================

function [b, t, p] = get_group2_coef(lm)
b = NaN; t = NaN; p = NaN;
coef = lm.Coefficients;
rn = coef.Properties.RowNames;

idx = find(strcmp(rn, 'group_2') | contains(rn,'group_2'), 1);
if ~isempty(idx)
    b = coef.Estimate(idx);
    t = coef.tStat(idx);
    p = coef.pValue(idx);
end
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
if isempty(x)
    x = NaN;
end
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
