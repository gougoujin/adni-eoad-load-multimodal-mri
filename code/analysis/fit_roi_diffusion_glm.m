% Repository usage summary
% Description: Fit age/sex-adjusted ROI-wise diffusion metric models.
% Usage: Run after extract_roi_diffusion_metrics has generated cohort matrices.
% Outputs: Outputs ROI-wise diffusion statistics and FDR-adjusted summaries.
% Note: Update input paths, toolboxes, and filenames for your local environment.

% Public repository version: update file paths, toolboxes, and local settings before running.
% This script/function was lightly sanitized for sharing and may require project-specific inputs.

function fit_roi_diffusion_glm()
clear; clc;

%% ===================== 路径 =====================
data_dir   = '';
cov_file   = '';          % 按实际修改
label_file = '';  % 按实际修改

fa_mat = fullfile(data_dir, 'fa_roi_all_subjects.mat');
md_mat = fullfile(data_dir, 'md_roi_all_subjects.mat');

out_fa_txt     = fullfile(data_dir, 'roi_level_fa_stats_age_sex_fixed.txt');
out_fa_fdr_txt = fullfile(data_dir, 'roi_level_fa_stats_age_sex_fixed_fdr.txt');

out_md_txt     = fullfile(data_dir, 'roi_level_md_stats_age_sex_fixed.txt');
out_md_fdr_txt = fullfile(data_dir, 'roi_level_md_stats_age_sex_fixed_fdr.txt');

%% ===================== 标签 =====================
labels = read_labels(label_file, 116);

%% ===================== 跑 FA =====================
S = load(fa_mat);
run_one_metric_fixed(S.subject_info, S.fa_all, cov_file, labels, 'FA', out_fa_txt, out_fa_fdr_txt);

%% ===================== 跑 MD =====================
S = load(md_mat);
run_one_metric_fixed(S.subject_info, S.md_all, cov_file, labels, 'MD', out_md_txt, out_md_fdr_txt);

fprintf('\n全部完成。\n');
end


function run_one_metric_fixed(subject_info, X, cov_file, labels, metric_name, out_txt, out_fdr_txt)

%% ========= 从提取结果中读 group/sub，作为唯一主索引 =========
nSub = size(subject_info, 1);
nROI = size(X, 2);

group_name = strings(nSub,1);
sub_id     = strings(nSub,1);
group_num  = nan(nSub,1);
tag_data   = strings(nSub,1);

for i = 1:nSub
    g = string(subject_info{i,1});   % 例如 group1
    s = string(subject_info{i,2});   % 例如 sub001

    group_name(i) = g;
    sub_id(i)     = s;

    tok = regexp(char(g), 'group(\d+)', 'tokens', 'once');
    if ~isempty(tok)
        group_num(i) = str2double(tok{1});
    end

    tag_data(i) = string(group_num(i)) + "_" + s;
end

%% ========= 读取协变量表 =========
% 默认：
% 第1列：group
% 第3列：sub id
% 第4列：age
% 第5列：sex
raw = readcell(cov_file);
raw = raw(2:end, :);

n = size(raw,1);

cov_group = nan(n,1);
cov_sub   = strings(n,1);
cov_age   = nan(n,1);
cov_sex_num = nan(n,1);
cov_sex_str = strings(n,1);

for i = 1:n
    % ---- group ----
    g = raw{i,1};
    if isnumeric(g)
        cov_group(i) = g;
    else
        tok = regexp(char(string(g)), '(\d+)', 'tokens', 'once');
        if ~isempty(tok)
            cov_group(i) = str2double(tok{1});
        end
    end

    % ---- sub id ----
    s = raw{i,3};
    if isnumeric(s)
        if ~isnan(s)
            cov_sub(i) = sprintf('sub%03d', round(s));
        end
    else
        tok = regexp(char(string(s)), '(\d+)', 'tokens', 'once');
        if ~isempty(tok)
            cov_sub(i) = sprintf('sub%03d', str2double(tok{1}));
        end
    end

    % ---- age ----
    a = raw{i,4};
    if isnumeric(a)
        cov_age(i) = a;
    else
        tmp = str2double(string(a));
        if ~isnan(tmp)
            cov_age(i) = tmp;
        end
    end

    % ---- sex ----
    sx = raw{i,5};
    if isnumeric(sx)
        cov_sex_num(i) = sx;
    else
        sxs = string(sx);
        cov_sex_str(i) = sxs;
        tmp = str2double(sxs);
        if ~isnan(tmp)
            cov_sex_num(i) = tmp;
        end
    end
end

C = table;
C.group_num = cov_group;
C.sub       = cov_sub;
C.age       = cov_age;

% sex 尽量转数值，否则保留分类
if sum(~isnan(cov_sex_num)) >= max(5, round(0.8*n))
    C.sex = cov_sex_num;
    sex_is_numeric = true;
else
    C.sex = categorical(cov_sex_str);
    sex_is_numeric = false;
end

C.tag = string(C.group_num) + "_" + C.sub;

keepC = ~isnan(C.group_num) & strlength(C.sub) > 0 & ~isnan(C.age);
if sex_is_numeric
    keepC = keepC & ~isnan(C.sex);
else
    keepC = keepC & ~ismissing(C.sex);
end
C = C(keepC,:);

%% ========= 和提取结果做唯一键匹配 =========
[common_tag, ia, ib] = intersect(tag_data, C.tag, 'stable');

X = X(ia, :);
group_name = group_name(ia);
sub_id     = sub_id(ia);
group_num  = group_num(ia);
Cm         = C(ib,:);

fprintf('\n[%s] Matched subjects = %d\n', metric_name, numel(common_tag));
fprintf('[%s] Group1 = %d, Group2 = %d, Group3 = %d\n', metric_name, ...
    sum(group_num==1), sum(group_num==2), sum(group_num==3));

%% ========= 额外检查：同一个 tag 是否重复 =========
if numel(unique(common_tag)) ~= numel(common_tag)
    warning('[%s] 匹配后仍有重复 tag，请检查协变量表。', metric_name);
end

%% ========= ROI-wise regression =========
roi_idx = (1:nROI)';
mean_g1 = nan(nROI,1);
mean_g2 = nan(nROI,1);
mean_g3 = nan(nROI,1);

n_g1 = nan(nROI,1);
n_g2 = nan(nROI,1);
n_g3 = nan(nROI,1);

p_anova = nan(nROI,1);

beta_g2 = nan(nROI,1);
t_g2    = nan(nROI,1);
p_g2    = nan(nROI,1);

beta_g3 = nan(nROI,1);
t_g3    = nan(nROI,1);
p_g3    = nan(nROI,1);

for r = 1:nROI
    y = X(:,r);

    if sex_is_numeric
        D = table(y, categorical(group_num), Cm.age, Cm.sex, ...
            'VariableNames', {'y','group','age','sex'});
        keep = isfinite(D.y) & isfinite(D.age) & isfinite(D.sex);
    else
        D = table(y, categorical(group_num), Cm.age, Cm.sex, ...
            'VariableNames', {'y','group','age','sex'});
        keep = isfinite(D.y) & isfinite(D.age) & ~ismissing(D.sex);
    end

    D = D(keep,:);

    if height(D) < 10
        continue;
    end

    % 每个 ROI 的真实样本数
    mean_g1(r) = mean(D.y(D.group=='1'), 'omitnan');
    mean_g2(r) = mean(D.y(D.group=='2'), 'omitnan');
    mean_g3(r) = mean(D.y(D.group=='3'), 'omitnan');

    n_g1(r) = sum(D.group=='1');
    n_g2(r) = sum(D.group=='2');
    n_g3(r) = sum(D.group=='3');

    % 至少两组有数据才有意义
    if sum([n_g1(r)>0, n_g2(r)>0, n_g3(r)>0]) < 2
        continue;
    end

    try
        lm = fitlm(D, 'y ~ group + age + sex');
        A = anova(lm, 'summary');

        ridx = find(strcmpi(A.Properties.RowNames, 'group'), 1);
        if isempty(ridx)
            ridx = find(contains(lower(A.Properties.RowNames), 'group'), 1);
        end
        if ~isempty(ridx)
            p_anova(r) = A.pValue(ridx);
        end

        coef = lm.Coefficients;
        rn = coef.Properties.RowNames;

        % 兼容不同命名风格
        idx2 = find(strcmp(rn, 'group_2') | strcmp(rn, 'group_2.0') | strcmp(rn, 'group_2 '), 1);
        idx3 = find(strcmp(rn, 'group_3') | strcmp(rn, 'group_3.0') | strcmp(rn, 'group_3 '), 1);

        if isempty(idx2)
            idx2 = find(contains(rn, 'group_2'), 1);
        end
        if isempty(idx3)
            idx3 = find(contains(rn, 'group_3'), 1);
        end

        if ~isempty(idx2)
            beta_g2(r) = coef.Estimate(idx2);
            t_g2(r)    = coef.tStat(idx2);
            p_g2(r)    = coef.pValue(idx2);
        end

        if ~isempty(idx3)
            beta_g3(r) = coef.Estimate(idx3);
            t_g3(r)    = coef.tStat(idx3);
            p_g3(r)    = coef.pValue(idx3);
        end

    catch ME
        fprintf('[WARN] %s ROI %d 拟合失败: %s\n', metric_name, r, ME.message);
    end
end

%% ========= FDR =========
q_anova = nan(nROI,1);
q_g2    = nan(nROI,1);
q_g3    = nan(nROI,1);

ok1 = isfinite(p_anova);
ok2 = isfinite(p_g2);
ok3 = isfinite(p_g3);

if any(ok1), q_anova(ok1) = mafdr(p_anova(ok1), 'BHFDR', true); end
if any(ok2), q_g2(ok2)    = mafdr(p_g2(ok2),    'BHFDR', true); end
if any(ok3), q_g3(ok3)    = mafdr(p_g3(ok3),    'BHFDR', true); end

%% ========= 写总表 =========
fid = fopen(out_txt, 'w');
if fid < 0
    error('无法写出文件: %s', out_txt);
end

fprintf(fid, ['ROI\tLabel\tMean_group1\tMean_group2\tMean_group3\t' ...
    'N_group1\tN_group2\tN_group3\t' ...
    'P_ANOVA\tBeta_g2_vs_g1\tT_g2_vs_g1\tP_g2_vs_g1\t' ...
    'Beta_g3_vs_g1\tT_g3_vs_g1\tP_g3_vs_g1\t' ...
    'Q_ANOVA\tQ_g2_vs_g1\tQ_g3_vs_g1\n']);

for r = 1:nROI
    fprintf(fid, '%d\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', ...
        roi_idx(r), labels{r}, ...
        num2str_nan(mean_g1(r)), num2str_nan(mean_g2(r)), num2str_nan(mean_g3(r)), ...
        round_nan(n_g1(r)), round_nan(n_g2(r)), round_nan(n_g3(r)), ...
        num2str_nan(p_anova(r)), ...
        num2str_nan(beta_g2(r)), num2str_nan(t_g2(r)), num2str_nan(p_g2(r)), ...
        num2str_nan(beta_g3(r)), num2str_nan(t_g3(r)), num2str_nan(p_g3(r)), ...
        num2str_nan(q_anova(r)), num2str_nan(q_g2(r)), num2str_nan(q_g3(r)));
end
fclose(fid);

%% ========= 写 FDR 通过表 =========
fid = fopen(out_fdr_txt, 'w');
if fid < 0
    error('无法写出文件: %s', out_fdr_txt);
end

fprintf(fid, 'ROI\tLabel\tContrast\tBeta\tT\tP\tQ\n');

for r = 1:nROI
    if isfinite(q_g2(r)) && q_g2(r) < 0.05
        fprintf(fid, '%d\t%s\tgroup2_vs_group1\t%.6f\t%.6f\t%.6g\t%.6g\n', ...
            roi_idx(r), labels{r}, beta_g2(r), t_g2(r), p_g2(r), q_g2(r));
    end
    if isfinite(q_g3(r)) && q_g3(r) < 0.05
        fprintf(fid, '%d\t%s\tgroup3_vs_group1\t%.6f\t%.6f\t%.6g\t%.6g\n', ...
            roi_idx(r), labels{r}, beta_g3(r), t_g3(r), p_g3(r), q_g3(r));
    end
end

fclose(fid);

fprintf('[%s] 完成：\n%s\n%s\n', metric_name, out_txt, out_fdr_txt);
end


function labels = read_labels(label_file, nROI)
labels = cell(nROI,1);
for i = 1:nROI
    labels{i} = sprintf('ROI_%d', i);
end

if ~exist(label_file, 'file')
    warning('找不到标签文件：%s', label_file);
    return;
end

fid = fopen(label_file, 'r');
if fid < 0
    warning('无法打开标签文件：%s', label_file);
    return;
end

k = 0;
while ~feof(fid) && k < nROI
    line = strtrim(fgetl(fid));
    if isempty(line), continue; end
    k = k + 1;
    labels{k} = line;
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


function x = round_nan(v)
if isnan(v)
    x = 0;
else
    x = round(v);
end
end
