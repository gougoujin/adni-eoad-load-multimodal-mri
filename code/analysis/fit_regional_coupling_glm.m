% Repository usage summary
% Description: Fit ROI-wise age/sex-adjusted group models for regional structural-functional coupling.
% Usage: Run after regional coupling has been computed.
% Outputs: Writes ROI-wise coefficient tables and FDR-adjusted summary files.
% Note: Update input paths, toolboxes, and filenames for your local environment.

% Public repository version: update file paths, toolboxes, and local settings before running.
% This script/function was lightly sanitized for sharing and may require project-specific inputs.

clear; clc;

%% ========= paths =========
roi_mat_file = '';
cov_file     = '';

out_summary_txt = '';
out_fdr_txt     = '';
out_mat         = '';

nROI = 116;

%% ========= load ROI coupling =========
load(roi_mat_file, 'matched', 'roi_coupling_all');

nSub = size(matched, 1);

group_name = strings(nSub,1);
sub_id     = strings(nSub,1);
group_num  = nan(nSub,1);

for i = 1:nSub
    group_name(i) = string(matched{i,1});
    sub_id(i)     = string(matched{i,2});

    tok = regexp(char(group_name(i)), 'group(\d+)', 'tokens', 'once');
    if ~isempty(tok)
        group_num(i) = str2double(tok{1});
    end
end

tag_scfc = string(group_num) + "_" + sub_id;

%% ========= read covariates by fixed columns =========
% col 1 = group
% col 3 = sub id
% col 4 = age
% col 5 = sex
raw = readcell(cov_file);
raw = raw(2:end, :);

n = size(raw,1);

cov_group = nan(n,1);
cov_sub   = strings(n,1);
cov_age   = nan(n,1);
cov_sex_num = nan(n,1);
cov_sex_str = strings(n,1);

for i = 1:n
    % group
    g = raw{i,1};
    if isnumeric(g)
        cov_group(i) = g;
    else
        tok = regexp(char(string(g)), '(\d+)', 'tokens', 'once');
        if ~isempty(tok)
            cov_group(i) = str2double(tok{1});
        end
    end

    % sub id
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

    % age
    a = raw{i,4};
    if isnumeric(a)
        cov_age(i) = a;
    else
        tmp = str2double(string(a));
        if ~isnan(tmp)
            cov_age(i) = tmp;
        end
    end

    % sex
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
C.sub = cov_sub;
C.age = cov_age;

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

%% ========= match subjects =========
[common_tag, ia, ib] = intersect(tag_scfc, C.tag, 'stable');

roi_coupling = roi_coupling_all(:, ia);
group_name = group_name(ia);
sub_id = sub_id(ia);
group_num = group_num(ia);
Cm = C(ib,:);

fprintf('Matched subjects for ROI stats: %d\n', numel(common_tag));

%% ========= ROI-wise regression =========
roi_idx = (1:nROI)';
n_valid = nan(nROI,1);

beta_g2 = nan(nROI,1);
p_g2    = nan(nROI,1);
beta_g3 = nan(nROI,1);
p_g3    = nan(nROI,1);
beta_age = nan(nROI,1);
p_age    = nan(nROI,1);
beta_sex = nan(nROI,1);
p_sex    = nan(nROI,1);

for r = 1:nROI
    y = roi_coupling(r, :)';

    D = table;
    D.y = y;
    D.group = categorical(group_num);   % reference = group1
    D.age = Cm.age;
    D.sex = Cm.sex;

    keep = ~isnan(D.y) & ~isnan(D.age);
    if isnumeric(D.sex)
        keep = keep & ~isnan(D.sex);
    else
        keep = keep & ~ismissing(D.sex);
    end

    D = D(keep,:);
    n_valid(r) = height(D);

    if height(D) < 10
        continue;
    end

    mdl = fitlm(D, 'y ~ group + age + sex');
    coef = mdl.Coefficients;
    rn = string(coef.Properties.RowNames);

    idx = find(rn == "group_2", 1);
    if ~isempty(idx)
        beta_g2(r) = coef{idx, "Estimate"};
        p_g2(r)    = coef{idx, "pValue"};
    end

    idx = find(rn == "group_3", 1);
    if ~isempty(idx)
        beta_g3(r) = coef{idx, "Estimate"};
        p_g3(r)    = coef{idx, "pValue"};
    end

    idx = find(rn == "age", 1);
    if ~isempty(idx)
        beta_age(r) = coef{idx, "Estimate"};
        p_age(r)    = coef{idx, "pValue"};
    end

    idx = find(rn == "sex", 1);
    if ~isempty(idx)
        beta_sex(r) = coef{idx, "Estimate"};
        p_sex(r)    = coef{idx, "pValue"};
    else
        idx = find(contains(rn, "sex"), 1);
        if ~isempty(idx)
            beta_sex(r) = coef{idx, "Estimate"};
            p_sex(r)    = coef{idx, "pValue"};
        end
    end
end

%% ========= FDR on group2 and group3 contrasts separately =========
[h_g2, q_g2] = local_bh_fdr(p_g2);
[h_g3, q_g3] = local_bh_fdr(p_g3);

%% ========= output tables =========
Tout = table;
Tout.ROI = roi_idx;
Tout.n_valid = n_valid;

Tout.beta_group2 = beta_g2;
Tout.p_group2 = p_g2;
Tout.fdr_pass_group2 = h_g2;
Tout.q_group2 = q_g2;

Tout.beta_group3 = beta_g3;
Tout.p_group3 = p_g3;
Tout.fdr_pass_group3 = h_g3;
Tout.q_group3 = q_g3;

Tout.beta_age = beta_age;
Tout.p_age = p_age;
Tout.beta_sex = beta_sex;
Tout.p_sex = p_sex;

writetable(Tout, out_summary_txt, 'Delimiter', '\t', 'FileType', 'text');
writetable(Tout, out_fdr_txt, 'Delimiter', '\t', 'FileType', 'text');

save(out_mat, 'Tout', 'roi_coupling', 'group_name', 'sub_id', 'group_num', 'Cm');

fprintf('\nOutput file/path/to/project', out_summary_txt, out_fdr_txt, out_mat);

%% ========= local function =========
function [h, q] = local_bh_fdr(p)
    p = p(:);
    h = false(size(p));
    q = nan(size(p));

    valid = ~isnan(p);
    pv = p(valid);

    if isempty(pv)
        return;
    end

    [ps, idx] = sort(pv);
    m = numel(ps);
    thr = (1:m)'/m * 0.05;

    below = ps <= thr;
    if any(below)
        k = find(below, 1, 'last');
        cutoff = ps(k);
        h(valid) = pv <= cutoff;
    end

    qv = ps .* m ./ (1:m)';
    for i = m-1:-1:1
        qv(i) = min(qv(i), qv(i+1));
    end
    qv(qv>1) = 1;

    tmp = nan(size(pv));
    tmp(idx) = qv;
    q(valid) = tmp;
end
