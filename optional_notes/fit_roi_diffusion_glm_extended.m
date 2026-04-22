% Repository usage summary
% Description: Extended ROI-wise diffusion GLM covering additional diffusion metrics.
% Usage: Optional extended diffusion statistics.
% Outputs: Use only if the expanded DTI metric set is needed.
% Note: Update input paths, toolboxes, and filenames for your local environment.

% Public repository version: update file paths, toolboxes, and local settings before running.
% This script/function was lightly sanitized for sharing and may require project-specific inputs.

function fit_roi_diffusion_glm_extended()
clear; clc;

data_file  = '';
cov_file   = '';          % 按实际修改
label_file = '';  % 按实际修改
out_dir    = '';

labels = read_labels(label_file, 116);
S = load(data_file);

run_one(S.subject_info, S.FA_all, 'FA', cov_file, labels, out_dir);
run_one(S.subject_info, S.MD_all, 'MD', cov_file, labels, out_dir);
run_one(S.subject_info, S.AD_all, 'AD', cov_file, labels, out_dir);
run_one(S.subject_info, S.RD_all, 'RD', cov_file, labels, out_dir);
run_one(S.subject_info, S.MO_all, 'MO', cov_file, labels, out_dir);

fprintf('\n全部完成。\n');
end

function run_one(subject_info, X, metric_name, cov_file, labels, out_dir)
raw = readcell(cov_file);
raw = raw(2:end,:);

n = size(raw,1);
cov_group = nan(n,1);
cov_sub   = strings(n,1);
cov_age   = nan(n,1);
cov_sex   = nan(n,1);

for i = 1:n
    cov_group(i) = to_num(raw{i,1});

    s = raw{i,3};
    if isnumeric(s) && ~isnan(s)
        cov_sub(i) = sprintf('sub%03d', round(s));
    else
        tok = regexp(char(string(s)), '(\d+)', 'tokens', 'once');
        if ~isempty(tok), cov_sub(i) = sprintf('sub%03d', str2double(tok{1})); end
    end

    cov_age(i) = to_num(raw{i,4});
    cov_sex(i) = to_num(raw{i,5});
end

C = table(cov_group, cov_sub, cov_age, cov_sex);
C.tag = string(C.cov_group) + "_" + C.cov_sub;
C = C(isfinite(C.cov_group) & strlength(C.cov_sub)>0 & isfinite(C.cov_age) & isfinite(C.cov_sex), :);

nSub = size(subject_info,1);
tag_data = strings(nSub,1);
group_num = nan(nSub,1);

for i = 1:nSub
    g = string(subject_info{i,1});
    s = string(subject_info{i,2});
    tok = regexp(char(g), 'group(\d+)', 'tokens', 'once');
    group_num(i) = str2double(tok{1});
    tag_data(i) = string(group_num(i)) + "_" + s;
end

[~, ia, ib] = intersect(tag_data, C.tag, 'stable');
X = X(ia,:);
group_num = group_num(ia);
C = C(ib,:);

out_txt = fullfile(out_dir, ['roi_level_' lower(metric_name) '_stats_age_sex.txt']);
out_fdr = fullfile(out_dir, ['roi_level_' lower(metric_name) '_stats_age_sex_fdr.txt']);

fid = fopen(out_txt,'w');
fprintf(fid,['ROI\tLabel\tMean_group1\tMean_group2\tMean_group3\tN_group1\tN_group2\tN_group3\t' ...
    'P_ANOVA\tBeta_g2_vs_g1\tT_g2_vs_g1\tP_g2_vs_g1\tBeta_g3_vs_g1\tT_g3_vs_g1\tP_g3_vs_g1\tQ_ANOVA\tQ_g2_vs_g1\tQ_g3_vs_g1\n']);

nROI = size(X,2);
P0 = nan(nROI,1); P2 = nan(nROI,1); P3 = nan(nROI,1);
rows = cell(nROI,18);

for r = 1:nROI
    y = X(:,r);
    keep = isfinite(y);
    D = table(y(keep), categorical(group_num(keep)), C.cov_age(keep), C.cov_sex(keep), ...
        'VariableNames', {'y','group','age','sex'});

    if height(D) < 10, continue; end

    lm = fitlm(D, 'y ~ group + age + sex');
    A  = anova(lm,'summary');
    coef = lm.Coefficients;

    ridx = find(contains(lower(A.Properties.RowNames), 'group'), 1);
    if ~isempty(ridx), P0(r) = A.pValue(ridx); end

    rn = coef.Properties.RowNames;
    i2 = find(contains(rn,'group_2'),1);
    i3 = find(contains(rn,'group_3'),1);

    b2 = NaN; t2 = NaN; p2 = NaN;
    b3 = NaN; t3 = NaN; p3 = NaN;

    if ~isempty(i2), b2 = coef.Estimate(i2); t2 = coef.tStat(i2); p2 = coef.pValue(i2); P2(r)=p2; end
    if ~isempty(i3), b3 = coef.Estimate(i3); t3 = coef.tStat(i3); p3 = coef.pValue(i3); P3(r)=p3; end

    m1 = mean(y(group_num==1), 'omitnan');
    m2 = mean(y(group_num==2), 'omitnan');
    m3 = mean(y(group_num==3), 'omitnan');

    n1 = sum(isfinite(y(group_num==1)));
    n2 = sum(isfinite(y(group_num==2)));
    n3 = sum(isfinite(y(group_num==3)));

    rows(r,:) = {r, labels{r}, m1,m2,m3, n1,n2,n3, P0(r), b2,t2,p2, b3,t3,p3, NaN,NaN,NaN};
end

Q0 = fill_q(P0); Q2 = fill_q(P2); Q3 = fill_q(P3);

for r = 1:nROI
    if isempty(rows{r,1}), continue; end
    rows{r,16} = Q0(r); rows{r,17} = Q2(r); rows{r,18} = Q3(r);
    fprintf(fid,'%d\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', ...
        rows{r,1}, rows{r,2}, ...
        ns(rows{r,3}), ns(rows{r,4}), ns(rows{r,5}), ...
        rows{r,6}, rows{r,7}, rows{r,8}, ...
        ns(rows{r,9}), ns(rows{r,10}), ns(rows{r,11}), ns(rows{r,12}), ...
        ns(rows{r,13}), ns(rows{r,14}), ns(rows{r,15}), ...
        ns(rows{r,16}), ns(rows{r,17}), ns(rows{r,18}));
end
fclose(fid);

fid = fopen(out_fdr,'w');
fprintf(fid,'ROI\tLabel\tContrast\tBeta\tT\tP\tQ\n');
for r = 1:nROI
    if isfinite(Q2(r)) && Q2(r)<0.05
        fprintf(fid,'%d\t%s\tgroup2_vs_group1\t%s\t%s\t%s\t%s\n', ...
            r, labels{r}, ns(rows{r,10}), ns(rows{r,11}), ns(rows{r,12}), ns(Q2(r)));
    end
    if isfinite(Q3(r)) && Q3(r)<0.05
        fprintf(fid,'%d\t%s\tgroup3_vs_group1\t%s\t%s\t%s\t%s\n', ...
            r, labels{r}, ns(rows{r,13}), ns(rows{r,14}), ns(rows{r,15}), ns(Q3(r)));
    end
end
fclose(fid);

fprintf('[%s] 完成\n', metric_name);
end

function q = fill_q(p)
q = nan(size(p));
ok = isfinite(p);
if any(ok), q(ok) = mafdr(p(ok), 'BHFDR', true); end
end

function x = to_num(v)
if isnumeric(v)
    x = v;
else
    x = str2double(string(v));
end
if isempty(x), x = NaN; end
end

function s = ns(x)
if isnan(x), s='NaN'; else, s=sprintf('%.6g',x); end
end

function labels = read_labels(label_file, nROI)
labels = arrayfun(@(x) sprintf('ROI_%d',x), 1:nROI, 'UniformOutput', false)';
if ~exist(label_file,'file'), return; end
fid = fopen(label_file,'r');
k = 0;
while ~feof(fid) && k<nROI
    line = strtrim(fgetl(fid));
    if isempty(line), continue; end
    k = k + 1;
    labels{k} = line;
end
fclose(fid);
end
