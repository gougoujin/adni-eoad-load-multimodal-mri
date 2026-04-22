% Repository usage summary
% Description: Fit age/sex-adjusted group models for mean structural-functional coupling.
% Usage: Provide coupling summary and covariate table paths, then run the script/function.
% Outputs: Writes matched design tables, model coefficients, and ANOVA summaries.
% Note: Update input paths, toolboxes, and filenames for your local environment.

% Public repository version: update file paths, toolboxes, and local settings before running.
% This script/function was lightly sanitized for sharing and may require project-specific inputs.

clear; clc;

% ========= paths =========
coupling_file = '';
cov_file      = '';

out_match_txt = '';
out_lm_txt    = '';
out_anova_txt = '';
out_mat       = '';

% ========= read coupling =========
T = readtable(coupling_file, 'FileType', 'text', 'Delimiter', '\t');

% keep valid subjects only
T = T(~isnan(T.mean_coupling) & T.valid_roi > 0, :);

T.group = string(T.group);
T.sub   = string(T.sub);

T.group_num = nan(height(T),1);
for i = 1:height(T)
    tok = regexp(char(T.group(i)), 'group(\d+)', 'tokens', 'once');
    if ~isempty(tok)
        T.group_num(i) = str2double(tok{1});
    end
end

T.tag = string(T.group_num) + "_" + T.sub;

% ========= read covariate excel by column index =========
% col 1 = group
% col 3 = sub id
% col 4 = age
% col 5 = sex
raw = readcell(cov_file);

if size(raw,2) < 5
    error('Covariate file has fewer than 5 columns.');
end

raw = raw(2:end, :);  % remove header row

n = size(raw,1);

cov_group = nan(n,1);
cov_sub   = strings(n,1);
cov_age   = nan(n,1);
cov_sex_num = nan(n,1);
cov_sex_str = strings(n,1);

for i = 1:n
    % ----- group: col 1 -----
    g = raw{i,1};
    if isnumeric(g)
        cov_group(i) = g;
    elseif isstring(g) || ischar(g)
        tok = regexp(char(string(g)), '(\d+)', 'tokens', 'once');
        if ~isempty(tok)
            cov_group(i) = str2double(tok{1});
        end
    end

    % ----- sub id: col 3 -----
    s = raw{i,3};
    if isnumeric(s)
        if ~isnan(s)
            cov_sub(i) = sprintf('sub%03d', round(s));
        end
    else
        s = string(s);
        tok = regexp(char(s), '(\d+)', 'tokens', 'once');
        if ~isempty(tok)
            cov_sub(i) = sprintf('sub%03d', str2double(tok{1}));
        end
    end

    % ----- age: col 4 -----
    a = raw{i,4};
    if isnumeric(a)
        cov_age(i) = a;
    elseif isstring(a) || ischar(a)
        tmp = str2double(string(a));
        if ~isnan(tmp)
            cov_age(i) = tmp;
        end
    end

    % ----- sex: col 5 -----
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

% build cov table
C = table;
C.group_num = cov_group;
C.sub = cov_sub;
C.age = cov_age;

% decide whether sex is numeric or text
if sum(~isnan(cov_sex_num)) >= max(5, round(0.8*n))
    C.sex = cov_sex_num;
    sex_is_numeric = true;
else
    C.sex = categorical(cov_sex_str);
    sex_is_numeric = false;
end

C.tag = string(C.group_num) + "_" + C.sub;

% remove rows with bad id/group
keepC = ~isnan(C.group_num) & strlength(C.sub) > 0;
if sex_is_numeric
    keepC = keepC & ~isnan(C.sex);
else
    keepC = keepC & ~ismissing(C.sex);
end
keepC = keepC & ~isnan(C.age);

C = C(keepC,:);

% ========= match =========
[common_tag, ia, ib] = intersect(T.tag, C.tag, 'stable');

Tm = T(ia,:);
Cm = C(ib,:);

fprintf('Valid coupling subjects: %d\n', height(T));
fprintf('Valid covariate subjects: %d\n', height(C));
fprintf('Matched subjects: %d\n', numel(common_tag));

% ========= model table =========
D = table;
D.mean_coupling = Tm.mean_coupling;
D.group = categorical(Tm.group_num);
D.age = Cm.age;
D.sex = Cm.sex;

% ========= export matched table =========
M = table;
M.group = Tm.group;
M.group_num = Tm.group_num;
M.sub = Tm.sub;
M.mean_coupling = Tm.mean_coupling;
M.valid_roi = Tm.valid_roi;
M.age = Cm.age;
M.sex = Cm.sex;

writetable(M, out_match_txt, 'Delimiter', '\t', 'FileType', 'text');

% ========= fit model =========
mdl = fitlm(D, 'mean_coupling ~ group + age + sex');

% ---- extract coefficients into plain table ----
coef_raw = mdl.Coefficients;
coef_names = coef_raw.Properties.RowNames;

coef_out = table;
coef_out.Term = string(coef_names);
coef_out.Estimate = coef_raw{:, 'Estimate'};
coef_out.SE = coef_raw{:, 'SE'};
coef_out.tStat = coef_raw{:, 'tStat'};
coef_out.pValue = coef_raw{:, 'pValue'};

% ---- extract ANOVA into plain table ----
anova_raw = anova(mdl, 'summary');

anova_out = table;
anova_out.Term = string(anova_raw.Properties.RowNames);

% 兼容不同 MATLAB 版本的列名
anova_vars = string(anova_raw.Properties.VariableNames);

if any(anova_vars == "SumSq")
    anova_out.SumSq = anova_raw{:, "SumSq"};
elseif any(anova_vars == "Sum_Sq")
    anova_out.SumSq = anova_raw{:, "Sum_Sq"};
else
    anova_out.SumSq = nan(height(anova_raw),1);
end

if any(anova_vars == "DF")
    anova_out.DF = anova_raw{:, "DF"};
else
    anova_out.DF = nan(height(anova_raw),1);
end

if any(anova_vars == "MeanSq")
    anova_out.MeanSq = anova_raw{:, "MeanSq"};
elseif any(anova_vars == "Mean_Sq")
    anova_out.MeanSq = anova_raw{:, "Mean_Sq"};
else
    anova_out.MeanSq = nan(height(anova_raw),1);
end

if any(anova_vars == "F")
    anova_out.F = anova_raw{:, "F"};
else
    anova_out.F = nan(height(anova_raw),1);
end

if any(anova_vars == "pValue")
    anova_out.pValue = anova_raw{:, "pValue"};
elseif any(anova_vars == "p_Value")
    anova_out.pValue = anova_raw{:, "p_Value"};
else
    anova_out.pValue = nan(height(anova_raw),1);
end

% ---- write outputs ----
writetable(M, out_match_txt, 'Delimiter', '\t', 'FileType', 'text');
writetable(coef_out, out_lm_txt, 'Delimiter', '\t', 'FileType', 'text');
writetable(anova_out, out_anova_txt, 'Delimiter', '\t', 'FileType', 'text');

save(out_mat, 'mdl', 'D', 'M', 'Tm', 'Cm', 'coef_out', 'anova_out');

disp('===== MODEL COEFFICIENTS =====');
disp(coef_out);

disp('===== ANOVA TABLE =====');
disp(anova_out);

fprintf('\nOutput file/path/to/project', ...
    out_match_txt, out_lm_txt, out_anova_txt, out_mat);
