% Repository usage summary
% Description: Run disease-only age-handling sensitivity analyses across principal positive domains.
% Usage: Run after the main coupling, GRETNA, and ALFF/fALFF result tables are available.
% Outputs: Outputs sensitivity summaries for quadratic-age and interaction models.
% Note: Update input paths, toolboxes, and filenames for your local environment.

% Public repository version: update file paths, toolboxes, and local settings before running.
% This script/function was lightly sanitized for sharing and may require project-specific inputs.

function run_patient_only_age_sensitivity_analyses()
% Disease-only age-handling sensitivity analyses for principal positive domains
%
% Models:
%   A) y ~ group + age + age^2 + sex
%   B) y ~ group + age + sex + group*age
%
% Disease-only sample:
%   group1 = EOAD
%   group2 = LOAD
%
% Domains analyzed:
%   1) SC-FC mean + 4 frontal ROIs
%   2) GRETNA representative graph metrics
%   3) Representative ALFF/fALFF ROIs
%
% Output:
%   /path/to/project
%       age2_summary_all.xlsx
%       interaction_summary_all.xlsx
%       stability_summary_all.xlsx

    clear; clc;

    %% ===================== paths =====================
    out_dir = '';
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end

    cov_file = '';
    aal_atlas_file = '';
    scfc_mat = '';

    gretna_root = '';
    dpabi_group_roots = { ...
        '', ...   % group1 = EOAD
        '', ...  % group2 = LOAD
        '' ...      % group3 = HC
        };

    if ~exist(cov_file, 'file')
        error('Missing covariate file: %s', cov_file);
    end
    if ~exist(aal_atlas_file, 'file')
        error('Missing atlas file: %s', aal_atlas_file);
    end

    %% ===================== fixed ROI indices =====================
    % GRETNA representative ROIs from your current main text / Table 2
    % 115 = Vermis_9
    % 79  = Heschl_L
    % 81  = Temporal_Sup_L
    % 97  = Cerebelum_4_5_L
    % 11  = Frontal_Inf_Oper_L
    % 23  = Frontal_Sup_Medial_L
    gretna_roi_idx = [115, 79, 81, 97, 11, 23];

    % ALFF representative ROIs
    % 38 = Hippocampus_R
    % 37 = Hippocampus_L
    % 78 = Thalamus_R
    % 72 = Caudate_R
    % 75 = Pallidum_L
    % 29 = Insula_L
    alff_roi_idx = [38, 37, 78, 72, 75, 29];

    % fALFF representative ROIs
    % 45 = Cuneus_L
    % 91 = Cerebelum_Crus1_L
    % 54 = Occipital_Inf_R
    % 85 = Temporal_Mid_L
    % 84 = Temporal_Pole_Sup_R
    % 67 = Precuneus_L
    falff_roi_idx = [45, 91, 54, 85, 84, 67];

    %% ===================== readable ROI names =====================
    gretna_roi_names = { ...
        'Vermis_9', ...
        'Heschl_L', ...
        'Temporal_Sup_L', ...
        'Cerebelum_4_5_L', ...
        'Frontal_Inf_Oper_L', ...
        'Frontal_Sup_Medial_L' ...
        };

    alff_roi_names = { ...
        'Hippocampus_R', ...
        'Hippocampus_L', ...
        'Thalamus_R', ...
        'Caudate_R', ...
        'Pallidum_L', ...
        'Insula_L' ...
        };

    falff_roi_names = { ...
        'Cuneus_L', ...
        'Cerebelum_Crus1_L', ...
        'Occipital_Inf_R', ...
        'Temporal_Mid_L', ...
        'Temporal_Pole_Sup_R', ...
        'Precuneus_L' ...
        };

    %% ===================== read covariates =====================
    Tcov = read_covariates(cov_file);

    %% ===================== containers =====================
    all_age2 = table();
    all_int  = table();

    %% =========================================================
    %% 1) SC-FC mean + 4 frontal ROIs
    %% =========================================================
    fprintf('\n========== SC-FC ==========\n');

    if exist(scfc_mat, 'file')
        S = load(scfc_mat);

        if ~isfield(S, 'roi_coupling') || ~isfield(S, 'group_num') || ~isfield(S, 'sub_id')
            error('SC-FC file missing roi_coupling/group_num/sub_id');
        end

        X = double(S.roi_coupling);     % ROI x subject
        group_num = double(S.group_num(:));
        sub_id = string(S.sub_id(:));

        Tsc = table;
        Tsc.group_num = group_num;
        Tsc.sub = sub_id;
        Tsc.tag = string(Tsc.group_num) + "_" + Tsc.sub;

        [Tf, Xf] = merge_subjects_keep_order(Tcov, Tsc, X);

        keep = Tf.group_num==1 | Tf.group_num==2;
        Tf = Tf(keep,:);
        Xf = Xf(:, keep);

        mean_scfc = mean(Xf, 1, 'omitnan')';

        roi_defs = { ...
            'SCFC_Mean',               NaN; ...
            'SCFC_Frontal_Mid_R',      8;  ...
            'SCFC_Frontal_Mid_Orb_R',  10; ...
            'SCFC_Frontal_Inf_Tri_R',  14; ...
            'SCFC_Frontal_Sup_Medial_L',23  ...
            };

        for i = 1:size(roi_defs,1)
            feat = roi_defs{i,1};
            roi  = roi_defs{i,2};

            if isnan(roi)
                y = mean_scfc;
            else
                y = Xf(roi,:)';
            end

            T = build_base_table(Tf, y);
            [A, B] = fit_two_models(T, 'SCFC', feat);
            all_age2 = [all_age2; A]; %#ok<AGROW>
            all_int  = [all_int;  B]; %#ok<AGROW>
        end
    else
        warning('SC-FC file not found: %s', scfc_mat);
    end

    %% =========================================================
    %% 2) GRETNA representative metrics
    %% =========================================================
    fprintf('\n========== GRETNA ==========\n');

    gretna_specs = { ...
        'BetweennessCentrality', 'aBc.txt'; ...
        'DegreeCentrality',      'aDc.txt'; ...
        'NodalEfficiency',       'aNe.txt'; ...
        'NodalLocalEfficiency',  'aNLe.txt' ...
        };

    Tg = table();

    for m = 1:size(gretna_specs,1)
        metric_dir_name = gretna_specs{m,1};
        metric_file_name = gretna_specs{m,2};

        fprintf('Reading %s (%s)\n', metric_dir_name, metric_file_name);

        Xmetric = table();

        for g = 1:3
            gfile = fullfile(gretna_root, metric_dir_name, sprintf('group%d', g), metric_file_name);
            if ~exist(gfile, 'file')
                warning('Missing GRETNA file: %s', gfile);
                continue;
            end

            A = readmatrix(gfile);
            A = double(A);
            if isempty(A)
                continue;
            end

            subs_g = Tcov.sub(Tcov.group_num == g);
            subs_g = sort_sub_ids(unique(subs_g));
            nsub = numel(subs_g);

            sa = size(A);
            if sa(2) == 116
                B = A;
            elseif sa(1) == 116
                B = A';
            else
                warning('Cannot recognize GRETNA matrix orientation: %s size=[%d %d]', gfile, sa(1), sa(2));
                continue;
            end

            n_use = min(size(B,1), nsub);
            B = B(1:n_use, :);
            subs_use = subs_g(1:n_use);

            Ttmp = table;
            Ttmp.group_num = repmat(g, n_use, 1);
            Ttmp.sub = subs_use(:);
            Ttmp.tag = string(Ttmp.group_num) + "_" + Ttmp.sub;

            for r = 1:numel(gretna_roi_idx)
                roi = gretna_roi_idx(r);
                vname = safe_var_name([metric_dir_name '_' gretna_roi_names{r}], sprintf('%s_ROI%d', metric_dir_name, roi));
                Ttmp.(vname) = B(:, roi);
            end

            if isempty(Xmetric)
                Xmetric = Ttmp;
            else
                Xmetric = [Xmetric; Ttmp]; %#ok<AGROW>
            end
        end

        if isempty(Xmetric)
            continue;
        end

        if isempty(Tg)
            Tg = Xmetric;
        else
            Tg = outerjoin(Tg, Xmetric, 'Keys', 'tag', 'MergeKeys', true);
        end
    end

    if ~isempty(Tg)
        [Tg2, ~] = merge_subjects_keep_order(Tcov, Tg, []);
        keep = Tg2.group_num==1 | Tg2.group_num==2;
        Tg2 = Tg2(keep,:);

        use_vars = {};
        all_vars = Tg2.Properties.VariableNames;

        for m = 1:size(gretna_specs,1)
            metric_dir_name = gretna_specs{m,1};
            for r = 1:numel(gretna_roi_idx)
                roi = gretna_roi_idx(r);
                vname = safe_var_name([metric_dir_name '_' gretna_roi_names{r}], sprintf('%s_ROI%d', metric_dir_name, roi));
                if ismember(vname, all_vars)
                    use_vars{end+1} = vname; %#ok<AGROW>
                end
            end
        end

        use_vars = unique(use_vars, 'stable');

        for i = 1:numel(use_vars)
            vn = use_vars{i};
            y = double(Tg2.(vn));
            T = build_base_table(Tg2, y);
            [A, B] = fit_two_models(T, 'Graph', vn);
            all_age2 = [all_age2; A]; %#ok<AGROW>
            all_int  = [all_int;  B]; %#ok<AGROW>
        end
    else
        warning('No GRETNA subject-level table assembled.');
    end

    %% =========================================================
    %% 3) ALFF / fALFF representative ROIs
    %% =========================================================
    fprintf('\n========== ALFF / fALFF ==========\n');

    atlas = double(read_nifti_auto(aal_atlas_file));
    Talff = table();

    for g = 1:3
        group_root = dpabi_group_roots{g};
        if ~exist(group_root, 'dir')
            warning('Missing DPABI group root: %s', group_root);
            continue;
        end

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
                if exist(gz, 'file')
                    alff_file = gz;
                else
                    alff_file = '';
                end
            end

            if ~exist(falff_file, 'file')
                gz = [falff_file '.gz'];
                if exist(gz, 'file')
                    falff_file = gz;
                else
                    falff_file = '';
                end
            end

            if isempty(alff_file) && isempty(falff_file)
                continue;
            end

            arow = nan(1, numel(alff_roi_idx));
            frow = nan(1, numel(falff_roi_idx));

            if ~isempty(alff_file)
                Aimg = double(read_nifti_auto(alff_file));
                if isequal(size(Aimg), size(atlas))
                    for r = 1:numel(alff_roi_idx)
                        roi = alff_roi_idx(r);
                        mask = atlas == roi;
                        arow(r) = mean_valid(Aimg(mask));
                    end
                end
            end

            if ~isempty(falff_file)
                Fimg = double(read_nifti_auto(falff_file));
                if isequal(size(Fimg), size(atlas))
                    for r = 1:numel(falff_roi_idx)
                        roi = falff_roi_idx(r);
                        mask = atlas == roi;
                        frow(r) = mean_valid(Fimg(mask));
                    end
                end
            end

            rows(end+1,:) = {g, sid, sprintf('%d_%s', g, sid)}; %#ok<AGROW>
            ALFF_vals(end+1,:) = arow; %#ok<AGROW>
            fALFF_vals(end+1,:) = frow; %#ok<AGROW>
        end

        if isempty(rows)
            continue;
        end

        X = table;
        X.group_num = cell2mat(rows(:,1));
        X.sub = string(rows(:,2));
        X.tag = string(rows(:,3));

        for r = 1:numel(alff_roi_idx)
            roi = alff_roi_idx(r);
            vname = safe_var_name(['zALFF_' alff_roi_names{r}], sprintf('zALFF_ROI%d', roi));
            X.(vname) = ALFF_vals(:,r);
        end

        for r = 1:numel(falff_roi_idx)
            roi = falff_roi_idx(r);
            vname = safe_var_name(['zfALFF_' falff_roi_names{r}], sprintf('zfALFF_ROI%d', roi));
            X.(vname) = fALFF_vals(:,r);
        end

        if isempty(Talff)
            Talff = X;
        else
            Talff = [Talff; X]; %#ok<AGROW>
        end
    end

    if ~isempty(Talff)
        [Talff2, ~] = merge_subjects_keep_order(Tcov, Talff, []);
        keep = Talff2.group_num==1 | Talff2.group_num==2;
        Talff2 = Talff2(keep,:);

        all_vars = Talff2.Properties.VariableNames;
        use_vars = {};

        for r = 1:numel(alff_roi_idx)
            roi = alff_roi_idx(r);
            vname = safe_var_name(['zALFF_' alff_roi_names{r}], sprintf('zALFF_ROI%d', roi));
            if ismember(vname, all_vars)
                use_vars{end+1} = vname; %#ok<AGROW>
            end
        end

        for r = 1:numel(falff_roi_idx)
            roi = falff_roi_idx(r);
            vname = safe_var_name(['zfALFF_' falff_roi_names{r}], sprintf('zfALFF_ROI%d', roi));
            if ismember(vname, all_vars)
                use_vars{end+1} = vname; %#ok<AGROW>
            end
        end

        use_vars = unique(use_vars, 'stable');

        for i = 1:numel(use_vars)
            vn = use_vars{i};
            y = double(Talff2.(vn));
            T = build_base_table(Talff2, y);
            [A, B] = fit_two_models(T, 'LocalFunction', vn);
            all_age2 = [all_age2; A]; %#ok<AGROW>
            all_int  = [all_int;  B]; %#ok<AGROW>
        end
    else
        warning('No ALFF/fALFF subject-level table assembled.');
    end

    %% ===================== write outputs =====================
    age2_file = fullfile(out_dir, 'age2_summary_all.xlsx');
    int_file  = fullfile(out_dir, 'interaction_summary_all.xlsx');
    stab_file = fullfile(out_dir, 'stability_summary_all.xlsx');

    writetable(all_age2, age2_file, 'FileType', 'spreadsheet');
    writetable(all_int,  int_file,  'FileType', 'spreadsheet');

    S = build_stability_summary(all_age2, all_int);
    writetable(S, stab_file, 'FileType', 'spreadsheet');

    fprintf('\nSave/path/to/project', age2_file, int_file, stab_file);
end

%% ===================== model fitting =====================

function [A, B] = fit_two_models(T, domain_name, feature_name)
    T = T(isfinite(T.y) & isfinite(T.age) & isfinite(T.sex), :);

    lm_age2 = fitlm(T, 'y ~ group + age_c + age_c2 + sex');
    lm_int  = fitlm(T, 'y ~ group + age_c + sex + group:age_c');

    A = extract_age2_summary(lm_age2, domain_name, feature_name);
    B = extract_interaction_summary(lm_int, domain_name, feature_name);

    fprintf('\n[%s | %s]\n', domain_name, feature_name);
    fprintf('  Age^2 model: group p = %.4g; age^2 p = %.4g\n', A.P_group, A.P_age2);
    fprintf('  Interaction model: group p = %.4g; interaction p = %.4g\n', B.P_group, B.P_interaction);
end

function A = extract_age2_summary(lm, domain_name, feature_name)
    coefTbl = lm.Coefficients;
    rn = string(coefTbl.Properties.RowNames);
    ci_all = coefCI(lm);
    [F_model, P_model] = extract_model_fp(lm);

    idx_group = find(rn=="group_LOAD", 1);
    idx_age   = find(rn=="age_c", 1);
    idx_age2  = find(rn=="age_c2", 1);
    idx_sex   = find(rn=="sex", 1);

    A = table;
    A.Domain  = string(domain_name);
    A.Feature = string(feature_name);
    A.Model   = "Age2_model";
    A.N       = height(lm.Variables);
    A.F_model = F_model;
    A.P_model = P_model;

    A.Beta_group = get_est(coefTbl, idx_group);
    A.T_group    = get_t(coefTbl, idx_group);
    A.P_group    = get_p(coefTbl, idx_group);
    A.CIlo_group = get_ci(ci_all, idx_group, 1);
    A.CIhi_group = get_ci(ci_all, idx_group, 2);

    A.Beta_age = get_est(coefTbl, idx_age);
    A.P_age    = get_p(coefTbl, idx_age);

    A.Beta_age2 = get_est(coefTbl, idx_age2);
    A.P_age2    = get_p(coefTbl, idx_age2);

    A.Beta_sex = get_est(coefTbl, idx_sex);
    A.P_sex    = get_p(coefTbl, idx_sex);
end

function B = extract_interaction_summary(lm, domain_name, feature_name)
    coefTbl = lm.Coefficients;
    rn = string(coefTbl.Properties.RowNames);
    ci_all = coefCI(lm);
    [F_model, P_model] = extract_model_fp(lm);

    idx_group = find(rn=="group_LOAD", 1);
    idx_age   = find(rn=="age_c", 1);
    idx_int   = find(contains(rn, "group_LOAD:age_c") | contains(rn, "age_c:group_LOAD"), 1);
    idx_sex   = find(rn=="sex", 1);

    B = table;
    B.Domain  = string(domain_name);
    B.Feature = string(feature_name);
    B.Model   = "Interaction_model";
    B.N       = height(lm.Variables);
    B.F_model = F_model;
    B.P_model = P_model;

    B.Beta_group = get_est(coefTbl, idx_group);
    B.T_group    = get_t(coefTbl, idx_group);
    B.P_group    = get_p(coefTbl, idx_group);
    B.CIlo_group = get_ci(ci_all, idx_group, 1);
    B.CIhi_group = get_ci(ci_all, idx_group, 2);

    B.Beta_age = get_est(coefTbl, idx_age);
    B.P_age    = get_p(coefTbl, idx_age);

    B.Beta_interaction = get_est(coefTbl, idx_int);
    B.T_interaction    = get_t(coefTbl, idx_int);
    B.P_interaction    = get_p(coefTbl, idx_int);

    B.Beta_sex = get_est(coefTbl, idx_sex);
    B.P_sex    = get_p(coefTbl, idx_sex);
end

function [F_model, P_model] = extract_model_fp(lm)
    F_model = NaN;
    P_model = NaN;
    try
        a = anova(lm, 'summary');
        rn = string(a.Properties.RowNames);
        idxm = find(strcmpi(rn, 'Model'), 1);
        if isempty(idxm)
            idxm = 1;
        end
        if ismember('F', a.Properties.VariableNames)
            F_model = a.F(idxm);
        end
        if ismember('pValue', a.Properties.VariableNames)
            P_model = a.pValue(idxm);
        end
    catch
    end
end

%% ===================== summary builder =====================

function S = build_stability_summary(all_age2, all_int)
    keys = unique(string(all_age2.Domain) + "||" + string(all_age2.Feature));
    S = table();

    for i = 1:numel(keys)
        key = keys(i);

        k1 = string(all_age2.Domain) + "||" + string(all_age2.Feature);
        k2 = string(all_int.Domain)  + "||" + string(all_int.Feature);

        A = all_age2(k1==key,:);
        B = all_int(k2==key,:);

        tmp = table;
        tmp.Domain = string(A.Domain(1));
        tmp.Feature = string(A.Feature(1));

        tmp.Age2_group_p = A.P_group(1);
        tmp.Age2_age2_p  = A.P_age2(1);

        tmp.Interaction_group_p = B.P_group(1);
        tmp.Interaction_p       = B.P_interaction(1);

        tmp.Direction_consistent = sign(A.Beta_group(1)) == sign(B.Beta_group(1));
        tmp.Stability_flag = stability_flag(A.P_group(1), B.P_group(1), B.P_interaction(1));

        S = [S; tmp]; %#ok<AGROW>
    end
end

function flag = stability_flag(p_age2_group, p_int_group, p_int)
    if isfinite(p_age2_group) && p_age2_group < 0.05 && ...
       isfinite(p_int_group)  && p_int_group  < 0.05 && ...
       (~isfinite(p_int) || p_int >= 0.05)
        flag = "robust";
    elseif (isfinite(p_age2_group) && p_age2_group < 0.05) || ...
           (isfinite(p_int_group)  && p_int_group  < 0.05)
        flag = "partially_stable";
    else
        flag = "attenuated";
    end
end

%% ===================== helpers =====================

function T = build_base_table(Tcov_sub, y)
    T = table;
    T.group_num = Tcov_sub.group_num;
    T.group = categorical(Tcov_sub.group_num, [1 2], {'EOAD','LOAD'});
    T.sub = Tcov_sub.sub;
    T.age = Tcov_sub.age;
    T.sex = Tcov_sub.sex;
    T.age_c = T.age - mean(T.age, 'omitnan');
    T.age_c2 = T.age_c.^2;
    T.y = y(:);
end

function [Tout, Xout] = merge_subjects_keep_order(Tcov, Tsub, X)
    [~, ia, ib] = intersect(Tcov.tag, Tsub.tag, 'stable');

    % keep covariates in stable order
    Tout = Tcov(ia,:);

    % when Tsub is a subject-level table with feature columns,
    % append those feature columns back to Tout
    vars_sub = Tsub.Properties.VariableNames;
    vars_keep = setdiff(vars_sub, {'group_num','sub','tag'});
    if ~isempty(vars_keep)
        Tout = [Tout, Tsub(ib, vars_keep)];
    end

    if isempty(X)
        Xout = [];
    else
        Xout = X(:, ib);
    end
end

function T = read_covariates(cov_file)
    raw = readcell(cov_file);
    raw = raw(2:end, :);
    n = size(raw,1);

    group_num = nan(n,1);
    sub_id = cell(n,1);
    age = nan(n,1);
    sex = nan(n,1);

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

    T = table;
    T.group_num = group_num;
    T.sub = string(sub_id);
    T.age = age;
    T.sex = sex;
    T.tag = string(T.group_num) + "_" + T.sub;
    T = T(isfinite(T.group_num) & strlength(T.sub)>0 & isfinite(T.age) & isfinite(T.sex), :);
end

function ids = sort_sub_ids(ids)
    ids = string(ids);
    nums = nan(size(ids));
    for i = 1:numel(ids)
        tok = regexp(char(ids(i)), '(\d+)', 'tokens', 'once');
        if ~isempty(tok)
            nums(i) = str2double(tok{1});
        end
    end
    [~, ord] = sort(nums);
    ids = ids(ord);
end

function out = safe_var_name(s, fallback)
    try
        out = matlab.lang.makeValidName(s);
        if isempty(out)
            out = fallback;
        end
    catch
        out = fallback;
    end
end

function x = to_num(v)
    if isnumeric(v)
        if isempty(v) || isnan(v)
            x = NaN;
        else
            x = double(v(1));
        end
    else
        x = str2double(string(v));
        if isnan(x)
            tok = regexp(char(string(v)), '[-+]?\d*\.?\d+', 'match', 'once');
            if isempty(tok)
                x = NaN;
            else
                x = str2double(tok);
            end
        end
    end
end

function y = mean_valid(x)
    x = x(isfinite(x));
    if isempty(x)
        y = NaN;
    else
        y = mean(x);
    end
end

function x = read_nifti_auto(fn)
    if endsWith(fn, '.gz', 'IgnoreCase', true)
        td = tempname;
        mkdir(td);
        gunzip(fn, td);
        d = dir(fullfile(td, '*.nii'));
        if isempty(d)
            error('gunzip failed for %s', fn);
        end
        nii = fullfile(d(1).folder, d(1).name);
        x = niftiread(nii);
        try, delete(nii); end %#ok<TRYNC>
        try, rmdir(td); end %#ok<TRYNC>
    else
        x = niftiread(fn);
    end
end

function x = get_est(tbl, idx)
    if isempty(idx)
        x = NaN;
    else
        x = tbl.Estimate(idx);
    end
end

function x = get_t(tbl, idx)
    if isempty(idx)
        x = NaN;
    else
        x = tbl.tStat(idx);
    end
end

function x = get_p(tbl, idx)
    if isempty(idx)
        x = NaN;
    else
        x = tbl.pValue(idx);
    end
end

function x = get_ci(ci, idx, col)
    if isempty(idx)
        x = NaN;
    else
        x = ci(idx,col);
    end
end
