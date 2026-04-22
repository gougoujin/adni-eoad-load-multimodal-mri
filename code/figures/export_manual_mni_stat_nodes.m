% Repository usage summary
% Description: Export BrainNet node files from manually curated MNI coordinates and statistics.
% Usage: Run after preparing a coordinate template and functional value table.
% Outputs: Writes node files, info tables, and range suggestions.
% Note: Update input paths, toolboxes, and filenames for your local environment.

% Public repository version: update file paths, toolboxes, and local settings before running.
% This script/function was lightly sanitized for sharing and may require project-specific inputs.

function export_manual_mni_stat_nodes(coord_file, functional_values_file, out_dir)
% v4: use already exported Figure6_alff_falff_values.csv directly
% This avoids any re-extraction / re-matching mismatch.
%
% Inputs:
%   coord_file            manual_mni_coords_template_filled.xlsx/csv
%   functional_values_file Figure6_alff_falff_values.csv
%   out_dir               output folder
%
% Outputs:
%   ALFF_manual_mni_stat_v4.node
%   fALFF_manual_mni_stat_v4.node
%   ALFF_manual_mni_stat_v4_info.csv
%   fALFF_manual_mni_stat_v4_info.csv
%   manual_mni_stat_v4_range_suggestion.txt
%
% Color = Beta_LOAD_vs_EOAD
% Size  = scaled -log10(P_LOAD_vs_EOAD)

    if nargin < 1 || isempty(coord_file)
        coord_file = '';
    end
    if nargin < 2 || isempty(functional_values_file)
        functional_values_file = '';
    end
    if nargin < 3 || isempty(out_dir)
        out_dir = '';
    end

    if ~isfile(coord_file)
        error(['缺少文件: manual_mni_coords_template_filled.xlsx/csv\n原路径应为: /path/to/project ' coord_file]);
    end
    if ~isfile(functional_values_file)
        error(['缺少文件: Figure6_alff_falff_values.csv\n原路径应为: /path/to/project ' functional_values_file]);
    end
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end

    % read coordinate template
    [~,~,ext] = fileparts(coord_file);
    if strcmpi(ext, '.xlsx') || strcmpi(ext, '.xls')
        C = readtable(coord_file);
    else
        C = readtable(coord_file);
    end
    C.Properties.VariableNames = matlab.lang.makeValidName(C.Properties.VariableNames);

    reqC = {'Group','ROI','X','Y','Z'};
    for i = 1:numel(reqC)
        if ~ismember(reqC{i}, C.Properties.VariableNames)
            error('坐标模板缺少列: %s', reqC{i});
        end
    end

    V = readtable(functional_values_file);
    V.Properties.VariableNames = matlab.lang.makeValidName(V.Properties.VariableNames);

    % tolerate join-renamed names from previous scripts
    V = local_fix_names(V);

    reqV = {'group_num','age','sex'};
    for i = 1:numel(reqV)
        if ~ismember(reqV{i}, V.Properties.VariableNames)
            error('功能值表缺少列: %s', reqV{i});
        end
    end

    Out = C;
    Out.VariableUsed = strings(height(C),1);
    Out.N_valid = zeros(height(C),1);
    Out.N_group1 = zeros(height(C),1);
    Out.N_group2 = zeros(height(C),1);
    Out.N_group3 = zeros(height(C),1);
    Out.Diagnostic = strings(height(C),1);
    Out.Beta_LOAD_vs_EOAD = nan(height(C),1);
    Out.T_LOAD_vs_EOAD = nan(height(C),1);
    Out.P_LOAD_vs_EOAD = nan(height(C),1);
    Out.ModelP = nan(height(C),1);
    Out.NodeColor = nan(height(C),1);
    Out.NodeSize = nan(height(C),1);

    for i = 1:height(Out)
        roi_name = string(Out.ROI(i));
        grp_name = string(Out.Group(i));

        if strcmpi(grp_name, 'ALFF')
            prefix = "zALFF_";
        else
            prefix = "zfALFF_";
        end

        vname = matlab.lang.makeValidName(prefix + roi_name);
        if ~ismember(vname, V.Properties.VariableNames)
            % fallback fuzzy match
            allvn = string(V.Properties.VariableNames);
            hit = find(strcmpi(allvn, vname), 1);
            if isempty(hit)
                hit = find(contains(allvn, roi_name), 1);
            end
            if isempty(hit)
                Out.Diagnostic(i) = "variable missing";
                continue;
            else
                vname = allvn(hit);
            end
        end

        y = double(V.(vname));
        g = double(V.group_num);
        a = double(V.age);
        s = double(V.sex);

        ok = isfinite(y) & isfinite(g) & isfinite(a) & isfinite(s);
        Out.VariableUsed(i) = vname;
        Out.N_valid(i) = sum(ok);
        Out.N_group1(i) = sum(ok & g == 1);
        Out.N_group2(i) = sum(ok & g == 2);
        Out.N_group3(i) = sum(ok & g == 3);

        if sum(ok) < 12
            Out.Diagnostic(i) = "too few valid";
            continue;
        end
        if numel(unique(g(ok))) < 2
            Out.Diagnostic(i) = "fewer than 2 groups";
            continue;
        end
        if sum(ok & g == 1) < 3 || sum(ok & g == 2) < 3
            Out.Diagnostic(i) = "too few EOAD/LOAD";
            continue;
        end

        D = table(y(ok), categorical(g(ok)), a(ok), s(ok), ...
            'VariableNames', {'y','group','age','sex'});
        lm = fitlm(D, 'y ~ group + age + sex');

        rn = string(lm.Coefficients.Properties.RowNames);
        idx2 = find(rn == "group_2", 1);
        if isempty(idx2)
            idx2 = find(contains(rn, "group_2"), 1);
        end
        if isempty(idx2)
            Out.Diagnostic(i) = "group_2 coef missing";
            continue;
        end

        Out.Beta_LOAD_vs_EOAD(i) = lm.Coefficients.Estimate(idx2);
        Out.T_LOAD_vs_EOAD(i) = lm.Coefficients.tStat(idx2);
        Out.P_LOAD_vs_EOAD(i) = lm.Coefficients.pValue(idx2);
        Out.NodeColor(i) = Out.Beta_LOAD_vs_EOAD(i);
        Out.NodeSize(i) = 2.5 + 2.0 * min(3, -log10(max(Out.P_LOAD_vs_EOAD(i), realmin)));
        Out.Diagnostic(i) = "ok";

        aov = anova(lm, 'summary');
        if ismember('Term', aov.Properties.VariableNames) && ismember('pValue', aov.Properties.VariableNames)
            idxm = strcmp(string(aov.Term), 'Model');
            if any(idxm)
                Out.ModelP(i) = aov.pValue(find(idxm,1));
            elseif ~isempty(aov.pValue)
                Out.ModelP(i) = aov.pValue(1);
            end
        end
    end

    local_write_group(Out(strcmpi(string(Out.Group), 'ALFF'), :), fullfile(out_dir, 'ALFF_manual_mni_stat_v4'));
    local_write_group(Out(strcmpi(string(Out.Group), 'fALFF'), :), fullfile(out_dir, 'fALFF_manual_mni_stat_v4'));

    writetable(Out, fullfile(out_dir, 'manual_mni_stat_v4_full_diagnostics.csv'));

    beta_all = Out.NodeColor(isfinite(Out.NodeColor));
    if ~isempty(beta_all)
        bmin = min(beta_all);
        bmax = max(beta_all);
        babs = max(abs(beta_all));
        fid = fopen(fullfile(out_dir, 'manual_mni_stat_v4_range_suggestion.txt'), 'w');
        fprintf(fid, 'Suggested color ranges for BrainNe/path/to/project');
        fprintf(fid, 'Raw range:   min = %.6f, max = %.6f\n', bmin, bmax);
        fprintf(fid, 'Symmetric range (recommended): min = %.6f, max = %.6f\n', -babs, babs);
        fprintf(fid, '\nColor = Beta_LOAD_vs_EOAD\n');
        fprintf(fid, 'Size  = scaled -log10(P_LOAD_vs_EOAD)\n');
        fclose(fid);
    end
end

function T = local_fix_names(T)
    vn = string(T.Properties.VariableNames);
    if ~ismember("group_num", vn)
        cands = vn(contains(vn, "group_num"));
        if ~isempty(cands), T.group_num = T.(cands(1)); end
    end
    if ~ismember("age", vn)
        cands = vn(strcmpi(vn,"age_T") | strcmpi(vn,"age_Xall") | contains(vn, "age"));
        if ~isempty(cands), T.age = T.(cands(1)); end
    end
    if ~ismember("sex", vn)
        cands = vn(strcmpi(vn,"sex_T") | strcmpi(vn,"sex_Xall") | contains(vn, "sex"));
        if ~isempty(cands), T.sex = T.(cands(1)); end
    end
end

function local_write_group(T, out_prefix)
    writetable(T, [out_prefix '_info.csv']);
    fid = fopen([out_prefix '.node'], 'w');
    for i = 1:height(T)
        c = T.NodeColor(i); if ~isfinite(c), c = 0; end
        s = T.NodeSize(i); if ~isfinite(s), s = 4.0; end
        fprintf(fid, '%.4f\t%.4f\t%.4f\t%.6f\t%.4f\t%s\n', ...
            T.X(i), T.Y(i), T.Z(i), c, s, string(T.ROI(i)));
    end
    fclose(fid);
end
