% Repository usage summary
% Description: Create a forest plot for representative ALFF/fALFF ROI effects.
% Usage: Run after the main ALFF/fALFF value table has been generated.
% Outputs: Exports figure images and a compact summary table.
% Note: Update input paths, toolboxes, and filenames for your local environment.

% Public repository version: update file paths, toolboxes, and local settings before running.
% This script/function was lightly sanitized for sharing and may require project-specific inputs.

function out_tif = figure_alff_falff_forest(values_file, out_dir)
% Figure 10: ALFF/fALFF representative ROI forest plot
% Uses Figure6_alff_falff_values.csv directly
%
% Output:
%   Figure10_alff_falff_forest.tif
%   Figure10_alff_falff_forest.png
%   Figure10_alff_falff_forest_table.csv
%
% Default:
%   values_file = ''
%   out_dir     = ''

    if nargin < 1 || isempty(values_file)
        values_file = '';
    end
    if nargin < 2 || isempty(out_dir)
        out_dir = '';
    end

    if ~isfile(values_file)
        error(['缺少文件: Figure6_alff_falff_values.csv\n' ...
               '原路径应为: /path/to/project' ...
               '当前路径: ' values_file]);
    end
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end

    T = readtable(values_file);
    T.Properties.VariableNames = matlab.lang.makeValidName(T.Properties.VariableNames);
    T = local_fix_names(T);

    req = {'group_num','age','sex'};
    for i = 1:numel(req)
        if ~ismember(req{i}, T.Properties.VariableNames)
            error('功能值表缺少列: %s', req{i});
        end
    end

    alff_rois = {'Hippocampus_R','Hippocampus_L','Thalamus_R','Caudate_R','Pallidum_L','Insula_L'};
    falff_rois = {'Cuneus_L','Cerebelum_Crus1_L','Occipital_Inf_R','Temporal_Mid_L','Temporal_Pole_Sup_R','Precuneus_L'};

    Roi = {};
    GroupType = {};
    Beta = [];
    SE = [];
    Tstat = [];
    Pval = [];
    CI_lo = [];
    CI_hi = [];
    Nvalid = [];

    % ALFF
    for i = 1:numel(alff_rois)
        vname = matlab.lang.makeValidName(['zALFF_' alff_rois{i}]);
        if ~ismember(vname, T.Properties.VariableNames)
            continue;
        end
        [b,se,tv,pv,nv] = local_fit_one(T, vname);
        Roi{end+1,1} = strrep(alff_rois{i}, '_', ' '); %#ok<AGROW>
        GroupType{end+1,1} = 'ALFF'; %#ok<AGROW>
        Beta(end+1,1) = b; %#ok<AGROW>
        SE(end+1,1) = se; %#ok<AGROW>
        Tstat(end+1,1) = tv; %#ok<AGROW>
        Pval(end+1,1) = pv; %#ok<AGROW>
        CI_lo(end+1,1) = b - 1.96*se; %#ok<AGROW>
        CI_hi(end+1,1) = b + 1.96*se; %#ok<AGROW>
        Nvalid(end+1,1) = nv; %#ok<AGROW>
    end

    % fALFF
    for i = 1:numel(falff_rois)
        vname = matlab.lang.makeValidName(['zfALFF_' falff_rois{i}]);
        if ~ismember(vname, T.Properties.VariableNames)
            continue;
        end
        [b,se,tv,pv,nv] = local_fit_one(T, vname);
        Roi{end+1,1} = strrep(falff_rois{i}, '_', ' '); %#ok<AGROW>
        GroupType{end+1,1} = 'fALFF'; %#ok<AGROW>
        Beta(end+1,1) = b; %#ok<AGROW>
        SE(end+1,1) = se; %#ok<AGROW>
        Tstat(end+1,1) = tv; %#ok<AGROW>
        Pval(end+1,1) = pv; %#ok<AGROW>
        CI_lo(end+1,1) = b - 1.96*se; %#ok<AGROW>
        CI_hi(end+1,1) = b + 1.96*se; %#ok<AGROW>
        Nvalid(end+1,1) = nv; %#ok<AGROW>
    end

    R = table(GroupType, Roi, Beta, SE, Tstat, Pval, CI_lo, CI_hi, Nvalid);

    % order by group then |beta|
    grp_order = zeros(height(R),1);
    grp_order(strcmp(R.GroupType,'ALFF')) = 1;
    grp_order(strcmp(R.GroupType,'fALFF')) = 2;
    R.grp_order = grp_order;
    R.absbeta = abs(R.Beta);
    R = sortrows(R, {'grp_order','absbeta'}, {'ascend','descend'});

    fig = figure('Color','w', 'Position', [80 80 1450 820]);
    ax = axes(fig); hold(ax, 'on');

    y = 1:height(R);
    colors = zeros(height(R),3);
    for i = 1:height(R)
        if strcmp(R.GroupType{i}, 'ALFF')
            colors(i,:) = [0.84 0.38 0.38];
        else
            colors(i,:) = [0.34 0.60 0.84];
        end
    end

    for i = 1:height(R)
        plot(ax, [R.CI_lo(i) R.CI_hi(i)], [y(i) y(i)], '-', 'Color', colors(i,:), 'LineWidth', 2);
        scatter(ax, R.Beta(i), y(i), 50, 'o', 'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', 'none');
        txt = sprintf('t=%.2f, p=%s', R.Tstat(i), local_p(R.Pval(i)));
        text(ax, R.CI_hi(i) + 0.02*range([min(R.CI_lo) max(R.CI_hi)]), y(i), txt, ...
            'FontName','Arial', 'FontSize',9, 'Color', [0.25 0.25 0.25], ...
            'HorizontalAlignment','left', 'VerticalAlignment','middle');
    end

    xline(ax, 0, '--', 'Color', [0.35 0.35 0.35], 'LineWidth', 1);

    % separator between ALFF and fALFF
    idx_sep = find(strcmp(R.GroupType,'ALFF'), 1, 'last');
    if ~isempty(idx_sep) && idx_sep < height(R)
        yline(ax, idx_sep + 0.5, '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
    end

    ax.YTick = y;
    ax.YTickLabel = R.Roi;
    ax.YDir = 'reverse';
    ax.Box = 'off';
    ax.TickDir = 'out';
    ax.FontName = 'Arial';
    ax.FontSize = 10;
    xlabel(ax, 'Adjusted \beta for LOAD vs EOAD (95% CI)', 'FontName','Arial', 'FontSize',11);
    title(ax, 'Figure 10. Representative ALFF/fALFF ROI effects', 'FontName','Arial', 'FontSize',13, 'FontWeight','bold');

    legend(ax, local_dummy_handles(), {'ALFF','fALFF'}, 'Location','best', 'Box','off', 'FontName','Arial', 'FontSize',10);

    out_tif = fullfile(out_dir, 'Figure10_alff_falff_forest.tif');
    out_png = fullfile(out_dir, 'Figure10_alff_falff_forest.png');
    exportgraphics(fig, out_tif, 'Resolution', 600, 'BackgroundColor','white', 'ContentType','image');
    exportgraphics(fig, out_png, 'Resolution', 300, 'BackgroundColor','white', 'ContentType','image');

    writetable(R(:, {'GroupType','Roi','Beta','SE','Tstat','Pval','CI_lo','CI_hi','Nvalid'}), ...
        fullfile(out_dir, 'Figure10_alff_falff_forest_table.csv'));

end

function [b,se,tv,pv,nv] = local_fit_one(T, vname)
    y = double(T.(vname));
    g = double(T.group_num);
    a = double(T.age);
    s = double(T.sex);

    ok = isfinite(y) & isfinite(g) & isfinite(a) & isfinite(s);
    nv = sum(ok);
    b = nan; se = nan; tv = nan; pv = nan;
    if nv < 12 || numel(unique(g(ok))) < 2
        return;
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
        return;
    end

    b = lm.Coefficients.Estimate(idx2);
    se = lm.Coefficients.SE(idx2);
    tv = lm.Coefficients.tStat(idx2);
    pv = lm.Coefficients.pValue(idx2);
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

function s = local_p(p)
    if ~isfinite(p)
        s = 'NA';
    elseif p < 0.001
        s = '<0.001';
    else
        s = sprintf('%.3f', p);
    end
end

function h = local_dummy_handles()
    h(1) = scatter(nan, nan, 50, 'o', 'MarkerFaceColor', [0.84 0.38 0.38], 'MarkerEdgeColor', 'none');
    h(2) = scatter(nan, nan, 50, 'o', 'MarkerFaceColor', [0.34 0.60 0.84], 'MarkerEdgeColor', 'none');
end
