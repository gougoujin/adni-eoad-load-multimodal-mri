% Repository usage summary
% Description: Create a dumbbell/forest plot for clinical association results.
% Usage: Run after AD-only clinical association tables have been generated.
% Outputs: Exports figure images and a plotting table.
% Note: Update input paths, toolboxes, and filenames for your local environment.

% Public repository version: update file paths, toolboxes, and local settings before running.
% This script/function was lightly sanitized for sharing and may require project-specific inputs.

function out_tif = figure_clinical_association_dumbbell(clin_file, out_dir, max_scales)
% Figure 12: clinical association dumbbell/forest plot
% Uses AD-only clinical association table from your pipeline.
%
% Inputs:
%   clin_file   = /path/to/project
%   out_dir     = /path/to/project
%   max_scales  = max number of scales to show (default 4, chosen by smallest p)
%
% Outputs:
%   Figure12_clinical_dumbbell.tif
%   Figure12_clinical_dumbbell.png
%   Figure12_clinical_dumbbell_table.csv

    if nargin < 1 || isempty(clin_file)
        clin_file = '';
    end
    if nargin < 2 || isempty(out_dir)
        out_dir = '';
    end
    if nargin < 3 || isempty(max_scales)
        max_scales = 4;
    end

    if ~isfile(clin_file)
        error(['缺少文件: clinical_association_AD_only.txt\n' ...
               '原路径应为: /path/to/project' ...
               '当前路径: ' clin_file]);
    end
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end

    T = readtable(clin_file, 'FileType','text', 'Delimiter','\t');
    T.Properties.VariableNames = matlab.lang.makeValidName(T.Properties.VariableNames);

    req = {'Imaging','Scale','Beta_raw','T_raw','P_raw','Q_raw'};
    for i = 1:numel(req)
        if ~ismember(req{i}, T.Properties.VariableNames)
            error('临床关联结果表缺少列: %s', req{i});
        end
    end

    % choose most informative scales by smallest raw p across all imaging vars
    scales = unique(string(T.Scale), 'stable');
    bestp = nan(numel(scales),1);
    for i = 1:numel(scales)
        idx = string(T.Scale) == scales(i);
        bestp(i) = min(T.P_raw(idx), [], 'omitnan');
    end
    [~, ord] = sort(bestp, 'ascend', 'MissingPlacement','last');
    scales = scales(ord);
    scales = scales(1:min(max_scales, numel(scales)));

    imaging_order = {'mean_coupling','ROI8','ROI10','ROI14','ROI23'};
    imaging_disp  = {'Mean coupling','ROI 8','ROI 10','ROI 14','ROI 23'};

    % collect rows
    rows = table;
    for s = 1:numel(scales)
        for i = 1:numel(imaging_order)
            idx = strcmp(string(T.Imaging), imaging_order{i}) & string(T.Scale) == scales(s);
            if any(idx)
                r = find(idx, 1, 'first');
                tmp = table;
                tmp.Scale = scales(s);
                tmp.Imaging = string(imaging_disp{i});
                tmp.Beta = T.Beta_raw(r);
                tmp.Tstat = T.T_raw(r);
                tmp.P = T.P_raw(r);
                tmp.Q = T.Q_raw(r);
                if ismember('SpearmanRho', T.Properties.VariableNames)
                    tmp.Rho = T.SpearmanRho(r);
                else
                    tmp.Rho = nan;
                end
                rows = [rows; tmp]; %#ok<AGROW>
            end
        end
    end

    if isempty(rows)
        error('clinical_association_AD_only.txt 中没有匹配到可画的数据。');
    end

    rows.SE = abs(rows.Beta ./ max(abs(rows.Tstat), 1e-6));
    rows.CI_lo = rows.Beta - 1.96 .* rows.SE;
    rows.CI_hi = rows.Beta + 1.96 .* rows.SE;

    % figure with one panel per selected scale
    nP = numel(scales);
    fig = figure('Color','w', 'Position', [80 80 420*nP 760]);
    tl = tiledlayout(fig, 1, nP, 'TileSpacing','compact', 'Padding','compact');

    colors = [0.18 0.18 0.18;
              0.80 0.35 0.35;
              0.85 0.55 0.35;
              0.35 0.60 0.84;
              0.45 0.72 0.50];

    for s = 1:nP
        ax = nexttile(tl, s); hold(ax, 'on');
        R = rows(rows.Scale == scales(s), :);

        y = 1:height(R);
        for i = 1:height(R)
            c = colors(i,:);
            plot(ax, [R.CI_lo(i) R.CI_hi(i)], [y(i) y(i)], '-', 'Color', c, 'LineWidth', 2);
            scatter(ax, R.Beta(i), y(i), 55, 'o', 'MarkerFaceColor', c, 'MarkerEdgeColor', 'none');

            txt = sprintf('p=%s', local_p(R.P(i)));
            if isfinite(R.Q(i))
                txt = sprintf('p=%s\nq=%s', local_p(R.P(i)), local_p(R.Q(i)));
            end
            text(ax, R.CI_hi(i) + 0.03*max(1, range([min(R.CI_lo) max(R.CI_hi)])), y(i), txt, ...
                'FontName','Arial', 'FontSize',8, 'Color', [0.25 0.25 0.25], ...
                'HorizontalAlignment','left', 'VerticalAlignment','middle');
        end

        xline(ax, 0, '--', 'Color', [0.35 0.35 0.35], 'LineWidth', 1);
        ax.YTick = y;
        ax.YTickLabel = cellstr(R.Imaging);
        ax.YDir = 'reverse';
        ax.Box = 'off';
        ax.TickDir = 'out';
        ax.FontName = 'Arial';
        ax.FontSize = 10;
        xlabel(ax, '\beta (95% CI)', 'FontName','Arial', 'FontSize',11);
        title(ax, char(strrep(scales(s), '_', ' ')), 'FontName','Arial', 'FontSize',12, 'FontWeight','bold');
    end

    sgtitle('Figure 12. Clinical associations in AD-only cohort', 'FontName','Arial', 'FontSize',13, 'FontWeight','bold');

    out_tif = fullfile(out_dir, 'Figure12_clinical_dumbbell.tif');
    out_png = fullfile(out_dir, 'Figure12_clinical_dumbbell.png');
    exportgraphics(fig, out_tif, 'Resolution', 600, 'BackgroundColor','white', 'ContentType','image');
    exportgraphics(fig, out_png, 'Resolution', 300, 'BackgroundColor','white', 'ContentType','image');

    writetable(rows, fullfile(out_dir, 'Figure12_clinical_dumbbell_table.csv'));
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
