% Repository usage summary
% Description: Create a forest plot summarizing representative SDI, FA, and MD effects.
% Usage: Run after SDI and diffusion ROI summary files are available.
% Outputs: Exports figure images and a combined plotting table.
% Note: Update input paths, toolboxes, and filenames for your local environment.

% Public repository version: update file paths, toolboxes, and local settings before running.
% This script/function was lightly sanitized for sharing and may require project-specific inputs.

function out_tif = figure_sdi_fa_md_forest(sdi_file, fa_file, md_file, out_dir, top_n)
% Figure 11: SDI / FA / MD forest plot
% v2 uses the real FA/MD column names:
%   Beta_g2_vs_g1, T_g2_vs_g1, P_g2_vs_g1
%
% Default:
%   sdi_file = /path/to/project
%   fa_file  = /path/to/project
%   md_file  = /path/to/project
%   out_dir  = /path/to/project

    if nargin < 1 || isempty(sdi_file), sdi_file = ''; end
    if nargin < 2 || isempty(fa_file),  fa_file  = ''; end
    if nargin < 3 || isempty(md_file),  md_file  = ''; end
    if nargin < 4 || isempty(out_dir),  out_dir  = ''; end
    if nargin < 5 || isempty(top_n),    top_n    = 6; end

    if ~isfile(sdi_file)
        error(['缺少文件: roi_level_sdi_stats_age_sex.txt\n原路径应为: /path/to/project ' sdi_file]);
    end
    if ~isfile(fa_file)
        error(['缺少文件: roi_level_fa_stats_age_sex_fixed.txt\n原路径应为: /path/to/project ' fa_file]);
    end
    if ~isfile(md_file)
        error(['缺少文件: roi_level_md_stats_age_sex_fixed.txt\n原路径应为: /path/to/project ' md_file]);
    end
    if ~exist(out_dir, 'dir'), mkdir(out_dir); end

    Ts = readtable(sdi_file, 'FileType','text', 'Delimiter','\t');
    Tf = readtable(fa_file,  'FileType','text', 'Delimiter','\t');
    Tm = readtable(md_file,  'FileType','text', 'Delimiter','\t');

    Ts.Properties.VariableNames = matlab.lang.makeValidName(Ts.Properties.VariableNames);
    Tf.Properties.VariableNames = matlab.lang.makeValidName(Tf.Properties.VariableNames);
    Tm.Properties.VariableNames = matlab.lang.makeValidName(Tm.Properties.VariableNames);

    Rs = local_pick_top_sdi(Ts, top_n);
    Rf = local_pick_top_famd(Tf, "FA", top_n);
    Rm = local_pick_top_famd(Tm, "MD", top_n);
    R = [Rs; Rf; Rm];

    fig = figure('Color','w', 'Position', [80 80 1550 900]);
    ax = axes(fig); hold(ax,'on');
    y = 1:height(R);

    colors = zeros(height(R),3);
    for i = 1:height(R)
        switch char(R.Modality(i))
            case 'SDI'
                colors(i,:) = [0.48 0.62 0.88];
            case 'FA'
                colors(i,:) = [0.86 0.52 0.40];
            case 'MD'
                colors(i,:) = [0.45 0.72 0.50];
        end
    end

    xr = [min(R.CI_lo, [], 'omitnan') max(R.CI_hi, [], 'omitnan')];
    if any(~isfinite(xr)) || xr(1) == xr(2), xr = [-1 1]; end
    xoff = 0.02 * range(xr);

    for i = 1:height(R)
        if isfinite(R.CI_lo(i)) && isfinite(R.CI_hi(i))
            plot(ax, [R.CI_lo(i) R.CI_hi(i)], [y(i) y(i)], '-', 'Color', colors(i,:), 'LineWidth', 2);
        end
        scatter(ax, R.Beta(i), y(i), 50, 'o', 'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', 'none');
        text(ax, R.Beta(i) + xoff, y(i), sprintf('p=%s', local_p(R.P(i))), ...
            'FontName','Arial', 'FontSize',9, 'Color',[0.25 0.25 0.25], ...
            'HorizontalAlignment','left', 'VerticalAlignment','middle');
    end

    xline(ax, 0, '--', 'Color', [0.35 0.35 0.35], 'LineWidth', 1);

    idx_s = find(R.Modality=="SDI", 1, 'last');
    idx_f = find(R.Modality=="FA", 1, 'last');
    if ~isempty(idx_s) && idx_s < height(R), yline(ax, idx_s+0.5, '-', 'Color', [0.72 0.72 0.72], 'LineWidth', 1); end
    if ~isempty(idx_f) && idx_f < height(R), yline(ax, idx_f+0.5, '-', 'Color', [0.72 0.72 0.72], 'LineWidth', 1); end

    ax.YTick = y;
    ax.YTickLabel = cellstr(R.Label);
    ax.YDir = 'reverse';
    ax.Box = 'off';
    ax.TickDir = 'out';
    ax.FontName = 'Arial';
    ax.FontSize = 10;
    xlabel(ax, 'Adjusted \beta for LOAD vs EOAD (95% CI when available)', 'FontName','Arial', 'FontSize',11);
    title(ax, 'Figure 11. SDI / FA / MD supplementary effect summary', 'FontName','Arial', 'FontSize',13, 'FontWeight','bold');

    out_tif = fullfile(out_dir, 'Figure11_sdi_fa_md_forest.tif');
    out_png = fullfile(out_dir, 'Figure11_sdi_fa_md_forest.png');
    exportgraphics(fig, out_tif, 'Resolution', 600, 'BackgroundColor','white', 'ContentType','image');
    exportgraphics(fig, out_png, 'Resolution', 300, 'BackgroundColor','white', 'ContentType','image');

    writetable(R, fullfile(out_dir, 'Figure11_sdi_fa_md_forest_table.csv'));
end

function R = local_pick_top_sdi(T, top_n)
req = {'ROI','beta_group2','p_group2'};
for i = 1:numel(req)
    if ~ismember(req{i}, T.Properties.VariableNames), error('SDI表缺少列: %s', req{i}); end
end
if ismember('t_group2', T.Properties.VariableNames), tval = T.t_group2; else, tval = nan(height(T),1); end

T.absbeta = abs(T.beta_group2);
T = sortrows(T, {'absbeta','p_group2'}, {'descend','ascend'});
T = T(1:min(top_n, height(T)), :);

R = table;
R.Modality = repmat("SDI", height(T), 1);
R.Label = "SDI | ROI " + string(T.ROI);
R.Beta = T.beta_group2;
R.P = T.p_group2;
R.Tstat = tval(1:min(top_n, numel(tval)));
R.SE = abs(R.Beta ./ max(abs(R.Tstat), 1e-6));
R.CI_lo = R.Beta - 1.96 .* R.SE;
R.CI_hi = R.Beta + 1.96 .* R.SE;
if ~ismember('t_group2', T.Properties.VariableNames), R.SE(:)=nan; R.CI_lo(:)=nan; R.CI_hi(:)=nan; end
end

function R = local_pick_top_famd(T, modality, top_n)
req = {'ROI','Beta_g2_vs_g1','P_g2_vs_g1'};
for i = 1:numel(req)
    if ~ismember(req{i}, T.Properties.VariableNames), error('%s表缺少列: %s', modality, req{i}); end
end
if ismember('T_g2_vs_g1', T.Properties.VariableNames), tval = T.T_g2_vs_g1; else, tval = nan(height(T),1); end

T.absbeta = abs(T.Beta_g2_vs_g1);
T = sortrows(T, {'absbeta','P_g2_vs_g1'}, {'descend','ascend'});
T = T(1:min(top_n, height(T)), :);

R = table;
R.Modality = repmat(string(modality), height(T), 1);
R.Label = string(modality) + " | ROI " + string(T.ROI);
R.Beta = T.Beta_g2_vs_g1;
R.P = T.P_g2_vs_g1;
R.Tstat = tval(1:min(top_n, numel(tval)));
R.SE = abs(R.Beta ./ max(abs(R.Tstat), 1e-6));
R.CI_lo = R.Beta - 1.96 .* R.SE;
R.CI_hi = R.Beta + 1.96 .* R.SE;
if ~ismember('T_g2_vs_g1', T.Properties.VariableNames), R.SE(:)=nan; R.CI_lo(:)=nan; R.CI_hi(:)=nan; end
end

function s = local_p(p)
if ~isfinite(p), s='NA'; elseif p<0.001, s='<0.001'; else, s=sprintf('%.3f', p); end
end
