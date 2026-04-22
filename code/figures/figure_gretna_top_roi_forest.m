% Repository usage summary
% Description: Create a forest/strip plot for the top recurrent GRETNA ROIs.
% Usage: Run after GRETNA detail summary tables are available.
% Outputs: Exports publication-style ROI summary figures.
% Note: Update input paths, toolboxes, and filenames for your local environment.

% Public repository version: update file paths, toolboxes, and local settings before running.
% This script/function was lightly sanitized for sharing and may require project-specific inputs.

function out_tif = figure_gretna_top_roi_forest(detail_file, out_dir, top_n)
% Figure 8: Top recurrent GRETNA ROI forest/strip plot
% Publication-style summary without uncertain spatial projection.
%
% Inputs:
%   detail_file = /path/to/project
%   out_dir     = /path/to/project
%   top_n       = number of recurrent ROIs to show (default 12)
%
% Outputs:
%   Figure8_gretna_toproi_forest.tif
%   Figure8_gretna_toproi_forest.png
%   Figure8_gretna_toproi_forest_table.csv

    if nargin < 1 || isempty(detail_file)
        detail_file = '';
    end
    if nargin < 2 || isempty(out_dir)
        out_dir = '';
    end
    if nargin < 3 || isempty(top_n)
        top_n = 12;
    end

    if ~isfile(detail_file)
        error(['缺少文件: gretna_bonferroni_detail_summary_autoTF.txt\n' ...
               '原路径应为: /path/to/project' ...
               '当前路径: ' detail_file]);
    end
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end

    T = readtable(detail_file, 'FileType','text', 'Delimiter','\t');
    T.Properties.VariableNames = matlab.lang.makeValidName(T.Properties.VariableNames);

    req = {'Metric','Contrast','ROI','Label','P','Stat'};
    for i = 1:numel(req)
        if ~ismember(req{i}, T.Properties.VariableNames)
            error('detail表缺少列: %s', req{i});
        end
    end

    roi_num = T.ROI;
    if iscell(roi_num) || isstring(roi_num)
        roi_num = str2double(string(roi_num));
    end
    keep = isfinite(roi_num) & roi_num >= 1 & roi_num <= 116;
    T = T(keep,:);

    if isempty(T)
        error('detail表中没有可用的 1-116 ROI 结果。');
    end

    % recurrent ROI ranking
    [u,~,ic] = unique(string(T.Label), 'stable');
    freq = accumarray(ic, 1);
    mag  = accumarray(ic, abs(T.Stat), [], @mean);
    score = freq * 100 + mag;
    [~, ord] = sort(score, 'descend');
    u = u(ord); freq = freq(ord); mag = mag(ord);

    top_n = min(top_n, numel(u));
    top_roi = u(1:top_n);

    % keep top roi rows
    TT = T(ismember(string(T.Label), top_roi), :);

    % order rows by recurrence
    row_order = containers.Map(cellstr(top_roi), num2cell(1:numel(top_roi)));

    % compact display
    TT.RowIdx = zeros(height(TT),1);
    for i = 1:height(TT)
        TT.RowIdx(i) = row_order(char(string(TT.Label(i))));
    end
    TT.ContrastDisp = strings(height(TT),1);
    TT.MetricDisp = strings(height(TT),1);
    for i = 1:height(TT)
        TT.ContrastDisp(i) = local_contrast_short(string(TT.Contrast(i)));
        TT.MetricDisp(i) = local_metric_short(string(TT.Metric(i)));
    end
    TT.LabelDisp = regexprep(string(TT.Label), '_', ' ');
    TT.NegLogP = -log10(TT.P);

    % figure
    fig = figure('Color','w', 'Position', [80 80 1500 860]);
    tl = tiledlayout(fig, 1, 2, 'Padding','compact', 'TileSpacing','compact');

    % left panel: recurrence + mean|stat|
    ax1 = nexttile(tl, 1); hold(ax1, 'on');
    yy = 1:top_n;
    barh(ax1, yy, freq(1:top_n), 'FaceColor', [0.72 0.78 0.88], 'EdgeColor', 'none');
    plot(ax1, mag(1:top_n), yy, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5);
    ax1.YDir = 'reverse';
    ax1.YTick = yy;
    ax1.YTickLabel = cellstr(regexprep(top_roi, '_', ' '));
    ax1.FontName = 'Arial';
    ax1.FontSize = 10;
    ax1.Box = 'off';
    ax1.TickDir = 'out';
    xlabel(ax1, 'Recurrence count  (black dot = mean |stat|)', 'FontName','Arial', 'FontSize',11);
    title(ax1, 'A  Most recurrent significant nodal ROIs', 'FontName','Arial', 'FontSize',12, 'FontWeight','bold');

    % right panel: strip / forest style
    ax2 = nexttile(tl, 2); hold(ax2, 'on');
    cmap = lines(6);
    contrast_list = unique(TT.ContrastDisp, 'stable');

    for i = 1:height(TT)
        y = TT.RowIdx(i);
        x = TT.Stat(i);
        psize = 20 + 18 * min(3, TT.NegLogP(i));
        c = local_pick_color(TT.MetricDisp(i), cmap);

        scatter(ax2, x, y, psize, 'o', ...
            'MarkerFaceColor', c, 'MarkerEdgeColor', 'none', ...
            'MarkerFaceAlpha', 0.85);

        % tiny label for contrast
        text(ax2, x, y+0.12, char(TT.ContrastDisp(i)), ...
            'FontName','Arial', 'FontSize',7, 'Color', [0.25 0.25 0.25], ...
            'HorizontalAlignment','center');
    end

    xline(ax2, 0, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1);
    ax2.YDir = 'reverse';
    ax2.YTick = yy;
    ax2.YTickLabel = cellstr(regexprep(top_roi, '_', ' '));
    ax2.FontName = 'Arial';
    ax2.FontSize = 10;
    ax2.Box = 'off';
    ax2.TickDir = 'out';
    xlabel(ax2, 'Signed statistic', 'FontName','Arial', 'FontSize',11);
    title(ax2, 'B  Metric-specific nodal effects', 'FontName','Arial', 'FontSize',12, 'FontWeight','bold');

    legend(ax2, local_dummy_handles(cmap), {'BC','DC','CC','LE','NSP','NE'}, ...
        'Location','bestoutside', 'Box','off', 'FontName','Arial', 'FontSize',9);

    sgtitle('Top recurrent graph-theory nodal findings', 'FontName','Arial', 'FontSize',13, 'FontWeight','bold');

    out_tif = fullfile(out_dir, 'Figure8_gretna_toproi_forest.tif');
    out_png = fullfile(out_dir, 'Figure8_gretna_toproi_forest.png');
    exportgraphics(fig, out_tif, 'Resolution', 600, 'BackgroundColor','white', 'ContentType','image');
    exportgraphics(fig, out_png, 'Resolution', 300, 'BackgroundColor','white', 'ContentType','image');

    % save plot table
    Tsave = sortrows(TT(:, {'Label','Metric','Contrast','P','Stat','NegLogP'}), {'Label','Metric','Contrast'});
    writetable(Tsave, fullfile(out_dir, 'Figure8_gretna_toproi_forest_table.csv'));
end

function s = local_contrast_short(c)
    if c == "group1_vs_group2"
        s = "EvsL";
    elseif c == "group1_vs_group3"
        s = "EvsH";
    elseif c == "group2_vs_group3"
        s = "LvsH";
    elseif c == "group_main_effect"
        s = "Main";
    else
        s = c;
    end
end

function s = local_metric_short(m)
    if m == "BetweennessCentrality"
        s = "BC";
    elseif m == "DegreeCentrality"
        s = "DC";
    elseif m == "NodalClustCoeff"
        s = "CC";
    elseif m == "NodalLocalEfficiency"
        s = "LE";
    elseif m == "NodalShortestPath"
        s = "NSP";
    elseif m == "NodalEfficiency"
        s = "NE";
    else
        s = m;
    end
end

function c = local_pick_color(metric_short, cmap)
    switch char(metric_short)
        case 'BC'
            c = cmap(1,:);
        case 'DC'
            c = cmap(2,:);
        case 'CC'
            c = cmap(3,:);
        case 'LE'
            c = cmap(4,:);
        case 'NSP'
            c = cmap(5,:);
        case 'NE'
            c = cmap(6,:);
        otherwise
            c = [0.4 0.4 0.4];
    end
end

function h = local_dummy_handles(cmap)
    h = gobjects(6,1);
    for i = 1:6
        h(i) = scatter(nan, nan, 50, 'o', 'MarkerFaceColor', cmap(i,:), 'MarkerEdgeColor', 'none');
    end
end
