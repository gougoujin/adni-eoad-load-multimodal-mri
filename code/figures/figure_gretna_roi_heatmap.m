% Repository usage summary
% Description: Create a heatmap summarizing recurrent GRETNA nodal findings.
% Usage: Run after summarize_gretna_bonferroni_results has produced master/detail tables.
% Outputs: Exports figure images and plotting tables.
% Note: Update input paths, toolboxes, and filenames for your local environment.

% Public repository version: update file paths, toolboxes, and local settings before running.
% This script/function was lightly sanitized for sharing and may require project-specific inputs.

function out_tif = figure_gretna_roi_heatmap(detail_file, master_file, out_dir, top_n_rois)
% Figure 3 (replacement): GRETNA nodal ROI x metric heatmap
% More publication-like than abstract summary plots.
% Default export: TIFF (lossless)
%
% Inputs:
%   detail_file : gretna_bonferroni_detail_summary_autoTF.txt
%   master_file : gretna_bonferroni_master_summary_autoTF.txt
%   out_dir     : output folder
%   top_n_rois  : number of recurrent ROIs to show (default 20)
%
% Example:
%   figure_gretna_roi_heatmap

    if nargin < 1 || isempty(detail_file)
        detail_file = '';
    end
    if nargin < 2 || isempty(master_file)
        master_file = '';
    end
    if nargin < 3 || isempty(out_dir)
        out_dir = '';
    end
    if nargin < 4 || isempty(top_n_rois)
        top_n_rois = 20;
    end

    if ~isfile(detail_file)
        error(['缺少文件: gretna_bonferroni_detail_summary_autoTF.txt\n' ...
               '原路径应为: /path/to/project' ...
               '当前路径: ' detail_file]);
    end
    if ~isfile(master_file)
        error(['缺少文件: gretna_bonferroni_master_summary_autoTF.txt\n' ...
               '原路径应为: /path/to/project' ...
               '当前路径: ' master_file]);
    end
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end

    Td = readtable(detail_file, 'FileType','text', 'Delimiter','\t');
    Tm = readtable(master_file, 'FileType','text', 'Delimiter','\t');

    Td.Properties.VariableNames = matlab.lang.makeValidName(Td.Properties.VariableNames);
    Tm.Properties.VariableNames = matlab.lang.makeValidName(Tm.Properties.VariableNames);

    reqD = {'Metric','Contrast','ROI','Label','P','Stat'};
    reqM = {'Metric','Contrast','N_Significant'};
    for i = 1:numel(reqD)
        if ~ismember(reqD{i}, Td.Properties.VariableNames)
            error('detail表缺少列: %s', reqD{i});
        end
    end
    for i = 1:numel(reqM)
        if ~ismember(reqM{i}, Tm.Properties.VariableNames)
            error('master表缺少列: %s', reqM{i});
        end
    end

    % Keep only AAL116 nodal entries
    roi_num = Td.ROI;
    if iscell(roi_num) || isstring(roi_num)
        roi_num = str2double(string(roi_num));
    end
    keep = isfinite(roi_num) & roi_num >= 1 & roi_num <= 116;
    Td = Td(keep, :);

    if isempty(Td)
        error('detail表中没有可用的 1-116 ROI 结果。');
    end

    % Keep only metric-contrast combinations with at least one significant node
    if ismember('IsNodal116', Tm.Properties.VariableNames)
        keepM = Tm.IsNodal116 == 1 & Tm.N_Significant > 0;
    else
        keepM = Tm.N_Significant > 0;
    end
    Tm = Tm(keepM, :);

    if isempty(Tm)
        error('master表中没有可用的显著 nodal metric/contrast 组合。');
    end

    % Order columns by contrast then metric
    contrast_order = {'group1_vs_group2','group1_vs_group3','group2_vs_group3','group_main_effect'};
    metric_order = {'BetweennessCentrality','DegreeCentrality','NodalClustCoeff', ...
                    'NodalLocalEfficiency','NodalShortestPath','NodalEfficiency'};

    Tm.cOrder = zeros(height(Tm),1);
    Tm.mOrder = zeros(height(Tm),1);

    for i = 1:height(Tm)
        ci = find(strcmp(contrast_order, char(Tm.Contrast{i})), 1);
        mi = find(strcmp(metric_order, char(Tm.Metric{i})), 1);
        if isempty(ci), ci = 99; end
        if isempty(mi), mi = 99; end
        Tm.cOrder(i) = ci;
        Tm.mOrder(i) = mi;
    end

    Tm = sortrows(Tm, {'cOrder','mOrder','N_Significant'}, {'ascend','ascend','descend'});

    col_labels = string(Tm.Contrast) + " | " + string(Tm.Metric);

    % Find most recurrent ROIs across all significant combinations
    roi_labels_all = string(Td.Label);
    [u_roi, ~, ic] = unique(roi_labels_all, 'stable');
    counts = accumarray(ic, 1);
    mean_abs_stat = accumarray(ic, abs(Td.Stat), [], @mean);

    score = counts * 100 + mean_abs_stat; % prioritize recurrence, break ties by magnitude
    [~, ord] = sort(score, 'descend');
    u_roi = u_roi(ord);
    counts = counts(ord);

    top_n_rois = min(top_n_rois, numel(u_roi));
    row_labels = u_roi(1:top_n_rois);

    % Build matrix
    M = nan(numel(row_labels), numel(col_labels));
    P = nan(numel(row_labels), numel(col_labels));

    for r = 1:numel(row_labels)
        for c = 1:numel(col_labels)
            contrast_str = string(Tm.Contrast(c));
            metric_str   = string(Tm.Metric(c));
            idx = string(Td.Label) == row_labels(r) & ...
                  string(Td.Contrast) == contrast_str & ...
                  string(Td.Metric) == metric_str;
            if any(idx)
                stats = Td.Stat(idx);
                pvals = Td.P(idx);
                [~, ix] = max(abs(stats));
                M(r,c) = stats(ix);
                P(r,c) = pvals(ix);
            end
        end
    end

    % Clean display labels
    row_disp = regexprep(row_labels, '^\d+\s+', '');
    row_disp = regexprep(row_disp, '_', ' ');

    contrast_short = strings(height(Tm),1);
    for i = 1:height(Tm)
        contrast_short(i) = local_contrast_short(string(Tm.Contrast(i)));
    end
    metric_short = strings(height(Tm),1);
    for i = 1:height(Tm)
        metric_short(i) = local_metric_short(string(Tm.Metric(i)));
    end
    col_disp = contrast_short + " | " + metric_short;

    % Figure
    fig = figure('Color','w', 'Position', [80 80 1800 980]);
    tl = tiledlayout(fig, 1, 2, 'TileSpacing','compact', 'Padding','compact');

    % Panel A: counts by metric-contrast
    ax1 = nexttile(tl, 1);
    barh(ax1, Tm.N_Significant, 'FaceColor', [0.35 0.45 0.78], 'EdgeColor', 'none');
    ax1.YDir = 'reverse';
    ax1.YTick = 1:numel(col_disp);
    ax1.YTickLabel = cellstr(col_disp);
    ax1.FontName = 'Arial';
    ax1.FontSize = 10;
    ax1.TickDir = 'out';
    ax1.Box = 'off';
    xlabel(ax1, 'Number of significant nodes', 'FontName','Arial','FontSize',11);
    title(ax1, 'A  Significant nodal findings by contrast and metric', ...
        'FontName','Arial','FontSize',13,'FontWeight','bold');

    % separator lines for contrasts
    hold(ax1, 'on');
    for i = 1:height(Tm)-1
        if ~strcmp(string(Tm.Contrast(i)), string(Tm.Contrast(i+1)))
            yline(ax1, i+0.5, '-', 'Color', [0.82 0.82 0.82], 'LineWidth', 1);
        end
    end

    % Panel B: heatmap
    ax2 = nexttile(tl, 2);
    imagesc(ax2, M);
    axis(ax2, 'tight');
    set(ax2, 'YDir','normal');
    ax2.XTick = 1:numel(col_disp);
    ax2.XTickLabel = cellstr(col_disp);
    ax2.YTick = 1:numel(row_disp);
    ax2.YTickLabel = cellstr(row_disp);
    xtickangle(ax2, 45);
    ax2.FontName = 'Arial';
    ax2.FontSize = 10;
    ax2.TickDir = 'out';
    ax2.Box = 'off';
    title(ax2, sprintf('B  Top %d recurrent significant ROIs', numel(row_disp)), ...
        'FontName','Arial','FontSize',13,'FontWeight','bold');

    % Diverging colormap
    clim = max(abs(M(:)), [], 'omitnan');
    if ~isfinite(clim) || clim == 0
        clim = 1;
    end
    caxis(ax2, [-clim clim]);
    colormap(ax2, local_bwr(256));
    cb = colorbar(ax2);
    cb.Label.String = 'Statistic value';
    cb.FontName = 'Arial';
    cb.FontSize = 10;

    % mark cells with p<0.01
    hold(ax2, 'on');
    for r = 1:size(P,1)
        for c = 1:size(P,2)
            if isfinite(P(r,c)) && P(r,c) < 0.01
                text(ax2, c, r, '•', 'HorizontalAlignment','center', ...
                    'VerticalAlignment','middle', 'FontWeight','bold', ...
                    'Color','k', 'FontSize', 12);
            end
        end
    end

    % separator lines for contrasts
    for i = 1:height(Tm)-1
        if ~strcmp(string(Tm.Contrast(i)), string(Tm.Contrast(i+1)))
            xline(ax2, i+0.5, '-', 'Color', [0.2 0.2 0.2], 'LineWidth', 1);
        end
    end

    annotation(fig, 'textbox', [0.12 0.01 0.78 0.04], 'String', ...
        'Color shows signed test statistic. Dots mark cells with p<0.01. Columns are grouped by contrast.', ...
        'EdgeColor','none', 'HorizontalAlignment','center', 'FontName','Arial', 'FontSize',10);

    out_tif = fullfile(out_dir, 'Figure3_gretna_roi_heatmap.tif');
    out_png = fullfile(out_dir, 'Figure3_gretna_roi_heatmap.png');
    exportgraphics(fig, out_tif, 'Resolution', 600, 'BackgroundColor','white', 'ContentType','image');
    exportgraphics(fig, out_png, 'Resolution', 300, 'BackgroundColor','white', 'ContentType','image');

    % save source matrix
    out_csv = fullfile(out_dir, 'Figure3_gretna_roi_heatmap_matrix.csv');
    Tsave = array2table(M, 'VariableNames', matlab.lang.makeValidName(cellstr(col_disp)));
    Tsave = addvars(Tsave, cellstr(row_disp), 'Before', 1, 'NewVariableNames', 'ROI');
    writetable(Tsave, out_csv);

    fprintf('Save/path/to/project', out_tif, out_png, out_csv);
end

function s = local_contrast_short(c)
    if c == "group1_vs_group2"
        s = "EOAD vs LOAD";
    elseif c == "group1_vs_group3"
        s = "EOAD vs HC";
    elseif c == "group2_vs_group3"
        s = "LOAD vs HC";
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

function cmap = local_bwr(n)
    if nargin < 1, n = 256; end
    n1 = floor(n/2);
    n2 = n - n1;
    b = [linspace(0,1,n1)' linspace(0,1,n1)' ones(n1,1)];
    r = [ones(n2,1) linspace(1,0,n2)' linspace(1,0,n2)'];
    cmap = [b; r];
end
