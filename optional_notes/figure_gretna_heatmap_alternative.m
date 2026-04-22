% Repository usage summary
% Description: Alternative GRETNA heatmap implementation kept for reference.
% Usage: Optional figure variant.
% Outputs: Not required for the main manuscript figure set.
% Note: Update input paths, toolboxes, and filenames for your local environment.

% Public repository version: update file paths, toolboxes, and local settings before running.
% This script/function was lightly sanitized for sharing and may require project-specific inputs.

function out_tif = figure_gretna_heatmap_alternative(detail_file, master_file, out_dir)
% Make Figure 2 from actual GRETNA summary table structure.
% Default export: TIFF (lossless)

    if nargin < 1 || isempty(detail_file)
        detail_file = '';
    end
    if nargin < 2 || isempty(master_file)
        master_file = '';
    end
    if nargin < 3 || isempty(out_dir)
        out_dir = '';
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

    % Actual columns in your files:
    % detail: Folder StatType Metric Contrast ROI Label P Stat ...
    % master: Folder StatType Metric Contrast ... N_Significant ...
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

    % Keep nodal rows only
    roi_num = Td.ROI;
    if iscell(roi_num) || isstring(roi_num)
        roi_num = str2double(string(roi_num));
    end
    keep = isfinite(roi_num) & roi_num >= 1 & roi_num <= 116;
    Td = Td(keep, :);

    if isempty(Td)
        error('detail表中没有可用的 1-116 ROI 结果。');
    end

    Td.metric_contrast = string(Td.Metric) + " | " + string(Td.Contrast);
    Tm.metric_contrast = string(Tm.Metric) + " | " + string(Tm.Contrast);

    % Order columns by N_Significant descending
    [~, ordM] = sort(Tm.N_Significant, 'descend', 'MissingPlacement','last');
    Tm = Tm(ordM, :);
    col_labels = unique(string(Tm.metric_contrast), 'stable');

    % Order rows by frequency then mean abs(stat)
    row_labels = unique(string(Td.Label), 'stable');
    freq = zeros(numel(row_labels),1);
    mag  = zeros(numel(row_labels),1);
    for i = 1:numel(row_labels)
        idx = string(Td.Label) == row_labels(i);
        freq(i) = sum(idx);
        mag(i) = mean(abs(Td.Stat(idx)), 'omitnan');
    end
    [~, ordR] = sortrows([freq mag], [-1 -2]);
    row_labels = row_labels(ordR);

    % Matrix
    M = nan(numel(row_labels), numel(col_labels));
    P = nan(numel(row_labels), numel(col_labels));

    for r = 1:numel(row_labels)
        for c = 1:numel(col_labels)
            idx = string(Td.Label) == row_labels(r) & string(Td.metric_contrast) == col_labels(c);
            if any(idx)
                stats = Td.Stat(idx);
                pvals = Td.P(idx);
                [~, ix] = max(abs(stats));
                M(r,c) = stats(ix);
                P(r,c) = pvals(ix);
            end
        end
    end

    keep_r = any(isfinite(M), 2);
    keep_c = any(isfinite(M), 1);
    M = M(keep_r, keep_c);
    P = P(keep_r, keep_c);
    row_labels = row_labels(keep_r);
    col_labels = col_labels(keep_c);

    % Count values from master
    count_vals = zeros(1, numel(col_labels));
    for c = 1:numel(col_labels)
        idx = string(Tm.metric_contrast) == col_labels(c);
        if any(idx)
            count_vals(c) = max(Tm.N_Significant(idx));
        else
            count_vals(c) = sum(isfinite(M(:,c)));
        end
    end

    % Shorter display labels
    row_disp = regexprep(row_labels, '_', ' ');
    col_disp = regexprep(col_labels, '_', ' ');

    fig = figure('Color','w', 'Position', [80 80 1900 950]);
    tl = tiledlayout(fig, 1, 2, 'Padding','compact', 'TileSpacing','compact');

    ax1 = nexttile(tl, 1);
    barh(ax1, count_vals, 'FaceColor', [0.45 0.45 0.75], 'EdgeColor', 'none');
    ax1.YDir = 'reverse';
    ax1.YTick = 1:numel(col_disp);
    ax1.YTickLabel = cellstr(col_disp);
    ax1.FontName = 'Arial';
    ax1.FontSize = 10;
    ax1.Box = 'off';
    ax1.TickDir = 'out';
    xlabel(ax1, 'Number of significant nodes', 'FontName','Arial', 'FontSize',11);
    title(ax1, 'A  Significant node counts by metric/contrast', 'FontName','Arial', 'FontSize',12, 'FontWeight','bold');

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
    ax2.FontSize = 9;
    ax2.TickDir = 'out';
    ax2.Box = 'off';
    title(ax2, 'B  Signed nodal statistics heatmap', 'FontName','Arial', 'FontSize',12, 'FontWeight','bold');
    colormap(ax2, parula(256));
    cb = colorbar(ax2);
    cb.Label.String = 'Statistic value';
    cb.FontName = 'Arial';
    cb.FontSize = 10;

    hold(ax2, 'on');
    for r = 1:size(P,1)
        for c = 1:size(P,2)
            if isfinite(P(r,c))
                if P(r,c) < 0.001
                    txt = '**';
                elseif P(r,c) < 0.01
                    txt = '*';
                else
                    txt = '';
                end
                if ~isempty(txt)
                    text(ax2, c, r, txt, 'HorizontalAlignment','center', ...
                        'VerticalAlignment','middle', 'FontWeight','bold', ...
                        'Color','w', 'FontSize',9);
                end
            end
        end
    end

    out_tif = fullfile(out_dir, 'Figure2_gretna_nodal_summary.tif');
    out_png = fullfile(out_dir, 'Figure2_gretna_nodal_summary.png');
    exportgraphics(fig, out_tif, 'Resolution', 600, 'BackgroundColor','white', 'ContentType','image');
    exportgraphics(fig, out_png, 'Resolution', 300, 'BackgroundColor','white', 'ContentType','image');

    out_csv = fullfile(out_dir, 'Figure2_gretna_nodal_summary_matrix.csv');
    Tsave = array2table(M, 'VariableNames', matlab.lang.makeValidName(cellstr(col_labels)));
    Tsave = addvars(Tsave, cellstr(row_labels), 'Before', 1, 'NewVariableNames', 'ROI_Label');
    writetable(Tsave, out_csv);

    fprintf('Save/path/to/project', out_tif, out_png, out_csv);
end
