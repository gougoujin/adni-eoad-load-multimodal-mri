% Repository usage summary
% Description: Create a whole-brain ROI rank plot for structural-functional coupling effects.
% Usage: Run after the ROI-level coupling statistics table is available.
% Outputs: Exports ranked effect plots and supporting CSV tables.
% Note: Update input paths, toolboxes, and filenames for your local environment.

% Public repository version: update file paths, toolboxes, and local settings before running.
% This script/function was lightly sanitized for sharing and may require project-specific inputs.

function out_tif = figure_roi_coupling_rankplot(stats_file, label_file, out_dir)
% Figure 9: whole-brain ROI-level coupling rank plot
% Highlights the 4 frontal ROIs used in Figure 1.
%
% Default export: TIFF (lossless)
%
% Inputs:
%   stats_file = /path/to/project
%   label_file = /path/to/project
%   out_dir    = /path/to/project
%
% Outputs:
%   Figure9_roi_coupling_rankplot.tif
%   Figure9_roi_coupling_rankplot.png
%   Figure9_roi_coupling_rankplot_table.csv

    if nargin < 1 || isempty(stats_file)
        stats_file = '';
    end
    if nargin < 2 || isempty(label_file)
        label_file = '';
    end
    if nargin < 3 || isempty(out_dir)
        out_dir = '';
    end

    if ~isfile(stats_file)
        error(['缺少文件: roi_level_coupling_stats_age_sex_fdr.txt\n' ...
               '原路径应为: /path/to/project' ...
               '当前路径: ' stats_file]);
    end
    if ~isfile(label_file)
        error(['缺少文件: AAL116_Labels.txt\n' ...
               '原路径应为: /path/to/project' ...
               '当前路径: ' label_file]);
    end
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end

    T = readtable(stats_file, 'FileType','text', 'Delimiter','\t');
    T.Properties.VariableNames = matlab.lang.makeValidName(T.Properties.VariableNames);

    req = {'ROI','beta_group2','p_group2','q_group2'};
    for i = 1:numel(req)
        if ~ismember(req{i}, T.Properties.VariableNames)
            error('统计表缺少列: %s', req{i});
        end
    end

    if ~ismember('fdr_pass_group2', T.Properties.VariableNames)
        T.fdr_pass_group2 = double(T.q_group2 < 0.05);
    end

    labels = local_read_labels(label_file, 116);

    roi = T.ROI;
    beta = T.beta_group2;
    p = T.p_group2;
    q = T.q_group2;
    pass = T.fdr_pass_group2;

    ok = isfinite(roi) & isfinite(beta) & isfinite(p) & isfinite(q);
    roi = roi(ok);
    beta = beta(ok);
    p = p(ok);
    q = q(ok);
    pass = pass(ok);

    [beta_sorted, ord] = sort(beta, 'descend');
    roi_sorted = roi(ord);
    p_sorted = p(ord);
    q_sorted = q(ord);
    pass_sorted = pass(ord);

    lbl = strings(numel(roi_sorted),1);
    for i = 1:numel(roi_sorted)
        rid = roi_sorted(i);
        if rid >= 1 && rid <= numel(labels)
            lbl(i) = string(labels{rid});
        else
            lbl(i) = "ROI_" + string(rid);
        end
    end

    % Figure 1 frontal set in your current pipeline
    frontal_ids = [3, 7, 15, 23];
    frontal_names = ["Frontal_Sup_Medial_L","Frontal_Mid_R","Frontal_Mid_Orb_R","Frontal_Inf_Tri_R"];
    is_frontal = ismember(roi_sorted, frontal_ids);

    x = 1:numel(beta_sorted);
    yq = -log10(q_sorted);

    fig = figure('Color','w', 'Position', [80 80 1550 820]);
    tl = tiledlayout(fig, 2, 1, 'TileSpacing','compact', 'Padding','compact');

    % Panel A
    ax1 = nexttile(tl, 1); hold(ax1, 'on');
    plot(ax1, x, beta_sorted, '-', 'Color', [0.78 0.78 0.78], 'LineWidth', 1);
    scatter(ax1, x, beta_sorted, 20, [0.62 0.62 0.62], 'filled', 'MarkerEdgeColor', 'none');

    sig_idx = find(pass_sorted == 1);
    if ~isempty(sig_idx)
        scatter(ax1, sig_idx, beta_sorted(sig_idx), 34, [0.88 0.33 0.33], ...
            'filled', 'MarkerEdgeColor', 'none');
    end

    idx4 = find(is_frontal);
    if ~isempty(idx4)
        scatter(ax1, idx4, beta_sorted(idx4), 72, 'o', ...
            'MarkerFaceColor', [0.15 0.15 0.15], ...
            'MarkerEdgeColor', [1 1 1], 'LineWidth', 0.8);
        for k = 1:numel(idx4)
            ii = idx4(k);
            nm = local_short_label(char(lbl(ii)));
            text(ax1, ii, beta_sorted(ii), ['  ' nm], ...
                'FontName','Arial', 'FontSize',10, 'Color', [0.1 0.1 0.1], ...
                'HorizontalAlignment','left', 'VerticalAlignment','middle');
        end
    end

    yline(ax1, 0, '--', 'Color', [0.35 0.35 0.35], 'LineWidth', 1);
    ax1.Box = 'off';
    ax1.TickDir = 'out';
    ax1.FontName = 'Arial';
    ax1.FontSize = 11;
    xlim(ax1, [0 numel(beta_sorted)+1]);
    ylabel(ax1, '\beta (LOAD vs EOAD)', 'FontName','Arial', 'FontSize',11);
    title(ax1, 'A  Whole-brain ROI-level coupling effects ranked by adjusted \beta', ...
        'FontName','Arial', 'FontSize',12, 'FontWeight','bold');

    % Panel B
    ax2 = nexttile(tl, 2); hold(ax2, 'on');
    bh = bar(ax2, x, yq, 1.0, 'FaceColor', [0.75 0.81 0.90], 'EdgeColor', 'none');
    if ~isempty(sig_idx)
        bar(ax2, sig_idx, yq(sig_idx), 1.0, 'FaceColor', [0.88 0.33 0.33], 'EdgeColor', 'none');
    end
    if ~isempty(idx4)
        scatter(ax2, idx4, yq(idx4)+0.03, 42, 'kv', 'filled');
    end
    yline(ax2, -log10(0.05), '--', 'Color', [0.35 0.35 0.35], 'LineWidth', 1);

    ax2.Box = 'off';
    ax2.TickDir = 'out';
    ax2.FontName = 'Arial';
    ax2.FontSize = 11;
    xlim(ax2, [0 numel(beta_sorted)+1]);
    xlabel(ax2, 'ROI rank', 'FontName','Arial', 'FontSize',11);
    ylabel(ax2, '-log_{10}(q)', 'FontName','Arial', 'FontSize',11);
    title(ax2, 'B  FDR evidence strength for LOAD vs EOAD contrast', ...
        'FontName','Arial', 'FontSize',12, 'FontWeight','bold');

    sgtitle('Whole-brain regional structural-functional coupling profile', ...
        'FontName','Arial', 'FontSize',13, 'FontWeight','bold');

    out_tif = fullfile(out_dir, 'Figure9_roi_coupling_rankplot.tif');
    out_png = fullfile(out_dir, 'Figure9_roi_coupling_rankplot.png');
    exportgraphics(fig, out_tif, 'Resolution', 600, 'BackgroundColor','white', 'ContentType','image');
    exportgraphics(fig, out_png, 'Resolution', 300, 'BackgroundColor','white', 'ContentType','image');

    Tout = table((1:numel(roi_sorted))', roi_sorted, lbl, beta_sorted, p_sorted, q_sorted, pass_sorted, ...
        'VariableNames', {'Rank','ROI','Label','Beta_group2','P_group2','Q_group2','FDR_pass'});
    writetable(Tout, fullfile(out_dir, 'Figure9_roi_coupling_rankplot_table.csv'));
end

function labels = local_read_labels(label_file, nROI)
    labels = cell(nROI,1);
    fid = fopen(label_file, 'r');
    if fid < 0
        error('无法打开标签文件: %s', label_file);
    end
    while ~feof(fid)
        line = strtrim(fgetl(fid));
        if isempty(line), continue; end
        parts = regexp(line, '\s+', 'split');
        if numel(parts) < 2, continue; end
        idx = str2double(parts{1});
        if isfinite(idx) && idx >= 1 && idx <= nROI
            labels{idx} = strtrim(parts{2});
        end
    end
    fclose(fid);
    for i = 1:nROI
        if isempty(labels{i})
            labels{i} = sprintf('ROI_%d', i);
        end
    end
end

function s = local_short_label(name_in)
    name_in = strrep(name_in, 'Frontal_Sup_Medial_L', 'SFGmed-L');
    name_in = strrep(name_in, 'Frontal_Mid_R', 'MFG-R');
    name_in = strrep(name_in, 'Frontal_Mid_Orb_R', 'ORBmid-R');
    name_in = strrep(name_in, 'Frontal_Inf_Tri_R', 'IFGtri-R');
    s = name_in;
end
