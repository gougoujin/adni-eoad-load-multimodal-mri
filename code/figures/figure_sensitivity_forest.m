% Repository usage summary
% Description: Create a forest/dumbbell plot for coupling sensitivity analyses.
% Usage: Run after the sensitivity summary table has been generated.
% Outputs: Exports TIFF/PNG figure files and optional plot tables.
% Note: Update input paths, toolboxes, and filenames for your local environment.

% Public repository version: update file paths, toolboxes, and local settings before running.
% This script/function was lightly sanitized for sharing and may require project-specific inputs.

function out_tif = figure_sensitivity_forest(sens_file, out_dir)
% Figure 4: sensitivity analysis forest/dumbbell plot
% Default export: TIFF (lossless)
%
% Input expected from:
%   run_coupling_sensitivity_and_clinical_models.txt
% Output file expected:
%   /path/to/project
%
% Example:
%   figure_sensitivity_forest

    if nargin < 1 || isempty(sens_file)
        sens_file = '';
    end
    if nargin < 2 || isempty(out_dir)
        out_dir = '';
    end

    if ~isfile(sens_file)
        error(['缺少文件: coupling_sensitivity_EOAD_vs_LOAD.txt\n' ...
               '原路径应为: /path/to/project' ...
               '当前路径: ' sens_file '\n' ...
               '请先运行: run_coupling_sensitivity_and_clinical_models']);
    end
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end

    T = readtable(sens_file, 'FileType','text', 'Delimiter','\t');
    T.Properties.VariableNames = matlab.lang.makeValidName(T.Properties.VariableNames);

    req = {'Metric','Beta_LOAD_vs_EOAD_raw','T_raw','P_raw','Beta_LOAD_vs_EOAD_clean','T_clean','P_clean'};
    for i = 1:numel(req)
        if ~ismember(req{i}, T.Properties.VariableNames)
            error('敏感性结果表缺少列: %s', req{i});
        end
    end

    metrics = string(T.Metric);
    metrics(metrics=="mean_coupling") = "Mean coupling";
    metrics(metrics=="ROI8")  = "ROI 8";
    metrics(metrics=="ROI10") = "ROI 10";
    metrics(metrics=="ROI14") = "ROI 14";
    metrics(metrics=="ROI23") = "ROI 23";

    b1 = T.Beta_LOAD_vs_EOAD_raw;
    t1 = T.T_raw;
    p1 = T.P_raw;
    b2 = T.Beta_LOAD_vs_EOAD_clean;
    t2 = T.T_clean;
    p2 = T.P_clean;

    % approximate SE and 95% CI from beta/t
    se1 = abs(b1 ./ t1);
    se2 = abs(b2 ./ t2);
    ci1_lo = b1 - 1.96 .* se1;
    ci1_hi = b1 + 1.96 .* se1;
    ci2_lo = b2 - 1.96 .* se2;
    ci2_hi = b2 + 1.96 .* se2;

    y = 1:height(T);

    fig = figure('Color','w', 'Position', [100 100 1150 520]);

    % left panel: dumbbell
    ax1 = subplot(1,2,1); hold(ax1, 'on');
    for i = 1:height(T)
        plot(ax1, [b1(i) b2(i)], [y(i) y(i)], '-', 'Color', [0.72 0.72 0.72], 'LineWidth', 2);
        plot(ax1, [ci1_lo(i) ci1_hi(i)], [y(i) y(i)], '-', 'Color', [0.85 0.42 0.42], 'LineWidth', 2);
        plot(ax1, [ci2_lo(i) ci2_hi(i)], [y(i)+0.12 y(i)+0.12], '-', 'Color', [0.30 0.56 0.84], 'LineWidth', 2);
        scatter(ax1, b1(i), y(i), 52, 'o', 'MarkerFaceColor', [0.85 0.42 0.42], 'MarkerEdgeColor', 'none');
        scatter(ax1, b2(i), y(i)+0.12, 52, 'o', 'MarkerFaceColor', [0.30 0.56 0.84], 'MarkerEdgeColor', 'none');
    end
    xline(ax1, 0, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1);
    ax1.YTick = y + 0.06;
    ax1.YTickLabel = cellstr(metrics);
    ax1.YDir = 'reverse';
    ax1.Box = 'off';
    ax1.TickDir = 'out';
    ax1.FontName = 'Arial';
    ax1.FontSize = 11;
    xlabel(ax1, 'Beta (LOAD vs EOAD)', 'FontName','Arial', 'FontSize',11);
    title(ax1, 'A  Effect stability before and after outlier removal', ...
        'FontName','Arial', 'FontSize',12, 'FontWeight','bold');
    legend(ax1, {'raw-clean link','raw 95% CI','clean 95% CI','raw','clean'}, ...
        'Location','best', 'Box','off', 'FontName','Arial', 'FontSize',9);

    % right panel: p-value comparison
    ax2 = subplot(1,2,2); hold(ax2, 'on');
    logp1 = -log10(p1);
    logp2 = -log10(p2);
    bw = 0.36;
    barh(ax2, y-bw/2, logp1, bw, 'FaceColor', [0.85 0.42 0.42], 'EdgeColor', 'none');
    barh(ax2, y+bw/2, logp2, bw, 'FaceColor', [0.30 0.56 0.84], 'EdgeColor', 'none');
    xline(ax2, -log10(0.05), '--', 'Color', [0.35 0.35 0.35], 'LineWidth', 1);
    ax2.YTick = y;
    ax2.YTickLabel = cellstr(metrics);
    ax2.YDir = 'reverse';
    ax2.Box = 'off';
    ax2.TickDir = 'out';
    ax2.FontName = 'Arial';
    ax2.FontSize = 11;
    xlabel(ax2, '-log_{10}(p)', 'FontName','Arial', 'FontSize',11);
    title(ax2, 'B  Significance before and after outlier removal', ...
        'FontName','Arial', 'FontSize',12, 'FontWeight','bold');
    legend(ax2, {'raw','clean','p=0.05'}, 'Location','best', 'Box','off', 'FontName','Arial', 'FontSize',9);

    sgtitle('Sensitivity analysis of structural-functional coupling', ...
        'FontName','Arial', 'FontSize',13, 'FontWeight','bold');

    out_tif = fullfile(out_dir, 'Figure4_sensitivity_forest.tif');
    out_png = fullfile(out_dir, 'Figure4_sensitivity_forest.png');
    exportgraphics(fig, out_tif, 'Resolution', 600, 'BackgroundColor','white', 'ContentType','image');
    exportgraphics(fig, out_png, 'Resolution', 300, 'BackgroundColor','white', 'ContentType','image');

    fprintf('Save/path/to/project', out_tif, out_png);
end
