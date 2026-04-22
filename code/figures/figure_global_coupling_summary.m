% Repository usage summary
% Description: Create the global structural-functional coupling summary figure.
% Usage: Run after global coupling GLM outputs have been generated.
% Outputs: Exports figure images and an accompanying figure table.
% Note: Update input paths, toolboxes, and filenames for your local environment.

% Public repository version: update file paths, toolboxes, and local settings before running.
% This script/function was lightly sanitized for sharing and may require project-specific inputs.

function out_tif = figure_global_coupling_summary(data_file, coef_file, anova_file, out_dir)
% Figure 5: Global structural-functional coupling main figure
% Left: group distributions
% Right: age/sex-adjusted coefficients and model statistics
% Default export: TIFF (lossless)
%
% Expected files:
%   /path/to/project
%   /path/to/project
%   /path/to/project
%
% Example:
%   figure_global_coupling_summary

    if nargin < 1 || isempty(data_file)
        data_file = '';
    end
    if nargin < 2 || isempty(coef_file)
        coef_file = '';
    end
    if nargin < 3 || isempty(anova_file)
        anova_file = '';
    end
    if nargin < 4 || isempty(out_dir)
        out_dir = '';
    end

    reqs = {
        data_file,  'coupling_covariate_matched_age_sex.txt', '';
        coef_file,  'coupling_ancova_coefficients_age_sex.txt', '';
        anova_file, 'coupling_ancova_table_age_sex.txt', '';
    };

    for i = 1:size(reqs,1)
        if ~isfile(reqs{i,1})
            error(['缺少文件: ' reqs{i,2} newline ...
                   '原路径应为: ' reqs{i,3} newline ...
                   '当前路径: ' reqs{i,1}]);
        end
    end

    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end

    T = readtable(data_file, 'FileType','text', 'Delimiter','\t');
    C = readtable(coef_file, 'FileType','text', 'Delimiter','\t');
    A = readtable(anova_file, 'FileType','text', 'Delimiter','\t');

    T.Properties.VariableNames = matlab.lang.makeValidName(T.Properties.VariableNames);
    C.Properties.VariableNames = matlab.lang.makeValidName(C.Properties.VariableNames);
    A.Properties.VariableNames = matlab.lang.makeValidName(A.Properties.VariableNames);

    reqT = {'group_num','mean_coupling','age','sex'};
    reqC = {'Term','Estimate','SE','tStat','pValue'};
    reqA = {'Term','F','pValue'};
    for i = 1:numel(reqT)
        if ~ismember(reqT{i}, T.Properties.VariableNames), error('数据表缺少列: %s', reqT{i}); end
    end
    for i = 1:numel(reqC)
        if ~ismember(reqC{i}, C.Properties.VariableNames), error('系数表缺少列: %s', reqC{i}); end
    end
    for i = 1:numel(reqA)
        if ~ismember(reqA{i}, A.Properties.VariableNames), error('ANOVA表缺少列: %s', reqA{i}); end
    end

    % Keep valid rows
    ok = isfinite(T.group_num) & isfinite(T.mean_coupling) & isfinite(T.age) & isfinite(T.sex);
    T = T(ok,:);

    % Group labels
    gtxt = strings(height(T),1);
    gtxt(T.group_num==1) = "EOAD";
    gtxt(T.group_num==2) = "LOAD";
    gtxt(T.group_num==3) = "HC";
    keep = gtxt ~= "";
    T = T(keep,:);
    gtxt = gtxt(keep);

    % Means
    y1 = T.mean_coupling(gtxt=="EOAD");
    y2 = T.mean_coupling(gtxt=="LOAD");
    y3 = T.mean_coupling(gtxt=="HC");

    % Fit model again for adjusted means, using the same formula as your pipeline
    D = table(T.mean_coupling, categorical(T.group_num), T.age, T.sex, ...
        'VariableNames', {'y','group','age','sex'});
    mdl = fitlm(D, 'y ~ group + age + sex');

    mu_age = mean(D.age, 'omitnan');
    mu_sex = mean(D.sex, 'omitnan');

    predT = table(categorical([1;2;3]), repmat(mu_age,3,1), repmat(mu_sex,3,1), ...
        'VariableNames', {'group','age','sex'});
    [yhat, yci] = predict(mdl, predT);

    % Read saved model stats
    idx_m = strcmp(string(A.Term), 'Model');
    modelF = NaN; modelP = NaN;
    if any(idx_m)
        modelF = A.F(find(idx_m,1));
        modelP = A.pValue(find(idx_m,1));
    end

    idx2 = strcmp(string(C.Term), 'group_2');
    idx3 = strcmp(string(C.Term), 'group_3');
    beta2 = NaN; se2 = NaN; t2 = NaN; p2 = NaN;
    beta3 = NaN; se3 = NaN; t3 = NaN; p3 = NaN;
    if any(idx2)
        beta2 = C.Estimate(find(idx2,1));
        se2   = C.SE(find(idx2,1));
        t2    = C.tStat(find(idx2,1));
        p2    = C.pValue(find(idx2,1));
    end
    if any(idx3)
        beta3 = C.Estimate(find(idx3,1));
        se3   = C.SE(find(idx3,1));
        t3    = C.tStat(find(idx3,1));
        p3    = C.pValue(find(idx3,1));
    end

    ci2 = [beta2 - 1.96*se2, beta2 + 1.96*se2];
    ci3 = [beta3 - 1.96*se3, beta3 + 1.96*se3];

    % Plot
    fig = figure('Color','w', 'Position', [80 80 1400 540]);
    tl = tiledlayout(fig, 1, 2, 'Padding','compact', 'TileSpacing','compact');

    % ========= Panel A =========
    ax1 = nexttile(tl, 1); hold(ax1, 'on');

    c1 = [0.86 0.39 0.39];
    c2 = [0.34 0.60 0.86];
    c3 = [0.45 0.75 0.50];

    local_scatter_box(ax1, 1, y1, c1);
    local_scatter_box(ax1, 2, y2, c2);
    local_scatter_box(ax1, 3, y3, c3);

    % adjusted means
    for k = 1:3
        plot(ax1, k, yhat(k), 'kd', 'MarkerFaceColor','k', 'MarkerSize',6);
        plot(ax1, [k k], [yci(k,1) yci(k,2)], 'k-', 'LineWidth', 1.2);
    end

    ax1.XLim = [0.5 3.5];
    ax1.XTick = [1 2 3];
    ax1.XTickLabel = {'EOAD','LOAD','HC'};
    ax1.Box = 'off';
    ax1.TickDir = 'out';
    ax1.FontName = 'Arial';
    ax1.FontSize = 11;
    ylabel(ax1, 'Mean structural-functional coupling', 'FontName','Arial', 'FontSize',11);
    title(ax1, 'A  Group distribution and adjusted means', ...
        'FontName','Arial', 'FontSize',12, 'FontWeight','bold');

    % ========= Panel B =========
    ax2 = nexttile(tl, 2); hold(ax2, 'on');

    % contrast plot
    plot(ax2, [0 0], [0.5 2.5], '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 1);
    plot(ax2, ci2, [1 1], '-', 'Color', c2, 'LineWidth', 2);
    scatter(ax2, beta2, 1, 55, 'o', 'MarkerFaceColor', c2, 'MarkerEdgeColor', 'none');
    plot(ax2, ci3, [2 2], '-', 'Color', c3, 'LineWidth', 2);
    scatter(ax2, beta3, 2, 55, 'o', 'o', 'MarkerFaceColor', c3, 'MarkerEdgeColor', 'none');

    ax2.YTick = [1 2];
    ax2.YTickLabel = {'LOAD vs EOAD', 'HC vs EOAD'};
    ax2.YLim = [0.5 2.5];
    ax2.Box = 'off';
    ax2.TickDir = 'out';
    ax2.FontName = 'Arial';
    ax2.FontSize = 11;
    xlabel(ax2, 'Adjusted beta (95% CI)', 'FontName','Arial', 'FontSize',11);
    title(ax2, 'B  Age/sex-adjusted group effects', ...
        'FontName','Arial', 'FontSize',12, 'FontWeight','bold');

    % stats text
    xlim_now = xlim(ax2);
    x_text = xlim_now(1) + 0.05*range(xlim_now);
    text(ax2, x_text, 2.42, sprintf('Model: F = %.2f, p = %s', modelF, local_p(modelP)), ...
        'FontName','Arial', 'FontSize',10, 'FontWeight','bold', 'HorizontalAlignment','left');

    text(ax2, x_text, 1.18, sprintf('\\beta = %.3f, t = %.3f, p = %s', beta2, t2, local_p(p2)), ...
        'FontName','Arial', 'FontSize',10, 'Color', [0.20 0.20 0.20], 'HorizontalAlignment','left');

    text(ax2, x_text, 2.18, sprintf('\\beta = %.3f, t = %.3f, p = %s', beta3, t3, local_p(p3)), ...
        'FontName','Arial', 'FontSize',10, 'Color', [0.20 0.20 0.20], 'HorizontalAlignment','left');

    sgtitle('Global structural-functional coupling', 'FontName','Arial', 'FontSize',13, 'FontWeight','bold');

    out_tif = fullfile(out_dir, 'Figure5_global_coupling_main.tif');
    out_png = fullfile(out_dir, 'Figure5_global_coupling_main.png');
    exportgraphics(fig, out_tif, 'Resolution', 600, 'BackgroundColor','white', 'ContentType','image');
    exportgraphics(fig, out_png, 'Resolution', 300, 'BackgroundColor','white', 'ContentType','image');

    fprintf('Save/path/to/project', out_tif, out_png);
end

function local_scatter_box(ax, xpos, y, c)
    if isempty(y), return; end
    boxchart(ax, ones(size(y))*xpos, y, ...
        'BoxFaceColor', c, ...
        'WhiskerLineColor', c*0.7, ...
        'MarkerStyle', 'none', ...
        'BoxFaceAlpha', 0.55, ...
        'LineWidth', 1.2);
    jitter = 0.10;
    scatter(ax, xpos + (rand(size(y))-0.5)*2*jitter, y, 18, c, 'filled', ...
        'MarkerFaceAlpha', 0.65, 'MarkerEdgeAlpha', 0.18);
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
