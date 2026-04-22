% Repository usage summary
% Description: Create the frontal regional coupling violin-panel figure used in the manuscript.
% Usage: Run after regional coupling and covariate files are available.
% Outputs: Exports manuscript-ready PNG/TIFF panels plus plot tables.
% Note: Update input paths, toolboxes, and filenames for your local environment.

% Public repository version: update file paths, toolboxes, and local settings before running.
% This script/function was lightly sanitized for sharing and may require project-specific inputs.

function figure_frontal_coupling_violin(data_dir, out_dir)
% MAKE_FIG1_FRONTAL_VIOLIN
% Build a manuscript-style 4-panel frontal ROI figure from existing project outputs.
%
% Required files in data_dir:
%   - roi_level_coupling_stats_age_sex.mat   (contains roi_coupling, group_num, sub_id, Cm)
%       OR roi_level_sc_fc_coupling_all.mat  (contains roi_coupling_all, matched)
%   - biao.xlsx
%
% Output:
%   - frontal_roi_values_for_plot.csv
%   - Figure1_panelB_frontal_violin.png
%   - Figure1_panelB_frontal_violin.tif
%
% Default ROIs follow the current significant frontal set:
%   AAL 007 Frontal_Mid_R
%   AAL 015 Frontal_Mid_Orb_R
%   AAL 023 Frontal_Inf_Tri_R
%   AAL 003 Frontal_Sup_Medial_L
%
% Example:
%   figure_frontal_coupling_violin('', '')
%
% MATLAB: R2023a+

if nargin < 1 || isempty(data_dir)
    data_dir = pwd;
end
if nargin < 2 || isempty(out_dir)
    out_dir = fullfile(data_dir, 'paper_figures');
end
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

roi_ids = [7 15 23 3];
roi_varnames = {'ROI1','ROI2','ROI3','ROI4'};
roi_titles = {
    'Right middle frontal gyrus'
    'Right orbital middle frontal gyrus'
    'Right inferior frontal triangular gyrus'
    'Left medial superior frontal gyrus'
};

c_eoad = [0.84 0.42 0.42];
c_load = [0.32 0.58 0.86];
c_hc   = [0.43 0.73 0.49];

stats_mat_1 = fullfile(data_dir, 'roi_level_coupling_stats_age_sex.mat');
stats_mat_2 = fullfile(data_dir, 'roi_level_sc_fc_coupling_all.mat');
cov_file    = fullfile(data_dir, 'biao.xlsx');

if ~exist(cov_file, 'file')
    error('Missing biao.xlsx: %s', cov_file);
end

% -------- load ROI matrix --------
if exist(stats_mat_1, 'file')
    S = load(stats_mat_1);
    if isfield(S, 'roi_coupling')
        X = S.roi_coupling; % nROI x nSub
    else
        error('roi_level_coupling_stats_age_sex.mat does not contain roi_coupling');
    end

    if isfield(S, 'group_num')
        group_num = double(S.group_num(:));
    else
        error('roi_level_coupling_stats_age_sex.mat missing group_num');
    end

    if isfield(S, 'sub_id')
        sub_id = string(S.sub_id(:));
    else
        error('roi_level_coupling_stats_age_sex.mat missing sub_id');
    end

elseif exist(stats_mat_2, 'file')
    S = load(stats_mat_2);
    if isfield(S, 'roi_coupling_all')
        X = S.roi_coupling_all;
    else
        error('roi_level_sc_fc_coupling_all.mat does not contain roi_coupling_all');
    end
    if ~isfield(S, 'matched')
        error('roi_level_sc_fc_coupling_all.mat missing matched');
    end

    matched = S.matched;
    nSub = size(matched,1);
    group_num = nan(nSub,1);
    sub_id = strings(nSub,1);
    for i = 1:nSub
        g = string(matched{i,1});
        tok = regexp(char(g), 'group(\d+)', 'tokens', 'once');
        if ~isempty(tok)
            group_num(i) = str2double(tok{1});
        end
        sub_id(i) = string(matched{i,2});
    end
else
    error('Need roi_level_coupling_stats_age_sex.mat or roi_level_sc_fc_coupling_all.mat in %s', data_dir);
end

if size(X,1) < max(roi_ids)
    error('ROI matrix has fewer rows than required ROI IDs.');
end

% -------- read covariates and match subjects --------
raw = readcell(cov_file);
raw = raw(2:end, :);
n = size(raw,1);

cov_group = nan(n,1);
cov_sub   = strings(n,1);
for i = 1:n
    g = raw{i,1};
    if isnumeric(g)
        cov_group(i) = g;
    else
        tok = regexp(char(string(g)), '(\d+)', 'tokens', 'once');
        if ~isempty(tok)
            cov_group(i) = str2double(tok{1});
        end
    end

    s = raw{i,3};
    if isnumeric(s)
        if ~isnan(s)
            cov_sub(i) = sprintf('sub%03d', round(s));
        end
    else
        tok = regexp(char(string(s)), '(\d+)', 'tokens', 'once');
        if ~isempty(tok)
            cov_sub(i) = sprintf('sub%03d', str2double(tok{1}));
        end
    end
end

Tcov = table(cov_group, cov_sub, 'VariableNames', {'group_num','sub'});
Tcov.tag = string(Tcov.group_num) + "_" + Tcov.sub;
Tcov = Tcov(~isnan(Tcov.group_num) & strlength(Tcov.sub) > 0, :);

Timg = table(group_num, sub_id, 'VariableNames', {'group_num','sub'});
Timg.tag = string(Timg.group_num) + "_" + Timg.sub;

[common_tag, ia, ib] = intersect(Timg.tag, Tcov.tag, 'stable'); %#ok<ASGLU>
X = X(:, ia);
Timg = Timg(ia,:);
Tcov = Tcov(ib,:);

% -------- build plotting table --------
Tplot = table;
Tplot.Group = repmat("", numel(common_tag), 1);
for i = 1:height(Tplot)
    switch Timg.group_num(i)
        case 1
            Tplot.Group(i) = "EOAD";
        case 2
            Tplot.Group(i) = "LOAD";
        case 3
            Tplot.Group(i) = "HC";
        otherwise
            Tplot.Group(i) = "";
    end
end
Tplot.sub = Timg.sub;

for k = 1:numel(roi_ids)
    Tplot.(roi_varnames{k}) = X(roi_ids(k), :)';
end

writetable(Tplot, fullfile(out_dir, 'frontal_roi_values_for_plot.csv'));

% -------- plot --------
fig = figure('Color','w', 'Position', [100 100 1360 360]);
letters = 'BCDE';
all_y = [];
for k = 1:4
    all_y = [all_y; Tplot.(roi_varnames{k})]; %#ok<AGROW>
end
all_y = all_y(isfinite(all_y));
ymin = min(all_y);
ymax = max(all_y);
pad = 0.08 * max(1e-6, ymax - ymin);
ylim_use = [ymin - pad, ymax + pad];

for k = 1:4
    ax = subplot(1,4,k); hold(ax,'on');
    y = double(Tplot.(roi_varnames{k}));
    g = string(Tplot.Group);

    vals = {y(g=="HC" & isfinite(y)), y(g=="EOAD" & isfinite(y)), y(g=="LOAD" & isfinite(y))};
    cols = {c_hc, c_eoad, c_load};
    xpos = [1 2 3];

    for j = 1:3
        v = vals{j};
        if isempty(v)
            continue;
        end

        try
            [f, yi] = ksdensity(v, 'NumPoints', 200);
            if max(f) > 0
                f = f ./ max(f) * 0.28;
                patch([xpos(j)-f, fliplr(xpos(j)+f)], [yi, fliplr(yi)], cols{j}, ...
                    'FaceAlpha', 0.18, 'EdgeColor', 'none');
            end
        catch
            % fallback: skip violin if ksdensity fails
        end

        q1 = prctile(v,25);
        q2 = median(v);
        q3 = prctile(v,75);
        iqr_v = q3 - q1;
        loww = max(min(v), q1 - 1.5*iqr_v);
        uppw = min(max(v), q3 + 1.5*iqr_v);

        rectangle('Position', [xpos(j)-0.14, q1, 0.28, max(eps, q3-q1)], ...
            'FaceColor', [cols{j} 0.35], 'EdgeColor', cols{j}*0.75, 'LineWidth', 1.0);
        plot([xpos(j)-0.14 xpos(j)+0.14], [q2 q2], 'k-', 'LineWidth', 1.4);
        plot([xpos(j) xpos(j)], [loww q1], '-', 'Color', cols{j}*0.75, 'LineWidth', 1.0);
        plot([xpos(j) xpos(j)], [q3 uppw], '-', 'Color', cols{j}*0.75, 'LineWidth', 1.0);
        plot([xpos(j)-0.07 xpos(j)+0.07], [loww loww], '-', 'Color', cols{j}*0.75, 'LineWidth', 1.0);
        plot([xpos(j)-0.07 xpos(j)+0.07], [uppw uppw], '-', 'Color', cols{j}*0.75, 'LineWidth', 1.0);

        xj = xpos(j) + (rand(size(v))-0.5)*0.18;
        scatter(xj, v, 18, cols{j}, 'filled', 'MarkerFaceAlpha', 0.58, 'MarkerEdgeAlpha', 0.12);
        plot(xpos(j), mean(v,'omitnan'), 'kd', 'MarkerFaceColor', 'k', 'MarkerSize', 5);
    end

    xlim([0.5 3.5]);
    ylim(ylim_use);
    xticks([1 2 3]);
    xticklabels({'HC','EOAD','LOAD'});
    set(ax, 'TickDir','out', 'Box','off', 'LineWidth',1, 'FontName','Arial', 'FontSize',10);
    title(roi_titles{k}, 'FontName','Arial', 'FontWeight','bold', 'FontSize',10);
    if k == 1
        ylabel('Structural-functional coupling', 'FontName','Arial', 'FontSize',11);
    end
    text(0.01, 0.98, letters(k), 'Units','normalized', 'FontSize',13, ...
        'FontWeight','bold', 'HorizontalAlignment','left', 'VerticalAlignment','top');
end

exportgraphics(fig, fullfile(out_dir, 'Figure1_panelB_frontal_violin.png'), 'Resolution', 600);
exportgraphics(fig, fullfile(out_dir, 'Figure1_panelB_frontal_violin.tif'), 'Resolution', 600);
close(fig);

fprintf('Done. Output folde/path/to/project', out_dir);
end
