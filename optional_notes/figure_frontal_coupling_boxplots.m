% Repository usage summary
% Description: Alternative frontal boxplot figure kept for reference.
% Usage: Optional figure variant.
% Outputs: Not required if the violin-panel figure is used.
% Note: Update input paths, toolboxes, and filenames for your local environment.

% Public repository version: update file paths, toolboxes, and local settings before running.
% This script/function was lightly sanitized for sharing and may require project-specific inputs.

function figure_frontal_coupling_boxplots()
% Plot 4 frontal ROI group distributions for EOAD / LOAD / HC
% Suitable for combining with the BrainNet node figure as a main manuscript figure.
%
% INPUT:
%   A CSV or XLSX file with one row per subject and at least these columns:
%   - Group   : 'EOAD', 'LOAD', 'HC'  (or 1,2,3; see mapping below)
%   - ROI1    : value for frontal ROI 1
%   - ROI2    : value for frontal ROI 2
%   - ROI3    : value for frontal ROI 3
%   - ROI4    : value for frontal ROI 4
%
% OUTPUT:
%   A high-resolution TIFF and PNG figure.
%
% MATLAB version: R2023a+

clc; clear; close all;

%% ===================== USER SETTINGS =====================
data_file = '';  % <- change to your file
sheet_name = 1;                                       % only for xlsx

out_dir = '';
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

% Column names in your file
group_col = 'Group';
roi_cols = {'ROI1', 'ROI2', 'ROI3', 'ROI4'};

% Titles shown above each panel
roi_titles = {
    'Right middle frontal gyrus'
    'Right orbital middle frontal gyrus'
    'Right inferior frontal triangular gyrus'
    'Left medial superior frontal gyrus'
};

% If Group is numeric, map here:
% 1 = EOAD, 2 = LOAD, 3 = HC
numeric_group_map = {'EOAD','LOAD','HC'};

% Colors (soft, manuscript-friendly)
c_eoad = [0.85 0.40 0.40];
c_load = [0.35 0.60 0.85];
c_hc   = [0.45 0.75 0.50];

% Y-axis label
y_label_text = 'Structural-functional coupling';

% Output name
out_base = 'Figure_frontal_roi_group_boxplots';
%% =======================================================

%% Read table
[~,~,ext] = fileparts(data_file);
if strcmpi(ext, '.csv')
    T = readtable(data_file);
else
    T = readtable(data_file, 'Sheet', sheet_name);
end

% Basic checks
vars_needed = [{group_col}, roi_cols];
for i = 1:numel(vars_needed)
    if ~ismember(vars_needed{i}, T.Properties.VariableNames)
        error('Missing required column: %s', vars_needed{i});
    end
end

%% Normalize group labels
G = T.(group_col);

if isnumeric(G)
    gtxt = strings(size(G));
    for i = 1:numel(G)
        if ismember(G(i), [1 2 3])
            gtxt(i) = numeric_group_map{G(i)};
        else
            gtxt(i) = "";
        end
    end
elseif iscell(G) || isstring(G) || iscategorical(G)
    gtxt = string(G);
else
    error('Unsupported group column type.');
end

gtxt = strtrim(upper(gtxt));
gtxt(gtxt=="1") = "EOAD";
gtxt(gtxt=="2") = "LOAD";
gtxt(gtxt=="3") = "HC";

valid = ismember(gtxt, ["EOAD","LOAD","HC"]);
T = T(valid,:);
gtxt = gtxt(valid);

%% Figure
fig = figure('Color','w','Position',[100 100 1200 320]);

for k = 1:4
    ax = subplot(1,4,k); hold(ax,'on');

    y = T.(roi_cols{k});
    y = double(y);

    idx1 = gtxt=="EOAD" & ~isnan(y);
    idx2 = gtxt=="LOAD" & ~isnan(y);
    idx3 = gtxt=="HC"   & ~isnan(y);

    y1 = y(idx1);
    y2 = y(idx2);
    y3 = y(idx3);

    % Boxplots
    x1 = ones(size(y1))*1;
    x2 = ones(size(y2))*2;
    x3 = ones(size(y3))*3;

    boxchart(x1, y1, 'BoxFaceColor', c_eoad, 'WhiskerLineColor', c_eoad*0.7, ...
        'MarkerStyle','none', 'BoxFaceAlpha',0.55, 'LineWidth',1.2);
    boxchart(x2, y2, 'BoxFaceColor', c_load, 'WhiskerLineColor', c_load*0.7, ...
        'MarkerStyle','none', 'BoxFaceAlpha',0.55, 'LineWidth',1.2);
    boxchart(x3, y3, 'BoxFaceColor', c_hc, 'WhiskerLineColor', c_hc*0.7, ...
        'MarkerStyle','none', 'BoxFaceAlpha',0.55, 'LineWidth',1.2);

    % Jittered points
    jitter = 0.08;
    scatter(1 + (rand(size(y1))-0.5)*2*jitter, y1, 20, c_eoad, 'filled', ...
        'MarkerFaceAlpha',0.70, 'MarkerEdgeAlpha',0.20);
    scatter(2 + (rand(size(y2))-0.5)*2*jitter, y2, 20, c_load, 'filled', ...
        'MarkerFaceAlpha',0.70, 'MarkerEdgeAlpha',0.20);
    scatter(3 + (rand(size(y3))-0.5)*2*jitter, y3, 20, c_hc, 'filled', ...
        'MarkerFaceAlpha',0.70, 'MarkerEdgeAlpha',0.20);

    % Mean markers
    if ~isempty(y1), plot(1, mean(y1,'omitnan'), 'kd', 'MarkerFaceColor','k', 'MarkerSize',5); end
    if ~isempty(y2), plot(2, mean(y2,'omitnan'), 'kd', 'MarkerFaceColor','k', 'MarkerSize',5); end
    if ~isempty(y3), plot(3, mean(y3,'omitnan'), 'kd', 'MarkerFaceColor','k', 'MarkerSize',5); end

    ax.XLim = [0.5 3.5];
    ax.XTick = [1 2 3];
    ax.XTickLabel = {'EOAD','LOAD','HC'};
    ax.TickDir = 'out';
    ax.Box = 'off';
    ax.LineWidth = 1;
    ax.FontName = 'Arial';
    ax.FontSize = 10;

    title(roi_titles{k}, 'FontName','Arial', 'FontSize',10, 'FontWeight','bold');
    if k == 1
        ylabel(y_label_text, 'FontName','Arial', 'FontSize',11);
    end

    % Panel letters
    text(0.03, 0.96, char('A'+k-1), 'Units','normalized', ...
        'FontName','Arial', 'FontSize',12, 'FontWeight','bold', ...
        'HorizontalAlignment','left', 'VerticalAlignment','top');
end

sgtitle('Frontal regional structural-functional coupling', ...
    'FontName','Arial', 'FontSize',12, 'FontWeight','bold');

% Save
png_file = fullfile(out_dir, [out_base '.png']);
tif_file = fullfile(out_dir, [out_base '.tif']);

exportgraphics(fig, png_file, 'Resolution', 300);
exportgraphics(fig, tif_file, 'Resolution', 600);

fprintf('Save/path/to/project', png_file, tif_file);
end
