% Repository usage summary
% Description: Export AAL ROI masks and BrainNet nodes for selected regions.
% Usage: Optional visualization helper.
% Outputs: Useful for mask/node generation outside the main pipeline.
% Note: Update input paths, toolboxes, and filenames for your local environment.

% Public repository version: update file paths, toolboxes, and local settings before running.
% This script/function was lightly sanitized for sharing and may require project-specific inputs.

function export_aal_roi_masks_and_nodes()
% 从 AAL116 图谱中提取显著 ROI
% 输出：
% 1) 每个ROI单独mask: ROI_xxx_name.nii
% 2) 多标签mask: significant_rois_multilabel.nii
% 3) BrainNet节点文件: significant_rois.node
% 4) 信息表: significant_rois_info.csv

clc; clear;

%% ===================== 用户需要修改的部分 =====================
aal_nii = '';
out_dir = '';
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

% 显著ROI编号
roi_ids = [3, 7, 15, 23];

% ROI名称
roi_names = {
    'Frontal_Sup_Medial_L'
    'Frontal_Mid_R'
    'Frontal_Mid_Orb_R'
    'Frontal_Inf_Tri_R'
};

% 效应值（建议 beta）
roi_beta = [0.18, 0.21, 0.19, 0.17];

% 显著性（建议 -log10(q)）
roi_logq = [1.8, 2.1, 1.9, 1.7];
%% ============================================================

if numel(roi_ids) ~= numel(roi_names) || numel(roi_ids) ~= numel(roi_beta) || numel(roi_ids) ~= numel(roi_logq)
    error('roi_ids, roi_names, roi_beta, roi_logq 长度必须一致');
end

%% 读取图谱
V = spm_vol(aal_nii);
Y = spm_read_vols(V);

nroi = numel(roi_ids);
node_xyz = zeros(nroi, 3);

% 多标签图：1,2,3,4... 对应不同ROI
multilabel = zeros(size(Y));

for i = 1:nroi
    roi_val = roi_ids(i);
    roi_mask = (Y == roi_val);

    idx = find(roi_mask);
    if isempty(idx)
        warning('ROI %d 在图谱中没有找到体素', roi_val);
        node_xyz(i, :) = [0 0 0];
        continue;
    end

    % ---------- 1. 单独ROI mask ----------
    single_mask = double(roi_mask);
    Vsingle = V;
    safe_name = regexprep(roi_names{i}, '[^a-zA-Z0-9_]', '_');
    Vsingle.fname = fullfile(out_dir, sprintf('ROI_%03d_%s.nii', roi_val, safe_name));
    spm_write_vol(Vsingle, single_mask);

    % ---------- 2. 多标签mask ----------
    multilabel(roi_mask) = i;

    % ---------- 3. 计算质心 ----------
    [x, y, z] = ind2sub(size(Y), idx);
    vox = [x y z ones(numel(x),1)]';
    mni = V.mat * vox;
    node_xyz(i, :) = mean(mni(1:3, :), 2)';
end

%% 保存多标签mask
Vmulti = V;
Vmulti.fname = fullfile(out_dir, 'significant_rois_multilabel.nii');
spm_write_vol(Vmulti, multilabel);

%% 为了让 BrainNet 颜色更明显，对 beta 做线性拉伸到 [1, 10]
beta_min = min(roi_beta);
beta_max = max(roi_beta);
if abs(beta_max - beta_min) < eps
    color_scaled = repmat(5, size(roi_beta));
else
    color_scaled = 1 + 9 * (roi_beta - beta_min) ./ (beta_max - beta_min);
end

% size 也做一个温和缩放
size_scaled = 2 + 4 * (roi_logq - min(roi_logq)) ./ max(eps, (max(roi_logq)-min(roi_logq)));

%% 保存 .node 文件
node_file = fullfile(out_dir, 'significant_rois.node');
fid = fopen(node_file, 'w');
if fid == -1
    error('无法创建 node 文件');
end

for i = 1:nroi
    fprintf(fid, '%.3f\t%.3f\t%.3f\t%.6f\t%.6f\t%s\n', ...
        node_xyz(i,1), node_xyz(i,2), node_xyz(i,3), ...
        color_scaled(i), size_scaled(i), roi_names{i});
end
fclose(fid);

%% 保存信息表
T = table(roi_ids(:), roi_names(:), ...
    node_xyz(:,1), node_xyz(:,2), node_xyz(:,3), ...
    roi_beta(:), roi_logq(:), color_scaled(:), size_scaled(:), ...
    'VariableNames', {'ROI_ID','ROI_Name','X','Y','Z','Beta','Log10Q','NodeColor','NodeSize'});
writetable(T, fullfile(out_dir, 'significant_rois_info.csv'));

fprintf('完成。\n');
fprintf('每个ROI单独nii已输出到: %s\n', out_dir);
fprintf('多标签nii: %s\n', fullfile(out_dir, 'significant_rois_multilabel.nii'));
fprintf('node文件: %s\n', node_file);
fprintf('信息表: %s\n', fullfile(out_dir, 'significant_rois_info.csv'));
end
