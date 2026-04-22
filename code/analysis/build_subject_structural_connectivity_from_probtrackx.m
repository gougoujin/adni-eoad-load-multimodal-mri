% Repository usage summary
% Description: Build a 116x116 structural connectivity matrix for one subject from probtrackx outputs.
% Usage: build_subject_structural_connectivity_from_probtrackx(subj_dir, nROI)
% Outputs: Writes subject-level structural connectivity matrices from tractography targets.
% Note: Update input paths, toolboxes, and filenames for your local environment.

% Public repository version: update file paths, toolboxes, and local settings before running.
% This script/function was lightly sanitized for sharing and may require project-specific inputs.

function build_subject_structural_connectivity_from_probtrackx(subj_dir, nROI)
% 用 seeds_to_targets.nii.gz + roi_dwi/roi_*.nii.gz
% 聚合成 116x116 结构连接矩阵
%
% 兼容：
% - 3D ROI mask
% - 4D ROI mask（第4维=1）
% - 4D ROI stack（第4维=nROI），自动取第 j 个 volume
%
% 用法：
%   build_subject_structural_connectivity_from_probtrackx('', 116)

if nargin < 2
    nROI = 116;
end

tract_dir = fullfile(subj_dir, 'tract');
roi_dir   = fullfile(subj_dir, 'roi_dwi');

if ~exist(tract_dir, 'dir')
    error('找不到 tract 目录: %s', tract_dir);
end
if ~exist(roi_dir, 'dir')
    error('找不到 roi_dwi 目录: %s', roi_dir);
end

% 找一个参考 seeds_to_targets，确定空间尺寸
ref_file = '';
for s = 1:nROI
    tmp = fullfile(tract_dir, sprintf('seed_%d', s), 'seeds_to_targets.nii.gz');
    if exist(tmp, 'file')
        ref_file = tmp;
        break;
    end
end
if isempty(ref_file)
    error('找不到任何 seeds_to_targets.nii.gz，无法确定参考尺寸');
end

Vref = double(niftiread(ref_file));
ref_size = size(Vref);
if numel(ref_size) < 3
    error('参考 seeds_to_targets 不是 3D/4D NIfTI');
end
ref_xyz = ref_size(1:3);

fprintf('参考 seeds_to_targets 尺寸: [%d %d %d]\n', ref_xyz(1), ref_xyz(2), ref_xyz(3));

% 预读取所有 target ROI mask
roi_masks = cell(nROI, 1);
roi_vox   = zeros(nROI, 1);

for j = 1:nROI
    roi_file = fullfile(roi_dir, sprintf('roi_%d.nii.gz', j));

    if ~exist(roi_file, 'file')
        roi_masks{j} = [];
        roi_vox(j) = 0;
        continue;
    end

    R = niftiread(roi_file);
    sz = size(R);

    % 统一转成 3D mask
    if numel(sz) == 3
        M = R;

    elseif numel(sz) == 4
        if sz(4) == 1
            M = R(:,:,:,1);
        elseif sz(4) == nROI
            % 如果 roi_j.nii.gz 实际上装着 116 个 volume，则取第 j 个 volume
            M = R(:,:,:,j);
        else
            error('roi_%d 的第4维=%d，无法判断该取哪一层', j, sz(4));
        end

    else
        error('roi_%d 维度异常: [%s]', j, num2str(sz));
    end

    M = M > 0;

    if ~isequal(size(M), ref_xyz)
        error('roi_%d 尺寸不一致: [%s] vs [%s]', ...
            j, num2str(size(M)), num2str(ref_xyz));
    end

    roi_masks{j} = M;
    roi_vox(j) = nnz(M);
end

SC = zeros(nROI, nROI);
missingSeeds = [];
doneSeeds = [];

for seed = 1:nROI
    seed_dir  = fullfile(tract_dir, sprintf('seed_%d', seed));
    seed_file = fullfile(seed_dir, 'seeds_to_targets.nii.gz');

    if ~exist(seed_file, 'file')
        missingSeeds(end+1) = seed; %#ok<AGROW>
        fprintf('seed_%d 缺失，整行保留 0\n', seed);
        continue;
    end

    try
        V = double(niftiread(seed_file));
    catch ME
        warning('seed_%d: 读取 seeds_to_targets.nii.gz 失败: %s', seed, ME.message);
        missingSeeds(end+1) = seed; %#ok<AGROW>
        continue;
    end

    if numel(size(V)) ~= 3
        error('seed_%d 的 seeds_to_targets 不是 3D: [%s]', seed, num2str(size(V)));
    end

    if ~isequal(size(V), ref_xyz)
        error('seed_%d 的 seeds_to_targets 尺寸异常: [%s] vs [%s]', ...
            seed, num2str(size(V)), num2str(ref_xyz));
    end

    if isempty(V) || nnz(V) == 0
        fprintf('seed_%d: seeds_to_targets 为空，整行保留 0\n', seed);
        missingSeeds(end+1) = seed; %#ok<AGROW>
        continue;
    end

    rowVec = zeros(1, nROI);

    for j = 1:nROI
        M = roi_masks{j};
        if isempty(M) || roi_vox(j) == 0
            continue;
        end
        rowVec(j) = sum(V(M));
    end

    % 自连接清零
    rowVec(seed) = 0;

    SC(seed, :) = rowVec;
    doneSeeds(end+1) = seed; %#ok<AGROW>

    fprintf('seed_%d 完成 | 非零目标数=%d | 总计数=%.0f\n', ...
        seed, nnz(rowVec), sum(rowVec));
end

% 对角线清零
SC(1:nROI+1:end) = 0;

% 对称化
SC_sym_mean = (SC + SC.') / 2;
SC_sym_sum  = SC + SC.';

% 按 target ROI 体素数归一化
SC_roiMean = zeros(size(SC));
for j = 1:nROI
    if roi_vox(j) > 0
        SC_roiMean(:, j) = SC(:, j) / roi_vox(j);
    end
end
SC_roiMean(1:nROI+1:end) = 0;
SC_roiMean_sym = (SC_roiMean + SC_roiMean.') / 2;

% 输出
out_txt         = fullfile(tract_dir, 'SC_counts_116x116.txt');
out_log         = fullfile(tract_dir, 'SC_log1p_116x116.txt');
out_mat         = fullfile(tract_dir, 'SC_counts_116x116.mat');
out_sym_mean    = fullfile(tract_dir, 'SC_counts_116x116_sym_mean.txt');
out_sym_sum     = fullfile(tract_dir, 'SC_counts_116x116_sym_sum.txt');
out_roiMean     = fullfile(tract_dir, 'SC_roiMean_116x116.txt');
out_roiMean_sym = fullfile(tract_dir, 'SC_roiMean_116x116_sym_mean.txt');

writematrix(SC, out_txt, 'Delimiter', 'tab');
writematrix(log1p(SC), out_log, 'Delimiter', 'tab');
writematrix(SC_sym_mean, out_sym_mean, 'Delimiter', 'tab');
writematrix(SC_sym_sum, out_sym_sum, 'Delimiter', 'tab');
writematrix(SC_roiMean, out_roiMean, 'Delimiter', 'tab');
writematrix(SC_roiMean_sym, out_roiMean_sym, 'Delimiter', 'tab');

save(out_mat, 'SC', 'SC_sym_mean', 'SC_sym_sum', ...
    'SC_roiMean', 'SC_roiMean_sym', ...
    'roi_vox', 'doneSeeds', 'missingSeeds', 'subj_dir');

fprintf('\n=== 完成 ===\n');
fprintf('输出文件：\n');
fprintf('%s\n', out_txt);
fprintf('%s\n', out_log);
fprintf('%s\n', out_sym_mean);
fprintf('%s\n', out_sym_sum);
fprintf('%s\n', out_roiMean);
fprintf('%s\n', out_roiMean_sym);
fprintf('%s\n', out_mat);
fprintf('完成 seed 数: %d\n', numel(doneSeeds));
fprintf('缺失 seed 数: %d\n', numel(missingSeeds));
fprintf('矩阵非零元素数: %d\n', nnz(SC));
fprintf('矩阵总和: %.0f\n', sum(SC(:)));
end
