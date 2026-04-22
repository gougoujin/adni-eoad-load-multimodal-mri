% Repository usage summary
% Description: Run two-sample SPM analyses for ALFF/fALFF/ReHo/VMHC.
% Usage: Optional voxel-wise analysis pipeline.
% Outputs: Not required for the curated manuscript code set.
% Note: Update input paths, toolboxes, and filenames for your local environment.

% Public repository version: update file paths, toolboxes, and local settings before running.
% This script/function was lightly sanitized for sharing and may require project-specific inputs.

%% ===== 4 个指标（ALFF / fALFF / ReHo / VMHC）两组 t 检验（带协变量）=====
clear; clc;

%% SPM 路径
% addpath('/path/to/required/toolbox');   % 改成你的 SPM 路径
spm('defaults','FMRI');
spm_jobman('initcfg');

%% ---------- 两组根目录 ----------
grp1_root = '';
grp2_root = '';

%% ---------- 指标 & 对应文件夹名 & 文件前缀 ----------
% 请确认下面 3 列是否和你真实情况一致，不一致就改成你自己的

metric_names   = {'ALFF',              'fALFF',               'ReHo',               'VMHC'};
metric_subdirs = {'ALFF_FunImgARCWF',  'fALFF_FunImgARCWF',   'ReHo_FunImgARCWF',   'VMHC_FunImgARCWFsymS'};
metric_prefix  = {'zALFFMap_sub',      'zfALFFMap_sub',       'zReHoMap_sub',       'zVMHCMap_sub'};

%% 协变量（行顺序：先组1，再组2）
cov_file  = '';

%% mask（AAL116）
mask_file = '';

%% 结果输出根目录
out_root  = '';

%% -------- 读协变量 --------
try
    C = readmatrix(cov_file);
catch
    try
        C = table2array(readtable(cov_file));
    catch
        [C,~,~] = xlsread(cov_file, '', '', 'basic');
    end
end
[nSub_all, nCov] = size(C);
fprintf('协变量矩阵：%d 被试 × %d 协变量\n', nSub_all, nCov);

%% -------- 用 ALFF 估计每组被试数量（只做一次） --------
idx_ref   = 1;  % 用第 1 个指标（ALFF）做参考
ref_dir1  = fullfile(grp1_root, metric_subdirs{idx_ref});
ref_dir2  = fullfile(grp2_root, metric_subdirs{idx_ref});
ref_pref  = metric_prefix{idx_ref};

pattern_ref_nii = [ref_pref '*.nii'];  % e.g. zALFFMap_sub*.nii
pattern_ref_img = [ref_pref '*.img'];

g1_ref = dir(fullfile(ref_dir1, pattern_ref_nii));
if isempty(g1_ref), g1_ref = dir(fullfile(ref_dir1, pattern_ref_img)); end

g2_ref = dir(fullfile(ref_dir2, pattern_ref_nii));
if isempty(g2_ref), g2_ref = dir(fullfile(ref_dir2, pattern_ref_img)); end

n1 = numel(g1_ref);
n2 = numel(g2_ref);
fprintf('参考指标 %s：组1 = %d，组2 = %d，总计 = %d\n', metric_names{idx_ref}, n1, n2, n1+n2);

if n1 + n2 ~= nSub_all
    error('协变量行数 (%d) ≠ 组1(%d)+组2(%d)，请检查 Excel 顺序或文件匹配模式。', ...
          nSub_all, n1, n2);
end

%% ===================== 指标循环 =====================
for m = 1:numel(metric_names)
    metric  = metric_names{m};
    subdir1 = metric_subdirs{m};
    subdir2 = metric_subdirs{m};
    prefix  = metric_prefix{m};

    fprintf('\n================ 指标：%s ==================\n', metric);

    grp1_dir = fullfile(grp1_root, subdir1);
    grp2_dir = fullfile(grp2_root, subdir2);

    pattern_nii = [prefix '*.nii'];
    pattern_img = [prefix '*.img'];

    % 组1文件
    g1_files = dir(fullfile(grp1_dir, pattern_nii));
    if isempty(g1_files), g1_files = dir(fullfile(grp1_dir, pattern_img)); end

    % 组2文件
    g2_files = dir(fullfile(grp2_dir, pattern_nii));
    if isempty(g2_files), g2_files = dir(fullfile(grp2_dir, pattern_img)); end

    if isempty(g1_files) || isempty(g2_files)
        warning('%s：有组没有找到图像（dir: %s / %s, pattern: %s），跳过该指标。', ...
                metric, grp1_dir, grp2_dir, prefix);
        continue;
    end

    if numel(g1_files) ~= n1 || numel(g2_files) ~= n2
        warning('%s：图像数量（G1=%d, G2=%d）与参考 n1=%d, n2=%d 不一致，跳过该指标。', ...
                metric, numel(g1_files), numel(g2_files), n1, n2);
        continue;
    end

    % 排序，保证被试顺序一致（按文件名）
    g1_names = sort({g1_files.name});
    g2_names = sort({g2_files.name});

    scans1 = cell(n1,1);
    for i = 1:n1
        scans1{i} = [fullfile(grp1_dir, g1_names{i}) ',1'];
    end

    scans2 = cell(n2,1);
    for i = 1:n2
        scans2{i} = [fullfile(grp2_dir, g2_names{i}) ',1'];
    end

    % 输出目录：/path/to/project 等
    out_dir = fullfile(out_root, metric);
    if ~exist(out_dir,'dir'), mkdir(out_dir); end

    %% =============== 1. 两样本 t 设计 =================
    matlabbatch = [];
    matlabbatch{1}.spm.stats.factorial_design.dir = {out_dir};

    matlabbatch{1}.spm.stats.factorial_design.des.t2.scans1 = scans1;
    matlabbatch{1}.spm.stats.factorial_design.des.t2.scans2 = scans2;
    matlabbatch{1}.spm.stats.factorial_design.des.t2.dept     = 0; % 组间独立
    matlabbatch{1}.spm.stats.factorial_design.des.t2.variance = 1; % 不等方差
    matlabbatch{1}.spm.stats.factorial_design.des.t2.gmsca    = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.t2.ancova   = 0;

    % 协变量
    matlabbatch{1}.spm.stats.factorial_design.cov = struct( ...
        'c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
    for c = 1:nCov
        matlabbatch{1}.spm.stats.factorial_design.cov(c).c     = C(:,c);
        matlabbatch{1}.spm.stats.factorial_design.cov(c).cname = sprintf('cov%d', c);
        matlabbatch{1}.spm.stats.factorial_design.cov(c).iCFI  = 1; % 不与组交互
        matlabbatch{1}.spm.stats.factorial_design.cov(c).iCC   = 1; % 总体中心化
    end

    % mask
    matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.im         = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.em         = {[mask_file ',1']};

    % 全局
    matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit         = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm        = 1;

    %% =============== 2. 参数估计 =================
    matlabbatch{2}.spm.stats.fmri_est.spmmat           = { fullfile(out_dir,'SPM.mat') };
    matlabbatch{2}.spm.stats.fmri_est.write_residuals  = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

    %% =============== 3. 对比：Group1>Group2 =================
    matlabbatch{3}.spm.stats.con.spmmat = { fullfile(out_dir,'SPM.mat') };
    matlabbatch{3}.spm.stats.con.delete = 0;

    matlabbatch{3}.spm.stats.con.consess{1}.tcon.name    = 'Group1>Group2';
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [1 -1 zeros(1,nCov)];
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';

    %% 运行
    try
        spm_jobman('run', matlabbatch);
        fprintf('指标 %s：完成。\n', metric);
    catch ME
        warning('指标 %s 运行出错：%s', metric, ME.message);
    end
end

fprintf('\n=== 全部指标完成，结果在 /path/to/project 下 ===\n');
