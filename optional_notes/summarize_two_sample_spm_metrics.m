% Repository usage summary
% Description: Summarize voxel-wise two-sample SPM results with AAL lookup.
% Usage: Optional voxel-wise summary pipeline.
% Outputs: Not required for the curated manuscript code set.
% Note: Update input paths, toolboxes, and filenames for your local environment.

% Public repository version: update file paths, toolboxes, and local settings before running.
% This script/function was lightly sanitized for sharing and may require project-specific inputs.

%% ===== 4 个指标（ALFF / fALFF / ReHo / VMHC）UN/FWE/FDR + AAL 汇总 =====
clear; clc;
% addpath('/path/to/required/toolbox');        % SPM
% addpath('/path/to/required/toolbox');         % lookup_aal_label, AAL 文件

metrics = {'ALFF','fALFF','ReHo','VMHC'};

root_dir      = '';               % 各指标结果根目录
extent_thresh = 10;                     % k>=10
STAT          = 'T';

p_unc = 0.001;   % UN
p_FWE = 0.05;    % voxel-level FWE
p_FDR = 0.05;    % voxel-level FDR

aal_mask_file  = '';
aal_label_file = '';

out_txt = fullfile(root_dir, 'Metrics_TwoSample_UN_FWE_FDR_AAL_Summary.txt');
fid = fopen(out_txt,'w');
if fid==-1, error('无法创建输出文件'); end

%% ---- 读 AAL 掩模 ----
Vaal   = spm_vol(aal_mask_file);
AALvol = spm_read_vols(Vaal);
[dimX,dimY,dimZ] = size(AALvol);

%% ---- 读 AAL 名称表 ----
AAL_names = containers.Map('KeyType','double','ValueType','char');
if exist(aal_label_file,'file')
    try
        Tlab = readtable(aal_label_file, ...
                         'FileType','text', ...
                         'ReadVariableNames',false, ...
                         'Delimiter',{' ','\t',','}, ...
                         'MultipleDelimsAsOne',true);
        for i = 1:size(Tlab,1)
            id   = Tlab{i,1};
            name = Tlab{i,2};
            if iscell(id),   id   = id{1};   end
            if iscell(name), name = name{1}; end
            AAL_names(double(id)) = char(name);
        end
        fprintf('已从 %s 读取 %d 个 AAL 标签名称。\n', aal_label_file, AAL_names.Count);
    catch ME
        fprintf('读取 AAL 标签文件失败：%s\n', ME.message);
        fprintf('仅输出 AAL_ID，名称用 "AAL_数字"。\n');
    end
else
    fprintf('未找到 AAL 标签文件 %s，仅输出 AAL_ID，名称用 "AAL_数字"。\n', aal_label_file);
end

%% ---- 输出表头：多了 Metric 列 ----
fprintf(fid, ['Metric\tCorrection\tClusterID\t',...
              'PeakX\tPeakY\tPeakZ\tTvalue\tP_uncorr_1tail\tNvoxels\tAAL_ID\tAAL_Name\n']);

%% ===================== 指标循环 =====================
for m = 1:numel(metrics)
    metric  = metrics{m};
    out_dir = fullfile(root_dir, metric);
    spm_mat = fullfile(out_dir,'SPM.mat');
    Tfile   = fullfile(out_dir,'spmT_0001.nii');

    if ~exist(spm_mat,'file') || ~exist(Tfile,'file')
        fprintf('\n指标 %s: 缺 SPM.mat 或 spmT_0001.nii，跳过。\n', metric);
        continue;
    end

    fprintf('\n========== 指标 %s ==========\n', metric);

    % 自由度 & 体积信息
    load(spm_mat,'SPM');
    if ~isfield(SPM,'xX') || ~isfield(SPM.xX,'erdf')
        fprintf('  %s: 无 SPM.xX.erdf，跳过。\n', metric);
        continue;
    end
    df_t = SPM.xX.erdf;

    if ~isfield(SPM,'xVol') || ~isfield(SPM.xVol,'R') || ~isfield(SPM.xVol,'S')
        fprintf('  %s: SPM.xVol 信息不全，跳过。\n', metric);
        continue;
    end
    R = SPM.xVol.R;
    S = SPM.xVol.S;

    % 读 T 图
    V = spm_vol(Tfile);
    Y = spm_read_vols(V);
    Yvec = Y(:);

    % 有效体素
    mask_valid = ~isnan(Yvec);
    idx_valid  = find(mask_valid);
    if isempty(idx_valid)
        fprintf('  %s: 有效体素为 0，跳过。\n', metric);
        continue;
    end
    Yv = Yvec(mask_valid);

    % 只看正向（Group1>Group2）
    pos_mask = Yv > 0;
    if ~any(pos_mask)
        fprintf('  %s: 所有 T<=0，无正向效应。\n', metric);
        continue;
    end
    idx_pos = idx_valid(pos_mask);
    T_pos   = Yv(pos_mask);

    % 一侧未校正 p
    p_pos_1tail = 1 - tcdf(T_pos, df_t);

    %% --- UN：p<0.001 ---
    Tthr_unc = tinv(1 - p_unc, df_t);
    fprintf('  %s: UN 阈值 T>%.3f (p<%.3g)\n', metric, Tthr_unc, p_unc);

    idx_thr_unc = idx_pos(T_pos > Tthr_unc);
    if ~isempty(idx_thr_unc)
        [ii,jj,kk] = ind2sub(V.dim, idx_thr_unc);
        XYZ_vox = [ii'; jj'; kk'];
        clusters  = spm_clusters(XYZ_vox);
        uClusters = unique(clusters);
        clusterCount = 0;

        for ci = 1:numel(uClusters)
            cid = uClusters(ci);
            vox_idx_cluster = find(clusters == cid);
            nVox = numel(vox_idx_cluster);
            if nVox < extent_thresh, continue; end

            vox_global_idx = idx_thr_unc(vox_idx_cluster);
            Tvals_cluster  = Yvec(vox_global_idx);

            [Tmax, imax]    = max(Tvals_cluster);
            peak_idx_global = vox_global_idx(imax);
            [pi,pj,pk] = ind2sub(V.dim, peak_idx_global);
            peak_mm = V.mat * [pi;pj;pk;1];

            p_uncorr = 1 - tcdf(Tmax, df_t);

            [aal_id, aal_name] = lookup_aal_label(peak_mm, Vaal, AALvol, ...
                                               dimX,dimY,dimZ, AAL_names);

            clusterCount = clusterCount + 1;
            fprintf(fid, '%s\tUN\t%d\t%.2f\t%.2f\t%.2f\t%.4f\t%.4g\t%d\t%d\t%s\n', ...
                metric, clusterCount, ...
                peak_mm(1), peak_mm(2), peak_mm(3), ...
                Tmax, p_uncorr, nVox, ...
                aal_id, aal_name);
        end
        if clusterCount==0
            fprintf('  %s: UN 有体素过阈，但 cluster<%d。\n', metric, extent_thresh);
        end
    else
        fprintf('  %s: UN 无体素通过。\n', metric);
    end

    %% --- FWE：voxel-level p<0.05 ---
    try
        df_vec   = [1 df_t];              % T 检验 df = [1 df_t]
        Tthr_FWE = spm_uc(p_FWE, df_vec, 'T', R, 1, S);
    catch ME
        fprintf('  %s: spm_uc(FWE) 失败：%s\n', metric, ME.message);
        Tthr_FWE = Inf;
    end
    fprintf('  %s: FWE T-threshold ≈ %.3f (p<%.2f)\n', metric, Tthr_FWE, p_FWE);

    idx_thr_fwe = idx_pos(T_pos > Tthr_FWE);
    if ~isempty(idx_thr_fwe) && isfinite(Tthr_FWE)
        [ii,jj,kk] = ind2sub(V.dim, idx_thr_fwe);
        XYZ_vox = [ii'; jj'; kk'];
        clusters  = spm_clusters(XYZ_vox);
        uClusters = unique(clusters);
        clusterCount = 0;

        for ci = 1:numel(uClusters)
            cid = uClusters(ci);
            vox_idx_cluster = find(clusters == cid);
            nVox = numel(vox_idx_cluster);
            if nVox < extent_thresh, continue; end

            vox_global_idx = idx_thr_fwe(vox_idx_cluster);
            Tvals_cluster  = Yvec(vox_global_idx);

            [Tmax, imax]    = max(Tvals_cluster);
            peak_idx_global = vox_global_idx(imax);
            [pi,pj,pk] = ind2sub(V.dim, peak_idx_global);
            peak_mm = V.mat * [pi;pj;pk;1];

            p_uncorr = 1 - tcdf(Tmax, df_t);

            [aal_id, aal_name] = lookup_aal_label(peak_mm, Vaal, AALvol, ...
                                               dimX,dimY,dimZ, AAL_names);

            clusterCount = clusterCount + 1;
            fprintf(fid, '%s\tFWE\t%d\t%.2f\t%.2f\t%.2f\t%.4f\t%.4g\t%d\t%d\t%s\n', ...
                metric, clusterCount, ...
                peak_mm(1), peak_mm(2), peak_mm(3), ...
                Tmax, p_uncorr, nVox, ...
                aal_id, aal_name);
        end
        if clusterCount==0
            fprintf('  %s: FWE 有体素过阈，但 cluster<%d。\n', metric, extent_thresh);
        end
    else
        fprintf('  %s: FWE 无体素通过。\n', metric);
    end

    %% --- FDR：voxel-level p<0.05 ---
    p_all = p_pos_1tail;
    [p_sorted, ~] = sort(p_all);
    mNum  = numel(p_sorted);
    crit  = (1:mNum)' * (p_FDR / mNum);    % BH
    below = find(p_sorted <= crit);

    if isempty(below)
        Tthr_FDR = Inf;
        fprintf('  %s: FDR 无体素 p<%.2f。\n', metric, p_FDR);
        idx_thr_fdr = [];
    else
        p_star   = p_sorted(max(below));
        Tthr_FDR = tinv(1 - p_star, df_t);
        fprintf('  %s: FDR T-threshold ≈ %.3f (p*≈%.4g)\n', ...
                metric, Tthr_FDR, p_star);
        idx_thr_fdr = idx_pos(T_pos > Tthr_FDR);
    end

    if ~isempty(idx_thr_fdr) && isfinite(Tthr_FDR)
        [ii,jj,kk] = ind2sub(V.dim, idx_thr_fdr);
        XYZ_vox = [ii'; jj'; kk'];
        clusters  = spm_clusters(XYZ_vox);
        uClusters = unique(clusters);
        clusterCount = 0;

        for ci = 1:numel(uClusters)
            cid = uClusters(ci);
            vox_idx_cluster = find(clusters == cid);
            nVox = numel(vox_idx_cluster);
            if nVox < extent_thresh, continue; end

            vox_global_idx = idx_thr_fdr(vox_idx_cluster);
            Tvals_cluster  = Yvec(vox_global_idx);

            [Tmax, imax]    = max(Tvals_cluster);
            peak_idx_global = vox_global_idx(imax);
            [pi,pj,pk] = ind2sub(V.dim, peak_idx_global);
            peak_mm = V.mat * [pi;pj;pk;1];

            p_uncorr = 1 - tcdf(Tmax, df_t);

            [aal_id, aal_name] = lookup_aal_label(peak_mm, Vaal, AALvol, ...
                                               dimX,dimY,dimZ, AAL_names);

            clusterCount = clusterCount + 1;
            fprintf(fid, '%s\tFDR\t%d\t%.2f\t%.2f\t%.2f\t%.4f\t%.4g\t%d\t%d\t%s\n', ...
                metric, clusterCount, ...
                peak_mm(1), peak_mm(2), peak_mm(3), ...
                Tmax, p_uncorr, nVox, ...
                aal_id, aal_name);
        end
        if clusterCount==0
            fprintf('  %s: FDR 有体素过阈，但 cluster<%d。\n', metric, extent_thresh);
        end
    end
end

fclose(fid);
fprintf('\n=== 汇总完成，结果写入：%s ===\n', out_txt);
