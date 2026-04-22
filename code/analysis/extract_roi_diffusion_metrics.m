% Repository usage summary
% Description: Extract ROI-level diffusion metrics such as FA and MD from subject images.
% Usage: Run after ROI masks and diffusion maps are available for all subjects.
% Outputs: Creates subject-by-ROI matrices and text summaries.
% Note: Update input paths, toolboxes, and filenames for your local environment.

% Public repository version: update file paths, toolboxes, and local settings before running.
% This script/function was lightly sanitized for sharing and may require project-specific inputs.

function extract_roi_diffusion_metrics()
clear; clc;

%% ===================== 路径配置 =====================
proc_root = '';
out_dir   = '';
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

nROI = 116;

out_fa_mat = fullfile(out_dir, 'fa_roi_all_subjects.mat');
out_md_mat = fullfile(out_dir, 'md_roi_all_subjects.mat');
out_fa_txt = fullfile(out_dir, 'fa_roi_summary.txt');
out_md_txt = fullfile(out_dir, 'md_roi_summary.txt');

%% ===================== 收集被试 =====================
rows = {};
fa_all = [];
md_all = [];

groups = dir(fullfile(proc_root, 'group*'));
groups = groups([groups.isdir]);

for g = 1:numel(groups)
    gname = groups(g).name;
    gdir  = fullfile(proc_root, gname);

    subs = dir(fullfile(gdir, 'sub*'));
    subs = subs([subs.isdir]);

    for s = 1:numel(subs)
        sid  = subs(s).name;
        sdir = fullfile(gdir, sid);

        dti_dir = fullfile(sdir, 'dti');
        roi_dir = fullfile(sdir, 'roi_dwi');

        fa_file = fullfile(dti_dir, 'dti_FA.nii.gz');
        md_file = fullfile(dti_dir, 'dti_MD.nii.gz');

        if ~exist(fa_file, 'file') || ~exist(md_file, 'file') || ~exist(roi_dir, 'dir')
            fprintf('[SKIP] 缺少 FA/MD/roi_dwi: %s/%s\n', gname, sid);
            continue;
        end

        try
            FA = double(niftiread(fa_file));
            MD = double(niftiread(md_file));
        catch ME
            fprintf('[SKIP] 读取失败 %s/%s: %s\n', gname, sid, ME.message);
            continue;
        end

        fa_roi = nan(1, nROI);
        md_roi = nan(1, nROI);

        for r = 1:nROI
            roi_file = fullfile(roi_dir, sprintf('roi_%d.nii.gz', r));
            if ~exist(roi_file, 'file')
                continue;
            end

            try
                M = niftiread(roi_file);
            catch
                continue;
            end

            sz = size(M);
            if numel(sz) == 4
                if sz(4) == 1
                    M = M(:,:,:,1);
                elseif sz(4) == nROI
                    M = M(:,:,:,r);
                else
                    warning('ROI 维度异常: %s', roi_file);
                    continue;
                end
            end

            mask = M > 0;

            if ~isequal(size(mask), size(FA)) || ~isequal(size(mask), size(MD))
                warning('空间尺寸不匹配: %s/%s ROI%d', gname, sid, r);
                continue;
            end

            fa_vals = FA(mask);
            md_vals = MD(mask);

            % 去掉无效值
            fa_vals = fa_vals(isfinite(fa_vals) & fa_vals > 0);
            md_vals = md_vals(isfinite(md_vals) & md_vals > 0);

            if ~isempty(fa_vals)
                fa_roi(r) = mean(fa_vals);
            end
            if ~isempty(md_vals)
                md_roi(r) = mean(md_vals);
            end
        end

        rows(end+1, :) = {gname, sid}; %#ok<AGROW>
        fa_all(end+1, :) = fa_roi; %#ok<AGROW>
        md_all(end+1, :) = md_roi; %#ok<AGROW>

        fprintf('[OK] %s/%s\n', gname, sid);
    end
end

%% ===================== 保存 mat =====================
subject_info = rows; %#ok<NASGU>
save(out_fa_mat, 'subject_info', 'fa_all');
save(out_md_mat, 'subject_info', 'md_all');

%% ===================== 保存 txt =====================
write_summary_txt(out_fa_txt, rows, fa_all, 'FA');
write_summary_txt(out_md_txt, rows, md_all, 'MD');

fprintf('\n完成。\nFA: %s\nMD: %s\n', out_fa_txt, out_md_txt);

end

function write_summary_txt(out_txt, rows, X, metric_name)
fid = fopen(out_txt, 'w');
if fid < 0
    error('无法写出文件: %s', out_txt);
end

fprintf(fid, 'group\tsubject');
for i = 1:size(X,2)
    fprintf(fid, '\t%s_ROI%03d', metric_name, i);
end
fprintf(fid, '\n');

for i = 1:size(X,1)
    fprintf(fid, '%s\t%s', rows{i,1}, rows{i,2});
    for j = 1:size(X,2)
        v = X(i,j);
        if isnan(v)
            fprintf(fid, '\tNaN');
        else
            fprintf(fid, '\t%.6f', v);
        end
    end
    fprintf(fid, '\n');
end

fclose(fid);
end
