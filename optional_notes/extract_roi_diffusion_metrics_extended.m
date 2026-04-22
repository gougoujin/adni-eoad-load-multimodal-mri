% Repository usage summary
% Description: Extended diffusion extraction script including additional metrics beyond FA/MD.
% Usage: Optional extended diffusion pipeline.
% Outputs: Use only if the expanded DTI metric set is needed.
% Note: Update input paths, toolboxes, and filenames for your local environment.

% Public repository version: update file paths, toolboxes, and local settings before running.
% This script/function was lightly sanitized for sharing and may require project-specific inputs.

function extract_roi_diffusion_metrics_extended()
clear; clc;

proc_root = '';
out_dir   = '';
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

nROI = 116;

rows = {};
FA_all = [];
MD_all = [];
AD_all = [];
RD_all = [];
MO_all = [];

groups = dir(fullfile(proc_root, 'group*'));
groups = groups([groups.isdir]);

for g = 1:numel(groups)
    gname = groups(g).name;
    gdir  = fullfile(proc_root, gname);

    subs = dir(fullfile(gdir, 'sub*'));
    subs = subs([subs.isdir]);

    for s = 1:numel(subs)
        sid = subs(s).name;
        sdir = fullfile(gdir, sid);

        dti_dir = fullfile(sdir, 'dti');
        roi_dir = fullfile(sdir, 'roi_dwi');

        fa_file = fullfile(dti_dir, 'dti_FA.nii.gz');
        md_file = fullfile(dti_dir, 'dti_MD.nii.gz');
        l1_file = fullfile(dti_dir, 'dti_L1.nii.gz');
        l2_file = fullfile(dti_dir, 'dti_L2.nii.gz');
        l3_file = fullfile(dti_dir, 'dti_L3.nii.gz');
        mo_file = fullfile(dti_dir, 'dti_MO.nii.gz');

        need_files = {fa_file, md_file, l1_file, l2_file, l3_file, mo_file};
        if ~all(cellfun(@(x) exist(x,'file')>0, need_files)) || ~exist(roi_dir,'dir')
            fprintf('[SKIP] 缺少 dti/roi_dwi: %s/%s\n', gname, sid);
            continue;
        end

        try
            FA = double(niftiread(fa_file));
            MD = double(niftiread(md_file));
            L1 = double(niftiread(l1_file));
            L2 = double(niftiread(l2_file));
            L3 = double(niftiread(l3_file));
            MO = double(niftiread(mo_file));
            AD = L1;
            RD = (L2 + L3) / 2;
        catch ME
            fprintf('[SKIP] 读取失败 %s/%s: %s\n', gname, sid, ME.message);
            continue;
        end

        FA_roi = nan(1,nROI);
        MD_roi = nan(1,nROI);
        AD_roi = nan(1,nROI);
        RD_roi = nan(1,nROI);
        MO_roi = nan(1,nROI);

        for r = 1:nROI
            roi_file = fullfile(roi_dir, sprintf('roi_%d.nii.gz', r));
            if ~exist(roi_file,'file'), continue; end

            M = niftiread(roi_file);
            if ndims(M)==4, M = M(:,:,:,1); end
            mask = M > 0;

            if ~isequal(size(mask), size(FA)), continue; end

            FA_roi(r) = mean_valid(FA(mask), true);
            MD_roi(r) = mean_valid(MD(mask), true);
            AD_roi(r) = mean_valid(AD(mask), true);
            RD_roi(r) = mean_valid(RD(mask), true);
            MO_roi(r) = mean_valid(MO(mask), false); % MO 可为负，不限制 >0
        end

        rows(end+1,:) = {gname, sid}; %#ok<AGROW>
        FA_all(end+1,:) = FA_roi; %#ok<AGROW>
        MD_all(end+1,:) = MD_roi; %#ok<AGROW>
        AD_all(end+1,:) = AD_roi; %#ok<AGROW>
        RD_all(end+1,:) = RD_roi; %#ok<AGROW>
        MO_all(end+1,:) = MO_roi; %#ok<AGROW>

        fprintf('[OK] %s/%s\n', gname, sid);
    end
end

subject_info = rows; %#ok<NASGU>
save(fullfile(out_dir, 'dti_roi_all_subjects.mat'), ...
    'subject_info','FA_all','MD_all','AD_all','RD_all','MO_all');

fprintf('\n完成：%s\n', fullfile(out_dir, 'dti_roi_all_subjects.mat'));
end

function m = mean_valid(x, positive_only)
x = x(isfinite(x));
if positive_only
    x = x(x > 0);
end
if isempty(x)
    m = NaN;
else
    m = mean(x);
end
end
