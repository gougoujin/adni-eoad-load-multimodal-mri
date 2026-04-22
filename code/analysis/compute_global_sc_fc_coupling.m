% Repository usage summary
% Description: Compute subject-level mean structural-functional coupling from matched SC and FC matrices.
% Usage: Run after subject matching to generate a cohort summary MAT/TXT for global coupling.
% Outputs: Outputs per-subject mean coupling and optional ROI-level intermediate values.
% Note: Update input paths, toolboxes, and filenames for your local environment.

% Public repository version: update file paths, toolboxes, and local settings before running.
% This script/function was lightly sanitized for sharing and may require project-specific inputs.

clear; clc;

load('', 'matched');

out_txt = '';
out_mat = '';

results = cell(size(matched,1), 6);
all_roi_coupling = [];

for i = 1:size(matched,1)
    group_name = matched{i,1};
    sid        = matched{i,2};
    sc_file    = matched{i,3};
    fc_file    = matched{i,4};

    % ----- 读 SC -----
    S = load(sc_file);
    if isfield(S, 'SC_roiMean_sym')
        SC = S.SC_roiMean_sym;
    elseif isfield(S, 'SC_sym_mean')
        SC = S.SC_sym_mean;
    elseif isfield(S, 'SC')
        SC = S.SC;
    else
        warning('%s %s: SC 文件缺少矩阵变量', group_name, sid);
        continue;
    end

    % ----- 读 FC -----
    [~,~,ext] = fileparts(fc_file);
    if strcmpi(ext, '.txt')
        FC = readmatrix(fc_file);
    elseif strcmpi(ext, '.mat')
        T = load(fc_file);
        fn = fieldnames(T);
        FC = [];
        for k = 1:numel(fn)
            X = T.(fn{k});
            if isnumeric(X) && isequal(size(X), [116 116])
                FC = X;
                break;
            end
        end
        if isempty(FC)
            warning('%s %s: mat 中没找到 116x116 FC', group_name, sid);
            continue;
        end
    else
        warning('%s %s: 不支持的 FC 文件类型', group_name, sid);
        continue;
    end

    if ~isequal(size(SC), [116 116]) || ~isequal(size(FC), [116 116])
        warning('%s %s: SC/FC 尺寸不是 116x116', group_name, sid);
        continue;
    end

    % 去对角线
    SC(1:117:end) = 0;
    FC(1:117:end) = 0;

    roi_coupling = nan(116,1);

    % 每个 ROI：该 ROI 到其它 115 个 ROI 的 SC 和 FC 做相关
    for r = 1:116
        x = SC(r, :);
        y = FC(r, :);

        keep = true(1,116);
        keep(r) = false;

        x = x(keep);
        y = y(keep);

        % 建议只在 SC 非零边上做 coupling
        nz = x > 0 & isfinite(x) & isfinite(y);

        if nnz(nz) >= 3
            roi_coupling(r) = corr(log1p(x(nz))', y(nz)', 'type', 'Spearman');
        end
    end

    subj_mean = mean(roi_coupling, 'omitnan');

    results{i,1} = group_name;
    results{i,2} = sid;
    results{i,3} = subj_mean;
    results{i,4} = nnz(~isnan(roi_coupling));
    results{i,5} = sc_file;
    results{i,6} = fc_file;

    all_roi_coupling(:, i) = roi_coupling; %#ok<SAGROW>

    fprintf('%s %s 完成 | mean coupling = %.4f | valid ROI = %d\n', ...
        group_name, sid, subj_mean, nnz(~isnan(roi_coupling)));
end

% 导出 summary
fid = fopen(out_txt, 'w');
fprintf(fid, 'group\tsub\tmean_coupling\tvalid_roi\tsc_file\tfc_file\n');
for i = 1:size(results,1)
    if isempty(results{i,1}), continue; end
    fprintf(fid, '%s\t%s\t%.6f\t%d\t%s\t%s\n', ...
        results{i,1}, results{i,2}, results{i,3}, results{i,4}, results{i,5}, results{i,6});
end
fclose(fid);

save(out_mat, 'results', 'all_roi_coupling');

fprintf('\n输出：\n%s\n%s\n', out_txt, out_mat);
