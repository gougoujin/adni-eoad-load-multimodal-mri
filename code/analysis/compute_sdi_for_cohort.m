% Repository usage summary
% Description: Compute structural decoupling index for all matched subjects.
% Usage: Run after subject matching and ROI time-series preparation.
% Outputs: Writes cohort-level SDI MAT/TXT outputs for downstream statistics.
% Note: Update input paths, toolboxes, and filenames for your local environment.

% Public repository version: update file paths, toolboxes, and local settings before running.
% This script/function was lightly sanitized for sharing and may require project-specific inputs.

clear; clc;

match_file = '';

ts_roots = {
    'group1', '';
    'group2', '';
    'group3', '';
};

out_mat = '';
out_txt = '';

nROI = 116;

load(match_file, 'matched');

nSub = size(matched,1);
SDI_all = nan(nROI, nSub);
summary_cell = cell(nSub, 5);

for i = 1:nSub
    group_name = matched{i,1};
    sid        = matched{i,2};
    sc_file    = matched{i,3};

    fprintf('\n============================\n');
    fprintf('Processing %s %s\n', group_name, sid);
    fprintf('============================\n');

    idxg = find(strcmp(ts_roots(:,1), group_name), 1);
    if isempty(idxg)
        warning('No TS root for %s', group_name);
        continue;
    end
    ts_root = ts_roots{idxg,2};

    S = load(sc_file);
    if isfield(S, 'SC_roiMean_sym')
        SC = S.SC_roiMean_sym;
    elseif isfield(S, 'SC_sym_mean')
        SC = S.SC_sym_mean;
    elseif isfield(S, 'SC')
        SC = S.SC;
    else
        warning('%s %s: no SC found', group_name, sid);
        continue;
    end

    mat_file = fullfile(ts_root, ['ROISignals_' char(sid) '.mat']);
    txt_file = fullfile(ts_root, ['ROISignals_' char(sid) '.txt']);

    if exist(mat_file, 'file')
        ts_file = mat_file;
    elseif exist(txt_file, 'file')
        ts_file = txt_file;
    else
        warning('%s %s: no ROISignals file found', group_name, sid);
        continue;
    end

    [~,~,ext] = fileparts(ts_file);
    X = [];

    if strcmpi(ext, '.txt')
        X = readmatrix(ts_file);
    elseif strcmpi(ext, '.mat')
        T = load(ts_file);
        fn = fieldnames(T);
        for k = 1:numel(fn)
            Y = T.(fn{k});
            if isnumeric(Y)
                sz = size(Y);
                if numel(sz) == 2 && (sz(2) == nROI || sz(1) == nROI)
                    X = Y;
                    break;
                end
            end
        end
    end

    if isempty(X)
        warning('%s %s: cannot read ROI timeseries from %s', group_name, sid, ts_file);
        continue;
    end

    if size(X,1) == nROI && size(X,2) ~= nROI
        X = X';
    end

    if size(X,2) ~= nROI
        warning('%s %s: ROI timeseries is not T x 116', group_name, sid);
        continue;
    end

    [SDI_roi, out] = compute_sdi_for_subject(SC, X); %#ok<NASGU>

    SDI_all(:, i) = SDI_roi;

    summary_cell{i,1} = group_name;
    summary_cell{i,2} = sid;
    summary_cell{i,3} = mean(SDI_roi, 'omitnan');
    summary_cell{i,4} = median(SDI_roi, 'omitnan');
    summary_cell{i,5} = nnz(~isnan(SDI_roi));

    fprintf('%s %s done | mean SDI = %.4f | valid ROI = %d\n', ...
        group_name, sid, summary_cell{i,3}, summary_cell{i,5});
end

save(out_mat, 'matched', 'SDI_all', 'summary_cell');

fid = fopen(out_txt, 'w');
fprintf(fid, 'group\tsub\tmean_SDI\tmedian_SDI\tvalid_roi\n');
for i = 1:nSub
    if isempty(summary_cell{i,1}), continue; end
    fprintf(fid, '%s\t%s\t%.6f\t%.6f\t%d\n', ...
        summary_cell{i,1}, summary_cell{i,2}, ...
        summary_cell{i,3}, summary_cell{i,4}, summary_cell{i,5});
end
fclose(fid);

fprintf('\nOutput file/path/to/project', out_mat, out_txt);
