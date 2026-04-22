% Repository usage summary
% Description: Compute ROI-level structural-functional coupling for all matched subjects.
% Usage: Run after subject matching to generate a 116 x N ROI coupling matrix.
% Outputs: Outputs regional coupling matrices and summary text files.
% Note: Update input paths, toolboxes, and filenames for your local environment.

% Public repository version: update file paths, toolboxes, and local settings before running.
% This script/function was lightly sanitized for sharing and may require project-specific inputs.

clear; clc;

%% ========= paths =========
match_file = '';
out_mat    = '';
out_txt    = '';

%% ========= load matched subjects =========
load(match_file, 'matched');

nROI = 116;
nSub = size(matched, 1);

roi_coupling_all = nan(nROI, nSub);
summary_cell = cell(nSub, 7);

for i = 1:nSub
    group_name = matched{i,1};
    sid        = matched{i,2};
    sc_file    = matched{i,3};
    fc_file    = matched{i,4};

    fprintf('\n============================\n');
    fprintf('Processing %s %s\n', group_name, sid);
    fprintf('============================\n');

    %% ----- load SC -----
    S = load(sc_file);

    if isfield(S, 'SC_roiMean_sym')
        SC = S.SC_roiMean_sym;
        sc_name = 'SC_roiMean_sym';
    elseif isfield(S, 'SC_sym_mean')
        SC = S.SC_sym_mean;
        sc_name = 'SC_sym_mean';
    elseif isfield(S, 'SC')
        SC = S.SC;
        sc_name = 'SC';
    else
        warning('%s %s: no SC matrix found', group_name, sid);
        continue;
    end

    %% ----- load FC -----
    [~,~,ext] = fileparts(fc_file);

    if strcmpi(ext, '.txt')
        FC = readmatrix(fc_file);
    elseif strcmpi(ext, '.mat')
        T = load(fc_file);
        fn = fieldnames(T);
        FC = [];
        for k = 1:numel(fn)
            X = T.(fn{k});
            if isnumeric(X) && isequal(size(X), [nROI nROI])
                FC = X;
                break;
            end
        end
        if isempty(FC)
            warning('%s %s: no 116x116 FC matrix found in mat', group_name, sid);
            continue;
        end
    else
        warning('%s %s: unsupported FC file type', group_name, sid);
        continue;
    end

    %% ----- sanity check -----
    if ~isequal(size(SC), [nROI nROI])
        warning('%s %s: SC size is not 116x116', group_name, sid);
        continue;
    end
    if ~isequal(size(FC), [nROI nROI])
        warning('%s %s: FC size is not 116x116', group_name, sid);
        continue;
    end

    %% ----- remove diagonal -----
    SC(1:nROI+1:end) = 0;
    FC(1:nROI+1:end) = 0;

    roi_coupling = nan(nROI,1);
    valid_edges_per_roi = zeros(nROI,1);

    %% ----- ROI-level coupling -----
    for r = 1:nROI
        x = SC(r, :);
        y = FC(r, :);

        keep = true(1, nROI);
        keep(r) = false;

        x = x(keep);
        y = y(keep);

        % only use nonzero SC edges with finite FC
        nz = (x > 0) & isfinite(x) & isfinite(y);

        valid_edges_per_roi(r) = nnz(nz);

        if nnz(nz) >= 3
            roi_coupling(r) = corr(log1p(x(nz))', y(nz)', 'type', 'Spearman');
        end
    end

    roi_coupling_all(:, i) = roi_coupling;

    summary_cell{i,1} = group_name;
    summary_cell{i,2} = sid;
    summary_cell{i,3} = mean(roi_coupling, 'omitnan');
    summary_cell{i,4} = median(roi_coupling, 'omitnan');
    summary_cell{i,5} = nnz(~isnan(roi_coupling));
    summary_cell{i,6} = mean(valid_edges_per_roi, 'omitnan');
    summary_cell{i,7} = sc_name;

    fprintf('%s %s done | mean ROI coupling = %.4f | valid ROI = %d\n', ...
        group_name, sid, summary_cell{i,3}, summary_cell{i,5});
end

%% ========= save MAT =========
save(out_mat, 'matched', 'roi_coupling_all', 'summary_cell');

%% ========= write summary TXT =========
fid = fopen(out_txt, 'w');
fprintf(fid, 'group\tsub\tmean_roi_coupling\tmedian_roi_coupling\tvalid_roi\tmean_valid_edges_per_roi\tsc_matrix_used\n');

for i = 1:nSub
    if isempty(summary_cell{i,1})
        continue;
    end
    fprintf(fid, '%s\t%s\t%.6f\t%.6f\t%d\t%.6f\t%s\n', ...
        summary_cell{i,1}, ...
        summary_cell{i,2}, ...
        summary_cell{i,3}, ...
        summary_cell{i,4}, ...
        summary_cell{i,5}, ...
        summary_cell{i,6}, ...
        summary_cell{i,7});
end
fclose(fid);

fprintf('\nOutput file/path/to/project', out_mat, out_txt);
