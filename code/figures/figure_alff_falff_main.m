% Repository usage summary
% Description: Create the main ALFF/fALFF manuscript figure from ROI-level values.
% Usage: Run after functional metric extraction and ROI selection are complete.
% Outputs: Exports figure images and value tables used by downstream plotting scripts.
% Note: Update input paths, toolboxes, and filenames for your local environment.

% Public repository version: update file paths, toolboxes, and local settings before running.
% This script/function was lightly sanitized for sharing and may require project-specific inputs.

function out_tif = figure_alff_falff_main(cov_file, aal_label_file, aal_atlas_file, dpabi_group_roots, out_dir)
% Figure 6: ALFF / fALFF main figure
% v5: fix outerjoin-renamed variables (group_num/age/sex/sub)

    if nargin < 1 || isempty(cov_file)
        cov_file = '';
    end
    if nargin < 2 || isempty(aal_label_file)
        aal_label_file = '';
    end
    if nargin < 3 || isempty(aal_atlas_file)
        aal_atlas_file = '';
    end
    if nargin < 4 || isempty(dpabi_group_roots)
        dpabi_group_roots = { ...
            '', ...
            '', ...
            '' ...
            };
    end
    if nargin < 5 || isempty(out_dir)
        out_dir = '';
    end

    if ~isfile(cov_file)
        error(['缺少文件: biao.xlsx\n原路径应为: /path/to/project ' cov_file]);
    end
    if ~isfile(aal_label_file)
        error(['缺少文件: AAL116_Labels.txt\n原路径应为: /path/to/project ' aal_label_file]);
    end
    if ~isfile(aal_atlas_file)
        error(['缺少文件: raal.nii\n原路径应为: /path/to/project ' aal_atlas_file]);
    end
    for g = 1:numel(dpabi_group_roots)
        if ~exist(dpabi_group_roots{g}, 'dir')
            error(['缺少目录: DPABI group root\n原路径应为: ' dpabi_group_roots{g}]);
        end
    end
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end

    alff_roi_labels = {'Hippocampus_R','Hippocampus_L','Thalamus_R','Caudate_R','Pallidum_L','Insula_L'};
    falff_roi_labels = {'Cuneus_L','Cerebelum_Crus1_L','Occipital_Inf_R','Temporal_Mid_L','Temporal_Pole_Sup_R','Precuneus_L'};

    aal_labels = local_read_labels(aal_label_file, 116);
    atlas = double(niftiread(aal_atlas_file));

    alff_roi_idx = local_labels_to_idx(alff_roi_labels, aal_labels);
    falff_roi_idx = local_labels_to_idx(falff_roi_labels, aal_labels);

    keep_alff = ~isnan(alff_roi_idx);
    keep_falff = ~isnan(falff_roi_idx);
    alff_roi_idx = alff_roi_idx(keep_alff);
    falff_roi_idx = falff_roi_idx(keep_falff);

    if isempty(alff_roi_idx) && isempty(falff_roi_idx)
        error('AAL116_Labels.txt 中一个 ROI 都没匹配上，请检查标签命名。');
    end

    raw = readcell(cov_file);
    raw = raw(2:end,:);
    n = size(raw,1);

    T = table;
    T.group_num = nan(n,1);
    T.sub = strings(n,1);
    T.age = nan(n,1);
    T.sex = nan(n,1);

    for i = 1:n
        T.group_num(i) = local_to_num(raw{i,1});
        s = raw{i,3};
        if isnumeric(s) && ~isnan(s)
            T.sub(i) = sprintf('sub%03d', round(s));
        else
            tok = regexp(char(string(s)), '(\d+)', 'tokens', 'once');
            if ~isempty(tok)
                T.sub(i) = sprintf('sub%03d', str2double(tok{1}));
            end
        end
        T.age(i) = local_to_num(raw{i,4});
        T.sex(i) = local_to_num(raw{i,5});
    end
    T.tag = string(T.group_num) + "_" + T.sub;
    T = T(isfinite(T.group_num) & strlength(T.sub)>0 & isfinite(T.age) & isfinite(T.sex), :);

    Xall = table();
    fprintf('\n========== 提取 zALFF / zfALFF ==========\n');

    for g = 1:3
        group_root = dpabi_group_roots{g};
        fprintf('扫描 group%d: %s\n', g, group_root);

        subs = dir(fullfile(group_root, 'sub*'));
        subs = subs([subs.isdir]);

        if isempty(subs)
            alff_dir = fullfile(group_root, 'Results', 'ALFF_FunImgARCW');
            ff = dir(fullfile(alff_dir, 'zALFFMap_sub*.nii*'));
            sid_list = strings(numel(ff),1);
            for ii = 1:numel(ff)
                tok = regexp(ff(ii).name, 'zALFFMap_(sub\d+)\.nii', 'tokens', 'once');
                if isempty(tok)
                    tok = regexp(ff(ii).name, 'zALFFMap_(sub\d+)\.nii\.gz', 'tokens', 'once');
                end
                if ~isempty(tok), sid_list(ii) = string(tok{1}); end
            end
            sid_list = unique(sid_list(strlength(sid_list)>0), 'stable');
            subs = struct('name', cellstr(sid_list), 'folder', repmat({group_root}, numel(sid_list),1), 'date', [], 'bytes', [], 'isdir', true, 'datenum', []);
        end

        rows = {};
        ALFF_vals = [];
        fALFF_vals = [];

        for s = 1:numel(subs)
            sid = subs(s).name;
            alff_file  = fullfile(group_root, 'Results', 'ALFF_FunImgARCW',  sprintf('zALFFMap_%s.nii', sid));
            falff_file = fullfile(group_root, 'Results', 'fALFF_FunImgARCW', sprintf('zfALFFMap_%s.nii', sid));

            if ~exist(alff_file, 'file')
                gz = [alff_file '.gz'];
                if exist(gz, 'file'), alff_file = gz; else, alff_file = ''; end
            end
            if ~exist(falff_file, 'file')
                gz = [falff_file '.gz'];
                if exist(gz, 'file'), falff_file = gz; else, falff_file = ''; end
            end
            if isempty(alff_file) && isempty(falff_file)
                continue;
            end

            arow = nan(1, numel(alff_roi_idx));
            frow = nan(1, numel(falff_roi_idx));

            if ~isempty(alff_file)
                A = double(local_read_nifti_auto(alff_file));
                if isequal(size(A), size(atlas))
                    for r = 1:numel(alff_roi_idx)
                        roi = alff_roi_idx(r);
                        mask = atlas == roi;
                        arow(r) = mean(A(mask & isfinite(A)), 'omitnan');
                    end
                end
            end

            if ~isempty(falff_file)
                F = double(local_read_nifti_auto(falff_file));
                if isequal(size(F), size(atlas))
                    for r = 1:numel(falff_roi_idx)
                        roi = falff_roi_idx(r);
                        mask = atlas == roi;
                        frow(r) = mean(F(mask & isfinite(F)), 'omitnan');
                    end
                end
            end

            rows(end+1,:) = {g, sid, sprintf('%d_%s', g, sid)}; %#ok<AGROW>
            ALFF_vals(end+1,:) = arow; %#ok<AGROW>
            fALFF_vals(end+1,:) = frow; %#ok<AGROW>
        end

        if isempty(rows), continue; end

        X = table;
        X.group_num = cell2mat(rows(:,1));
        X.sub = string(rows(:,2));
        X.tag = string(rows(:,3));

        for r = 1:numel(alff_roi_idx)
            roi = alff_roi_idx(r);
            vname = matlab.lang.makeValidName(['zALFF_' aal_labels{roi}]);
            X.(vname) = ALFF_vals(:,r);
        end

        for r = 1:numel(falff_roi_idx)
            roi = falff_roi_idx(r);
            vname = matlab.lang.makeValidName(['zfALFF_' aal_labels{roi}]);
            X.(vname) = fALFF_vals(:,r);
        end

        if isempty(Xall)
            Xall = X;
        else
            Xall = [Xall; X]; %#ok<AGROW>
        end
    end

    if isempty(Xall)
        error(['未提取到任何 ALFF/fALFF 数据。\n' ...
               '已按原脚本方式查找:\n' ...
               '  <group_root>\Results\ALFF_FunImgARCW\zALFFMap_subXXX.nii(.gz)\n' ...
               '  <group_root>\Results\fALFF_FunImgARCW\zfALFFMap_subXXX.nii(.gz)\n']);
    end

    Tm = outerjoin(T, Xall, 'Keys', 'tag', 'MergeKeys', true);
    Tm = local_fix_join_names(Tm);
    writetable(Tm, fullfile(out_dir, 'Figure6_alff_falff_values.csv'));

    alff_vars = cell(numel(alff_roi_idx),1);
    alff_titles = cell(numel(alff_roi_idx),1);
    for r = 1:numel(alff_roi_idx)
        roi = alff_roi_idx(r);
        alff_vars{r} = matlab.lang.makeValidName(['zALFF_' aal_labels{roi}]);
        alff_titles{r} = strrep(aal_labels{roi}, '_', ' ');
    end

    falff_vars = cell(numel(falff_roi_idx),1);
    falff_titles = cell(numel(falff_roi_idx),1);
    for r = 1:numel(falff_roi_idx)
        roi = falff_roi_idx(r);
        falff_vars{r} = matlab.lang.makeValidName(['zfALFF_' aal_labels{roi}]);
        falff_titles{r} = strrep(aal_labels{roi}, '_', ' ');
    end

    n1 = numel(alff_vars);
    n2 = numel(falff_vars);
    ncol = max([n1 n2 1]);

    fig = figure('Color','w', 'Position', [60 60 320*ncol 950]);
    tl = tiledlayout(fig, 2, ncol, 'TileSpacing','compact', 'Padding','compact');

    c1 = [0.85 0.40 0.40];
    c2 = [0.35 0.60 0.85];
    c3 = [0.45 0.75 0.50];
    letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';

    for k = 1:ncol
        ax = nexttile(tl, k); hold(ax, 'on');
        if k <= n1
            local_plot_one_panel(ax, Tm, alff_vars{k}, alff_titles{k}, c1, c2, c3, letters(k), 'zALFF');
        else
            axis(ax, 'off');
        end
    end

    for k = 1:ncol
        ax = nexttile(tl, ncol+k); hold(ax, 'on');
        if k <= n2
            local_plot_one_panel(ax, Tm, falff_vars{k}, falff_titles{k}, c1, c2, c3, letters(ncol+k), 'zfALFF');
        else
            axis(ax, 'off');
        end
    end

    annotation(fig, 'textbox', [0.18 0.965 0.25 0.03], 'String', 'ALFF representative ROIs', ...
        'EdgeColor','none', 'FontName','Arial', 'FontSize',12, 'FontWeight','bold', 'HorizontalAlignment','center');
    annotation(fig, 'textbox', [0.60 0.965 0.25 0.03], 'String', 'fALFF representative ROIs', ...
        'EdgeColor','none', 'FontName','Arial', 'FontSize',12, 'FontWeight','bold', 'HorizontalAlignment','center');

    out_tif = fullfile(out_dir, 'Figure6_alff_falff_main.tif');
    out_png = fullfile(out_dir, 'Figure6_alff_falff_main.png');
    exportgraphics(fig, out_tif, 'Resolution', 600, 'BackgroundColor','white', 'ContentType','image');
    exportgraphics(fig, out_png, 'Resolution', 300, 'BackgroundColor','white', 'ContentType','image');

    fprintf('Save/path/to/project', out_tif, out_png);
end

function Tm = local_fix_join_names(Tm)
    vn = string(Tm.Properties.VariableNames);

    if ~ismember("group_num", vn)
        cands = vn(contains(vn, "group_num"));
        if ~isempty(cands)
            Tm.group_num = Tm.(cands(1));
        end
    end
    if ~ismember("age", vn)
        cands = vn(strcmpi(vn, "age_T") | strcmpi(vn, "age_Xall") | contains(vn, "age"));
        if ~isempty(cands)
            Tm.age = Tm.(cands(1));
        end
    end
    if ~ismember("sex", vn)
        cands = vn(strcmpi(vn, "sex_T") | strcmpi(vn, "sex_Xall") | contains(vn, "sex"));
        if ~isempty(cands)
            Tm.sex = Tm.(cands(1));
        end
    end
    if ~ismember("sub", vn)
        cands = vn(contains(vn, "sub"));
        if ~isempty(cands)
            Tm.sub = string(Tm.(cands(1)));
        end
    end
end

function local_plot_one_panel(ax, T, varname, ttl, c1, c2, c3, panel_letter, ylab)
    if ~ismember(varname, T.Properties.VariableNames)
        text(ax, 0.5, 0.5, 'Missing variable', 'Units','normalized', 'HorizontalAlignment','center', 'FontName','Arial');
        axis(ax, 'off');
        return;
    end
    if ~ismember('group_num', T.Properties.VariableNames) || ~ismember('age', T.Properties.VariableNames) || ~ismember('sex', T.Properties.VariableNames)
        text(ax, 0.5, 0.5, 'Missing covariates after join', 'Units','normalized', 'HorizontalAlignment','center', 'FontName','Arial');
        axis(ax, 'off');
        return;
    end

    y = double(T.(varname));
    g = double(T.group_num);

    idx1 = g==1 & isfinite(y) & isfinite(T.age) & isfinite(T.sex);
    idx2 = g==2 & isfinite(y) & isfinite(T.age) & isfinite(T.sex);
    idx3 = g==3 & isfinite(y) & isfinite(T.age) & isfinite(T.sex);

    local_scatter_box(ax, 1, y(idx1), c1);
    local_scatter_box(ax, 2, y(idx2), c2);
    local_scatter_box(ax, 3, y(idx3), c3);

    ok = isfinite(y) & isfinite(g) & isfinite(T.age) & isfinite(T.sex);
    if sum(ok) >= 12 && numel(unique(g(ok))) >= 2
        D = table(y(ok), categorical(g(ok)), T.age(ok), T.sex(ok), 'VariableNames', {'y','group','age','sex'});
        lm = fitlm(D, 'y ~ group + age + sex');
        a = anova(lm, 'summary');
        p_model = NaN;
        if ismember('Term', a.Properties.VariableNames) && ismember('pValue', a.Properties.VariableNames)
            idxm = strcmp(string(a.Term), 'Model');
            if any(idxm), p_model = a.pValue(find(idxm,1)); end
        end
        text(ax, 0.03, 0.93, sprintf('%c  p=%s', panel_letter, local_p(p_model)), ...
            'Units','normalized', 'FontName','Arial', 'FontSize',9, 'FontWeight','bold', ...
            'HorizontalAlignment','left', 'VerticalAlignment','top');
    else
        text(ax, 0.03, 0.93, sprintf('%c', panel_letter), ...
            'Units','normalized', 'FontName','Arial', 'FontSize',9, 'FontWeight','bold', ...
            'HorizontalAlignment','left', 'VerticalAlignment','top');
    end

    set(ax, 'XLim', [0.5 3.5], 'XTick', [1 2 3], 'XTickLabel', {'EOAD','LOAD','HC'}, ...
        'FontName','Arial', 'FontSize',9, 'Box','off', 'TickDir','out');
    title(ax, ttl, 'FontName','Arial', 'FontSize',10, 'FontWeight','normal', 'Interpreter','none');
    ylabel(ax, ylab, 'FontName','Arial', 'FontSize',9);
end

function local_scatter_box(ax, xpos, y, c)
    if isempty(y), return; end
    boxchart(ax, ones(size(y))*xpos, y, 'BoxFaceColor', c, 'WhiskerLineColor', c*0.7, ...
        'MarkerStyle', 'none', 'BoxFaceAlpha', 0.55, 'LineWidth', 1.0);
    scatter(ax, xpos + (rand(size(y))-0.5)*0.20, y, 12, c, 'filled', 'MarkerFaceAlpha', 0.60, 'MarkerEdgeAlpha', 0.15);
end

function labels = local_read_labels(label_file, nROI)
    labels = cell(nROI,1);
    fid = fopen(label_file, 'r');
    if fid < 0, error('无法打开标签文件: %s', label_file); end
    while ~feof(fid)
        line = strtrim(fgetl(fid));
        if isempty(line), continue; end
        parts = regexp(line, '\s+', 'split');
        if numel(parts) < 2, continue; end
        idx = str2double(parts{1});
        if isfinite(idx) && idx >= 1 && idx <= nROI
            labels{idx} = strtrim(parts{2});
        end
    end
    fclose(fid);
    for i = 1:nROI
        if isempty(labels{i}), labels{i} = sprintf('ROI_%d', i); end
    end
end

function idx = local_labels_to_idx(target_labels, aal_labels)
    idx = nan(numel(target_labels),1);
    aal = string(aal_labels);
    for i = 1:numel(target_labels)
        m = find(aal == string(target_labels{i}), 1);
        if ~isempty(m), idx(i) = m; end
    end
end

function x = local_to_num(v)
    if isnumeric(v), x = double(v); else, x = str2double(string(v)); end
end

function V = local_read_nifti_auto(f)
    if endsWith(lower(f), '.gz')
        tmpdir = tempname; mkdir(tmpdir); gunzip(f, tmpdir);
        dd = dir(fullfile(tmpdir, '*.nii'));
        if isempty(dd), error('无法解压 NIfTI: %s', f); end
        nii = fullfile(tmpdir, dd(1).name);
        V = niftiread(nii);
        delete(nii); rmdir(tmpdir);
    else
        V = niftiread(f);
    end
end

function s = local_p(p)
    if ~isfinite(p), s = 'NA';
    elseif p < 0.001, s = '<0.001';
    else, s = sprintf('%.3f', p);
    end
end
