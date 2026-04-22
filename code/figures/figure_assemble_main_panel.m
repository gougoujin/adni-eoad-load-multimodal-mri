% Repository usage summary
% Description: Assemble a composite main figure panel from previously exported images.
% Usage: Run after the component image panels have been created.
% Outputs: Exports a combined figure image for manuscript use.
% Note: Update input paths, toolboxes, and filenames for your local environment.

% Public repository version: update file paths, toolboxes, and local settings before running.
% This script/function was lightly sanitized for sharing and may require project-specific inputs.

function out_file = figure_assemble_main_panel(brain_img_file, panelB_img_file, out_file)
% Assemble Figure 1 main panel (fixed imwrite export).
%
% Example:
% figure_assemble_main_panel( ...
%   '', ...
%   '', ...
%   '');

    if nargin < 1 || isempty(brain_img_file)
        brain_img_file = '';
    end
    if nargin < 2 || isempty(panelB_img_file)
        panelB_img_file = '';
    end
    if nargin < 3 || isempty(out_file)
        out_file = '';
    end

    if ~isfile(brain_img_file)
        error(['缺少文件: Figure1_panelA_brainnet.png\n原路径应为: /path/to/project ' brain_img_file]);
    end
    if ~isfile(panelB_img_file)
        error(['缺少文件: Figure1_panelB_frontal_violin.png\n原路径应为: /path/to/project ' panelB_img_file]);
    end

    out_dir = fileparts(out_file);
    if ~isempty(out_dir) && ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end

    A = imread(brain_img_file);
    B = imread(panelB_img_file);

    if size(A,3) == 1, A = repmat(A, [1 1 3]); end
    if size(B,3) == 1, B = repmat(B, [1 1 3]); end

    target_h = max(size(A,1), size(B,1));
    A2 = imresize(A, [target_h NaN], 'bicubic');
    B2 = imresize(B, [target_h NaN], 'bicubic');

    gap = 40;
    pad = 40;
    canvas_h = target_h + 2*pad;
    canvas_w = size(A2,2) + size(B2,2) + gap + 2*pad;
    canvas = uint8(255 * ones(canvas_h, canvas_w, 3));

    yA = pad + 1; xA = pad + 1;
    canvas(yA:yA+size(A2,1)-1, xA:xA+size(A2,2)-1, :) = A2;

    yB = pad + 1; xB = xA + size(A2,2) + gap;
    canvas(yB:yB+size(B2,1)-1, xB:xB+size(B2,2)-1, :) = B2;

    % PNG export: do not pass numeric Resolution directly in this MATLAB path
    imwrite(canvas, out_file);

    [p,f,~] = fileparts(out_file);
    tif_file = fullfile(p, [f '.tif']);
    imwrite(canvas, tif_file, 'Compression', 'none');

    fprintf('Save/path/to/project', out_file, tif_file);
end
