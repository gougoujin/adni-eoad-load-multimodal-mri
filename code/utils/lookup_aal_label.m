% Repository usage summary
% Description: Look up an AAL ROI label from a peak coordinate.
% Usage: [aal_id, aal_name] = lookup_aal_label(...)
% Outputs: Returns the numeric AAL label and human-readable region name.
% Note: Update input paths, toolboxes, and filenames for your local environment.

% Public repository version: update file paths, toolboxes, and local settings before running.
% This script/function was lightly sanitized for sharing and may require project-specific inputs.

function [aal_id, aal_name] = lookup_aal_label(peak_mm, Vaal, AALvol, dimX, dimY, dimZ, AAL_names)
% 根据峰值 MNI 坐标，在 AAL 掩模里查标签 ID 和名称
%
% 输入：
%   peak_mm   - 3×1 或 4×1 MNI 坐标（mm）
%   Vaal      - spm_vol 读入的 AAL 掩模头
%   AALvol    - AAL 掩模体数据（与 Vaal 对应）
%   dimX/Y/Z  - 掩模维度大小
%   AAL_names - containers.Map，键是 AAL_ID，值是名称
%
% 输出：
%   aal_id    - AAL 标签数字（0 表示无效）
%   aal_name  - 标签名称（Unknown 或 AAL_xx）

    % MNI -> AAL 体素坐标
    coord_aal = inv(Vaal.mat) * [peak_mm(1); peak_mm(2); peak_mm(3); 1];
    xA = round(coord_aal(1));
    yA = round(coord_aal(2));
    zA = round(coord_aal(3));

    aal_id   = 0;
    aal_name = 'Unknown';

    % 判断是否在体素范围内
    if xA>=1 && xA<=dimX && yA>=1 && yA<=dimY && zA>=1 && zA<=dimZ
        lbl = round(AALvol(xA,yA,zA));
        if lbl > 0
            aal_id = lbl;
            if ~isempty(AAL_names) && isKey(AAL_names, double(lbl))
                aal_name = AAL_names(double(lbl));
            else
                aal_name = sprintf('AAL_%d', lbl);
            end
        end
    end
end
