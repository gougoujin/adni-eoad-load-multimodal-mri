% Repository usage summary
% Description: Compute the structural decoupling index for one subject from SC and ROI time series.
% Usage: [SDI_roi,out] = compute_sdi_for_subject(SC, X)
% Outputs: Returns ROI-wise SDI and optional spectral decomposition outputs.
% Note: Update input paths, toolboxes, and filenames for your local environment.

% Public repository version: update file paths, toolboxes, and local settings before running.
% This script/function was lightly sanitized for sharing and may require project-specific inputs.

function [SDI_roi, out] = compute_sdi_for_subject(SC, X)
% SC: 116x116 structural connectivity
% X : T x 116 ROI time series

nROI = size(SC,1);

if size(SC,2) ~= nROI
    error('SC must be square.');
end

if size(X,2) ~= nROI
    error('X must be T x %d.', nROI);
end

SC = (SC + SC.') / 2;
SC(1:nROI+1:end) = 0;
SC(~isfinite(SC)) = 0;
SC(SC < 0) = 0;

deg = sum(SC, 2);
Dinv2 = diag(1 ./ sqrt(deg + eps));
L = eye(nROI) - Dinv2 * SC * Dinv2;

[U, Lambda] = eig(L);
[lambda, idx] = sort(diag(Lambda), 'ascend');
U = U(:, idx);

X = double(X);
X = X - mean(X, 1, 'omitnan');
X(isnan(X)) = 0;

Xg = X * U;

kcut = floor(nROI / 2);
Xg_low  = Xg(:, 1:kcut);
Xg_high = Xg(:, kcut+1:end);

X_low  = Xg_low  * U(:, 1:kcut)';
X_high = Xg_high * U(:, kcut+1:end)';

E_low  = sqrt(sum(X_low.^2, 1))';
E_high = sqrt(sum(X_high.^2, 1))';

SDI_roi = log2((E_high + eps) ./ (E_low + eps));

out = struct();
out.U = U;
out.lambda = lambda;
out.X_low = X_low;
out.X_high = X_high;
out.E_low = E_low;
out.E_high = E_high;
out.kcut = kcut;
end
