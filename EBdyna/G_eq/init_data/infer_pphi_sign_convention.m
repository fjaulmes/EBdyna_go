function m = infer_pphi_sign_convention(m)
%INFER_PPHI_SIGN_CONVENTION Infer the psi sign used in canonical P_phi.
%
% sign = +1:
%   BR = -(1/R) dpsi/dZ
%   BZ = +(1/R) dpsi/dR
%   R*Aphi = +psi
%   Pphi/e = (m/e)*R*vphi + Z*psi
%
% sign = -1:
%   BR = +(1/R) dpsi/dZ
%   BZ = -(1/R) dpsi/dR
%   R*Aphi = -psi
%   Pphi/e = (m/e)*R*vphi - Z*psi
%
% Required fields:
%   m.R_grid
%   m.Z_grid
%   m.psi_XZ
%   m.Br_XZ and m.Bz_XZ
%   or m.B_2D.BR and m.B_2D.BZ
%
% Added fields:
%   m.pphi_sign_conv
%   m.pphi_sign_check

%% Check required fields

required_fields = {'R_grid', 'Z_grid', 'psi_XZ'};

for i = 1:numel(required_fields)
    if ~isfield(m, required_fields{i})
        error('infer_pphi_sign_convention:MissingField', ...
            'Missing required field m.%s.', required_fields{i});
    end
end

R   = double(m.R_grid(:));
Z   = double(m.Z_grid(:));
psi = double(m.psi_XZ);

nR = numel(R);
nZ = numel(Z);

%% Get poloidal magnetic-field maps

if isfield(m, 'Br_XZ') && isfield(m, 'Bz_XZ')

    BR = double(m.Br_XZ);
    BZ = double(m.Bz_XZ);

elseif isfield(m, 'B_2D') ...
        && isfield(m.B_2D, 'BR') ...
        && isfield(m.B_2D, 'BZ')

    BR = double(m.B_2D.BR);
    BZ = double(m.B_2D.BZ);

else
    error('infer_pphi_sign_convention:MissingField', ...
        ['Cannot find poloidal-field maps. Expected either ' ...
         'm.Br_XZ/m.Bz_XZ or m.B_2D.BR/m.B_2D.BZ.']);
end

%% Standardize all maps to [nR x nZ]

[psi, psi_was_transposed] = orient_RZ_map(psi, nR, nZ, 'psi_XZ');
[BR,  BR_was_transposed]  = orient_RZ_map(BR,  nR, nZ, 'BR');
[BZ,  BZ_was_transposed]  = orient_RZ_map(BZ,  nR, nZ, 'BZ');

%% Calculate flux derivatives

% psi has dimensions [R,Z].
%
% gradient's first output differentiates along columns, hence Z.
% gradient's second output differentiates along rows, hence R.

[dpsi_dZ, dpsi_dR] = gradient(psi, Z, R);

RR = repmat(R, 1, nZ);

%% Reconstruct B for both conventions

% Convention +1: R*Aphi = +psi
BR_plus = -dpsi_dZ ./ RR;
BZ_plus = +dpsi_dR ./ RR;

% Convention -1: R*Aphi = -psi
BR_minus = +dpsi_dZ ./ RR;
BZ_minus = -dpsi_dR ./ RR;

%% Build valid comparison mask

mask = isfinite(psi) ...
     & isfinite(BR) ...
     & isfinite(BZ) ...
     & isfinite(BR_plus) ...
     & isfinite(BZ_plus);

% Exclude boundaries because gradient uses one-sided differences there.
if nR > 2
    mask([1 end], :) = false;
end

if nZ > 2
    mask(:, [1 end]) = false;
end

% Prefer the plasma region if normalized flux is available.
if isfield(m, 'psi_n_XZ')

    psi_n = double(m.psi_n_XZ);
    psi_n = orient_RZ_map(psi_n, nR, nZ, 'psi_n_XZ');

    plasma_mask = isfinite(psi_n) ...
                & psi_n >= 0 ...
                & psi_n <= 1;

    if nnz(mask & plasma_mask) >= 100
        mask = mask & plasma_mask;
        log_msg('debug', ...
            'Using %d plasma-region points for Pphi sign inference.', ...
            nnz(mask));
    else
        log_msg('warn', ...
            ['Normalized-flux plasma mask contains too few points. ' ...
             'Using the full valid map region instead.']);
    end
end

if nnz(mask) < 10
    error('infer_pphi_sign_convention:InsufficientData', ...
        'Too few valid points are available for Pphi sign inference.');
end

%% Compare both conventions

B_file  = [BR(mask);       BZ(mask)];
B_plus  = [BR_plus(mask);  BZ_plus(mask)];
B_minus = [BR_minus(mask); BZ_minus(mask)];

field_rms = sqrt(mean(B_file.^2));

if field_rms <= eps
    error('infer_pphi_sign_convention:ZeroField', ...
        'The poloidal magnetic field is too small for sign inference.');
end

error_plus = sqrt(mean((B_plus-B_file).^2)) / field_rms;
error_minus = sqrt(mean((B_minus-B_file).^2)) / field_rms;

correlation_plus = local_correlation(B_plus, B_file);
correlation_minus = local_correlation(B_minus, B_file);

%% Select sign convention

if error_plus < error_minus
    m.pphi_sign_conv = +1;
else
    m.pphi_sign_conv = -1;
end

best_error = min(error_plus, error_minus);
worst_error = max(error_plus, error_minus);
error_ratio = worst_error / max(best_error, eps);

if error_ratio < 10
    log_msg('warn', ...
        ['Pphi sign inference is weakly separated: ' ...
         'error ratio = %.3g. Check map orientation and normalization.'], ...
        error_ratio);
end

%% Store diagnostics

m.pphi_sign_check = struct( ...
    'error_plus', error_plus, ...
    'error_minus', error_minus, ...
    'correlation_plus', correlation_plus, ...
    'correlation_minus', correlation_minus, ...
    'error_ratio', error_ratio, ...
    'number_of_points', nnz(mask), ...
    'psi_was_transposed', psi_was_transposed, ...
    'BR_was_transposed', BR_was_transposed, ...
    'BZ_was_transposed', BZ_was_transposed);

%% Log result

log_msg('debug', ...
    ['Pphi sign check: error(+1) = %.6e, error(-1) = %.6e, ' ...
     'corr(+1) = %.9f, corr(-1) = %.9f.'], ...
    error_plus, error_minus, ...
    correlation_plus, correlation_minus);

if m.pphi_sign_conv == +1

    log_msg('info', ...
        ['Pphi convention inferred from 2D maps: sign = +1. ' ...
         'R*Aphi = +psi and Pphi/e = (m/e)*R*vphi + Z*psi.']);

else

    log_msg('info', ...
        ['Pphi convention inferred from 2D maps: sign = -1. ' ...
         'R*Aphi = -psi and Pphi/e = (m/e)*R*vphi - Z*psi.']);

end

end


function [A, was_transposed] = orient_RZ_map(A, nR, nZ, field_name)
%ORIENT_RZ_MAP Return map with dimensions [nR x nZ].

was_transposed = false;

if isequal(size(A), [nR, nZ])

    return

elseif isequal(size(A), [nZ, nR])

    A = A.';
    was_transposed = true;

else

    error('infer_pphi_sign_convention:UnexpectedDimensions', ...
        '%s has size %s; expected [%d %d] or [%d %d].', ...
        field_name, mat2str(size(A)), ...
        nR, nZ, nZ, nR);

end

end


function c = local_correlation(a, b)
%LOCAL_CORRELATION Calculate correlation without toolbox-specific options.

a = a(:);
b = b(:);

a = a - mean(a);
b = b - mean(b);

denominator = sqrt(sum(a.^2) * sum(b.^2));

if denominator <= eps
    c = NaN;
else
    c = sum(a .* b) / denominator;
end

end