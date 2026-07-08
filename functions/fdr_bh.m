function [sig, pcrit] = fdr_bh(pvals, q)
%FDR_BH Benjamini-Hochberg false-discovery-rate correction.
%
%   [SIG, PCRIT] = FDR_BH(PVALS, Q) returns a logical array SIG indicating
%   which p-values survive Benjamini-Hochberg FDR correction at level Q.
%   Non-finite p-values are ignored. PCRIT is the largest p-value declared
%   significant, or NaN if no p-values survive.
%
%   This compact implementation is included to avoid adding another toolbox
%   dependency for a standard multiple-comparison correction.

if nargin < 2 || isempty(q)
    q = 0.05;
end

original_size = size(pvals);
p = pvals(:);
valid = isfinite(p);
sig_vector = false(size(p));
pcrit = NaN;

if ~any(valid)
    sig = reshape(sig_vector, original_size);
    return;
end

p_valid = p(valid);
p_sorted = sort(p_valid);
m = numel(p_sorted);
thresholds = (1:m)' / m * q;
last_sig = find(p_sorted <= thresholds, 1, 'last');

if ~isempty(last_sig)
    pcrit = p_sorted(last_sig);
    sig_vector(valid) = p_valid <= pcrit;
end

sig = reshape(sig_vector, original_size);
end
