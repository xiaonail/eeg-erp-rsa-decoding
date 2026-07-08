function [rho, p] = compare_rdms(neural_rdm, model_rdm)
%COMPARE_RDMS Spearman correlation of the two RDM upper triangles.

mask = triu(true(size(neural_rdm)), 1);
[rho, p] = corr(neural_rdm(mask), model_rdm(mask), 'type', 'Spearman');
end
