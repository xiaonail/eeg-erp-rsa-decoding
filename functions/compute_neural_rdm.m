function rdm = compute_neural_rdm(patterns)
%COMPUTE_NEURAL_RDM RDM (1 - Pearson r) from a condition-by-feature matrix.

rdm = 1 - corr(patterns', 'type', 'Pearson');
end
