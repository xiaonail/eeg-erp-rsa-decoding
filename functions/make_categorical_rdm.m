function model_rdm = make_categorical_rdm(categories)
%MAKE_CATEGORICAL_RDM Create a binary model RDM from category labels.
%
%   MODEL_RDM = MAKE_CATEGORICAL_RDM(CATEGORIES) returns 0 for pairs from the
%   same category and 1 for pairs from different categories. CATEGORIES should
%   be ordered in the same way as the neural condition patterns.

categories = categories(:);
model_rdm = double(categories ~= categories');
end
