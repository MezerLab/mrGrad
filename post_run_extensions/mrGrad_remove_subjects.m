function rg_new = mrGrad_remove_subjects(rg,idx_rm)
% This function is just using mrGrad_subgroup
rg_new = mrGrad_subgroup(rg,~idx_rm);

