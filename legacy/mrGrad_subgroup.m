function RG_sub = mrGrad_subgroup(RG, subIdx)

warning(['mrGrad v2.0: The function mrGrad_subgroup() is no longer supported.', newline, ...
       'Please use mrGrad_subset() instead.']);
RG_sub = mrGrad_subset(RG, subIdx);