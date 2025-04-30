function asymmetry_scores = mrGrad_asymmetry(data_left,data_right,absolute_flag,norm_flag)
% Compute gradient asymmetry as defined in Drori et al. 2022
% Asymmetry per segment: (L-R) / ((L+r)/2), then averaged across segments
% of interest

% Non-defaults:
% if absolute_flag is set to 1: returns absolute asymmetry value.
% if norm_flag is set to 0: asymmetry is not normalized by mean(L,R)

data_mean = (data_left + data_right)/2;
data_diff = data_left - data_right;

asymmetry_scores = data_diff;

% NORMALIZE ASYMMETRY
if exist('norm_flag','var') && ~norm_flag
    divide_by_mean = 0;
    warning('User choice: asymmetry not normalized.');
%     w = warning('query','last');
%     warning('off',w.identifier);
else
    divide_by_mean = 1;
end
if divide_by_mean
    asymmetry_scores = asymmetry_scores./data_mean;
end

% absolute asymmetry?
if ~exist('absolute_flag','var')
    absolute_flag = 0;
end
if absolute_flag
    asymmetry_scores = abs(asymmetry_scores);
    warning('asymmetry computed as absolute difference');
end
