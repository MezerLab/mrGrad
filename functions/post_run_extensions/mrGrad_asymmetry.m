function asymmetry_scores = mrGrad_asymmetry(data_left,data_right,absolute_flag,norm_flag)
% mrGrad_asymmetry computes gradient asymmetry based on Drori et al. (2022).
%
% Asymmetry is computed for each segment as: 
%   (L - R) / ((L + R) / 2), where L = data_left and R = data_right.
%
% INPUTS:
%   - data_left: Vector or matrix of left hemisphere data.
%   - data_right: Vector or matrix of right hemisphere data.
%   - absolute_flag (optional): boolean, returns absolute asymmetry values. Default is false.
%   - norm_flag (optional): boolean, normalizes asymmetry by the mean of (L, R). Default is true.
%
% OUTPUT:
%   - asymmetry_scores: Vector of computed asymmetry scores for each segment.

absolute_flag = exist('absolute_flag','var') && ~isempty(absolute_flag) && absolute_flag;
norm_flag = ~exist('norm_flag','var') || isempty(norm_flag) || exist('norm_flag','var') && norm_flag;

% Compute the mean and difference between left and right data
data_mean = (data_left + data_right)/2;
data_diff = data_left - data_right;

% Initialize asymmetry scores with the difference between left and right
asymmetry_scores = data_diff;

% Normalize the asymmetry if norm_flag is set
if norm_flag
    asymmetry_scores = asymmetry_scores./data_mean;
end

% Compute absolute asymmetry if absolute_flag is set
if absolute_flag
    asymmetry_scores = abs(asymmetry_scores);
end
