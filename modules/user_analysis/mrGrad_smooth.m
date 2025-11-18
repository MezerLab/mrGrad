function RG_smooth = mrGrad_smooth(RG, smoothing_window)
% mrGrad_smooth
%
% DESCRIPTION:
%   Smooths the mrGrad results in each ROI of an mrGrad structure using
%   Gaussian smoothing along the first dimension (axis position).
%
% SYNTAX:
%   RG_smooth = mrGrad_smooth(RG, smoothing_window)
%
% INPUTS:
%   RG : struct
%       mrGrad-style struct with fields:
%           RG.(ROI).Results.(subfield).y   - table of numeric values
%
%   smoothing_window : scalar
%       Window length for smoothing (passed to smoothdata).
%
% OUTPUT:
%   RG_smooth : struct
%       Same as RG, but all `.Results.(subfield).y` tables have been smoothed.
%
% NOTES:
%   - Applies smoothdata(...,2,"gaussian",smoothing_window) row-wise.
% -------------------------------------------------------------------------
% Input Validation
% -------------------------------------------------------------------------
if ~isstruct(RG)
    error('Input RG must be a struct.');
end
if nargin < 2 || isempty(smoothing_window)
    error('You must provide a smoothing_window value.');
end

RG_smooth = RG;

% -------------------------------------------------------------------------
% Iterate through each ROI
% -------------------------------------------------------------------------
roi_names = fieldnames(RG);
for i = 1:numel(roi_names)
    roi = roi_names{i};

    subfields = fieldnames(RG.(roi).Results);
    for j = 1:numel(subfields)
        subf = subfields{j};
        y = RG.(roi).Results.(subf).y;

        % Smooth across columns (i.e., along dimension 2)
        y_smooth = y;
        y_smooth{:,:} = smoothdata(y{:,:}, 2, "gaussian", smoothing_window);

        RG_smooth.(roi).Results.(subf).y = y_smooth;

        % Re-compute mean, std, sem
        RG_smooth.(roi).Results.(subf).y_mean = mean(RG_smooth.(roi).Results.(subf).y{:,:},1,"omitnan");
        RG_smooth.(roi).Results.(subf).y_std = std(RG_smooth.(roi).Results.(subf).y{:,:},[],1,"omitnan");
        n_actual = sum(~isnan(RG_smooth.(roi).Results.(subf).y{:,:}),1);
        RG_smooth.(roi).Results.(subf).y_sem = RG_smooth.(roi).Results.(subf).y_std ./ sqrt(n_actual);


    end
end
end
