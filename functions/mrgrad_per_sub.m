function vout = mrgrad_per_sub(images,mask,varargin)
%--------------------------------------------------------------------------
% INPUTS:
%
%   qmap:   masked 3D qmri map
%
%   mask:   ROI mask
%
%--------------------------------------------------------------------------
% OPTIONAL INPUTS:
%
%   'PC':            followed by a number of principal component (e.g. 1)
%                    default: 1
%
%   'n_segments':    followed by the number of desired segments
%                    default: 7
%
%   'figures':       followed by 1 or 0, to either produce figures or not.
%                    default: 1
%
%   'stat':        followed by 'mean' or 'median'
%                  Specifies the statistic to obtain from each segment.
%
%   'segmentingMethod': followed by 'equidistance' (default) or
%               'equivolume', to specify whether to use equally-spaced
%               segments or segemnts of equal voxel count.
%
%   'BL_normalize': followed by 1 or 0 (default). Substract baseline from
%                   gradients. basline = median value in the structure.
%                   This intends to eliminate absolute differences between
%                   subjects, if needed.
%
% (C) Mezer lab, the Hebrew University of Jerusalem, Israel, Copyright 2021
%--------------------------------------------------------------------------

varforward = varargin;

[found, sampling_method, varargin] = argParse(varargin, 'sampling_method');
if ~found; sampling_method = 'equidistance'; end
[found, stat, varargin] = argParse(varargin, 'stat');
if ~found; stat = 'median'; end

%==========================================================================
% (1) COMPUTE MRGRAD ROI AXES
%==========================================================================
% MAIN FUNCTION CALL mrgrad_axes.m
%--------------------------------------------------------------------------
if strcmpi(sampling_method,'equidistance')
    axes_data = mrgrad_axes_equidistance(mask,varforward{:});
elseif strcmpi(sampling_method,'equivolume')
    axes_data = mrgrad_axes_equivolume(mask,varforward{:});
end
%==========================================================================
% (2) COMPUTE IMAGE INTENSITY PROFILE (GRADIENT)
%==========================================================================
linearInd = axes_data.segment_inds_linear;

switch lower(stat)
    case 'mean'
        stat_func = @mean;
        error_func = @std;
        error_name = 'std';
    case 'median'
        stat_func = @median;
        error_func = @mad;
        error_name = 'mad';
end
rg = single(nan(axes_data.N_segments,numel(images)));
sub_std = rg;
global_stat = single(nan(numel(images),2));
for jj = 1:numel(images)
    im = images{jj};
    rg(:,jj) = cellfun(@(x) stat_func(im(x)), linearInd);
    sub_std(:,jj) = cellfun(@(x) error_func(im(x)), linearInd);

    global_stat(jj,:) = [stat_func(axes_data.all_inds_linear), error_func(axes_data.all_inds_linear)];
end

% OUTPUT
vout = axes_data;
vout.function = rg;
vout.error = sub_std;
vout.err_type = error_name;
vout.whole_roi_stats = global_stat;
end