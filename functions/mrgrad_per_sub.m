function vout = mrgrad_per_sub(qmap,mask,varargin)
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
%   'Nsegs':         followed by a number of wanted segments
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

% % specify PC
% [found, pc, varargin] = argParse(varargin, 'PC');
% if ~found; pc = 1; end
[found, sampling_method, varargin] = argParse(varargin, 'sampling_method');
if ~found; sampling_method = 'equidistance'; end
[found, Nsegs, varargin] = argParse(varargin, 'Nsegs');
if ~found; Nsegs = 7; end
[found, stat, varargin] = argParse(varargin, 'stat');
if ~found; stat = 'median'; end
[found, BL_normalize, varargin] = argParse(varargin, 'BL_normalize');
if ~found; BL_normalize = false; end
[found, isfigs, varargin] = argParse(varargin, 'figures');
if ~found; isfigs = 0; end

%% calculate average qMRI (e.g. T1) in each segment
%==========================================================================
% MAIN FUNCTION EXECUTION
%--------------------------------------------------------------------------
if strcmpi(sampling_method,'equidistance')
    axes_data = RG_axes(mask,varforward{:});
elseif strcmpi(sampling_method,'equivolume')
    axes_data = RG_axes_equalVol(mask,varforward{:});
end

%==========================================================================
linearInd = axes_data.segment_inds_linear;

if strcmpi(stat,'mean')
    rg = cellfun(@(x) mean(qmap(x)), linearInd);
    sub_std = cellfun(@(x) std(qmap(x)), linearInd);
    error_type = 'std';
elseif strcmpi(stat,'median')    
    rg = cellfun(@(x) median(qmap(x)), linearInd);
    sub_std = cellfun(@(x) mad(qmap(x)), linearInd);
    error_type = 'mad';
    
    % subtract median value of structure
    if BL_normalize        
        % save mean value for normalization option
        m_qmap = qmap(mask>0);
        m_qmap = median(m_qmap);
        rg = rg-m_qmap;
        warning('individual baselines are subtracted from gradients');
    end
end
%% FIG gradient
if isfigs
figure;
plot(1:Nsegs, rg,'--o');
h1=title('MRI gradient alog an axis of an ROI, single subject');
xlabel(['Segments along PC',num2str(1)],'FontSize',22);
ylabel('mean qMRI','FontSize',22);
h1.FontSize=22;
end
%% OUTPUT
vout = axes_data;
vout.function = rg;
vout.error = sub_std;
vout.err_type = error_type;
end