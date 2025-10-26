function mask = ROImask(segmentation,roi_labels,use_erosion)
%--------------------------------------------------------------------------
% INPUTS:
%
%   segmentation             path to segmentation file
%
%   roi_labels      segmentation's labels to mask  e.g. [12,51]
%   
%   use_erosion     logical specifying whether to erode outer shel of ROI
%                   (optional. default: 0) 
%
% OUTPUT:
%       
%   MASK    binary mask of ROI
%
% (C) Elior Drori, Mezer lab, the Hebrew University of Jerusalem, Israel, Copyright 2021
%--------------------------------------------------------------------------

if any(arrayfun(@(type) isa(segmentation,type),["char", "string", "cell"]))
    segmentation = niftiread(string(segmentation));
end

% generate mask from segmentation and labels
mask = ismember(single(segmentation),single(roi_labels));

% erode mask to aviod pve
use_erosion = exist('use_erosion','var') && use_erosion;
if use_erosion
    mask = imerode(mask,strel("sphere",1));
end