function mask = ROImask(seg,ROIs,use_erosion)
%--------------------------------------------------------------------------
% INPUTS:
%
%   SEG             path to segmentation file
%
%   ROIS            ROI labels (freesurfer)     e.g. 12
%   
%   USE_EROSION     logical specifying whether to erode outer shel of ROI
%                   (optional. default: 0) 
%
% OUTPUT:
%       
%   MASK    binary mask of ROI
%
% (C) Elior Drori, Mezer lab, the Hebrew University of Jerusalem, Israel, Copyright 2021
%--------------------------------------------------------------------------

if ~exist('readFileNifti.m','file')
    error('please download ''vistasoft'' package and add to MATLAB path (https://github.com/vistalab/vistasoft)');
end

if ~exist('use_erosion','var')
    use_erosion = false;
end

if ischar(seg)
    seg = readFileNifti(seg);
end
if isstruct(seg)
    seg = seg.data;
end
seg = double(seg);
mask = zeros(size(seg));

for r = 1:length(ROIs)
    mask(seg == ROIs(r)) = 1;
end

% MASK EROSION TO AVOID PVE
if use_erosion
    Nerosions = 1; % control the number of erosions (1 is enough);

    a = [0,0,0;...
         0,1,0;...
         0,0,0];

    b = [0,1,0;...
         1,1,1;...
         0,1,0];

    SE = cat(3,a,b,a);
    for j = 1:Nerosions
        mask = imerode(mask,SE);
    end
    if Nerosions>1
        warning('used multiple (%d) erosions on ROI',Nerosions)
    end
end
mask = logical(mask);
end
