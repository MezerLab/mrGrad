function [map_list, seg_list] = huji_MapAndSeg(subjects, param,varargin)
%--------------------------------------------------------------------------
% INPUT:
%--------------------------------------------------------------------------
%   subjects    cell of HUJI subject names
%
%   param       parameter name (e.g. 'T1')    
%--------------------------------------------------------------------------
% Optional inputs:
%--------------------------------------------------------------------------
%   'freesurfer'    choose freesurfer segmentation instead of default new
%                   FSL FIRST segmentations
%
%   'BSatlas'       choose brainstem atlas segmentation instead of default
%
%   'noEdge'        0 or 1 (default) use eroded segmentations
%--------------------------------------------------------------------------
% OUTPUT:
%--------------------------------------------------------------------------
%   map_list        cell of MRI images paths
%
%   seg_list        cell of segmentation files paths    
%--------------------------------------------------------------------------

% eroded segmentation?
[found, noEdge, varargin] = parseparam(varargin, 'noEdge');
if ~found; noEdge = true; end

% special segmentations?
freesurfer = any(cellfun(@(x) ischar(x) && strcmpi(x, 'freesurfer'), varargin));
BSatlas = any(cellfun(@(x) ischar(x) && strcmpi(x, 'BSatlas'), varargin));
if freesurfer && BSatlas
    error('please choose only one segmentation')
end

% mrQ version?
ver = 1;
mrqVers = {'mrQ', 'mrQ_fixbias'};

mrqDir = mrqVers{ver};
if ~isequal(mrqDir,'mrQ')
    warning('Using mrQ version: %s',mrqDir);
end

if ismember(upper(param), {'R2S','R2STAR'})
    param = 'R2*';
elseif ismember(upper(param), {'MTVC'})
    param = 'MTV_c';
end

% P = readtable('/ems/elsc-labs/mezer-a/elior.drori/Code/Striatum/Datasets/huji/HUJI_maps_table.csv');

% MRI parameter look-up table
HighRes= [181 217 181];
LowRes = [121 145 121];
P = table({},{},[],{},'VariableNames',{'Name','Units','Resolution','Path'});
P(1,:) = {'T1','s',HighRes,fullfile(mrqDir,'OutPutFiles_1/BrainMaps/T1_map_Wlin.nii.gz')};
P(2,:) = {'R1','s^{-1}',HighRes,fullfile(mrqDir,'OutPutFiles_1/BrainMaps/R1_map_Wlin.nii.gz')};
P(3,:) = {'MTV','fraction',HighRes,fullfile(mrqDir,'OutPutFiles_1/BrainMaps/TV_map.nii.gz')};
P(4,:) = {'MTV_c','fraction',HighRes,fullfile(mrqDir,'OutPutFiles_1/BrainMaps/TV_correctedForT2s.nii.gz')};
P(5,:) = {'B1','s',HighRes,fullfile(mrqDir,'OutPutFiles_1/BiasMap/B1_Map.nii.gz')};
P(6,:) = {'R2*','ms^{-1}',HighRes,'multiecho_flash_R2s/R2_mean_2TV.nii.gz'};
P(7,:) = {'MTsat','p.u.',HighRes,'MT/MT_sat.nii.gz'};
P(8,:) = {'MD','mm^2/s',LowRes,'Dif_fsl_preprocessed/eddy/aligned2T1/dtiInit/dti94trilin/bin/MD.nii.gz'};
P(9,:) = {'FA','unitless',LowRes,'Dif_fsl_preprocessed/eddy/aligned2T1/dtiInit/dti94trilin/bin/FA.nii.gz'};
% P(10,:) = {'T1_2DTI','s',LowRes,'Dif_fsl_preprocessed/eddy/aligned2T1/dtiInit/dti94trilin/bin/T1_map_Wlin_2DTI_resamp.nii.gz'};
% P(11,:) = {'R1slope','slope',HighRes,'MDMmap/R1slope.nii.gz'};
% P(12,:) = {'MTsatslope','slope',HighRes,'MDMmap/MTsatslope.nii.gz'};
% P(13,:) = {'R1_rsqr','slope',HighRes,'MDMmap/R1slope_Rsqr.nii.gz'};
% P(14,:) = {'MTsat_rsqr','slope',HighRes,'MDMmap/MTsatslope_Rsqr.nii.gz'};

% MRI IMAGE FILES
idx = strcmpi(P.Name,param);
im_path = P.Path{idx};
im_res = P.Resolution(idx,:);
analysisDir = '/ems/elsc-labs/mezer-a/Mezer-Lab/analysis/HUJI/Calibration/Human/';
map_list = cellfun(@(x) fullfile(analysisDir,x,im_path), subjects,'UniformOutput', false);

% SEGMENTATION FILES
seg_list_fs = cellfun(@(x) fullfile(analysisDir,x,'freesurfer/seg.nii.gz'), subjects,'un', 0);
seg_list_fsl = cellfun(@(x) fullfile(analysisDir,x,'FSLFIRST_v2/first_all_none_firstseg.nii.gz'), subjects,'un', 0);
seg_list_bs = cellfun(@(x) fullfile(analysisDir,x,'BSAtlasSpace/MASSP/massp_massp-label_2sub.nii.gz'), subjects,'un', 0);

if freesurfer
    seg_list = seg_list_fs;
elseif BSatlas
    seg_list = seg_list_bs;
else
    seg_list = seg_list_fsl;
end

if isequal(im_res,[121 145 121])
    S = isequal(seg_list,seg_list_fs);
    seg_path = 'Dif_fsl_preprocessed/eddy/aligned2T1/dtiInit/dti94trilin/bin/';
    seg_list_LowRes = cellfun(@(x) fullfile(analysisDir,x,seg_path,'freesurfer/seg_2dti.nii.gz'), subjects,'un', 0);
    seg_list = seg_list_LowRes;
    if ~S
        warning('Ignoring user''s preference: using freesurfer segmentation resampled to DTI');
    end
end

if noEdge && ~BSatlas
    seg_list = cellfun(@(x) strrep(x,'.nii','_noEdge.nii'),seg_list,'un',0);
%     warning('using no-edge segmentation');
end




