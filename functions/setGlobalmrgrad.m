function mrgrad_defs = setGlobalmrgrad(varargin)
mrgrad_defs = [];
[found, ROI, varargin] = argParse(varargin, 'ROI');
if ~found
    error('please specify ROI label');
end

[found, roi_names, varargin] = argParse(varargin, 'roi_names');
if ~found
    roi_names = ROI_name(ROI);
end
NROIs = numel(ROI);

[found, segmentingMethod, varargin] = argParse(varargin, 'segmentingMethod');
if ~found; segmentingMethod = 'equidistance'; end
segmentingMethod = lower(segmentingMethod);

if strcmpi(segmentingMethod,'spacing')
    segmentingMethod = 'equidistance';
elseif strcmpi(segmentingMethod,'VoxN')
    segmentingMethod = 'equivolume';
end


[found, stat, varargin] = argParse(varargin, 'stat');
if ~found; stat = 'median'; end


[found, PC, varargin] = argParse(varargin, 'PC');
if ~found; PC = 1:3; end

[found, Nsegs, varargin] = argParse(varargin, 'Nsegs');
if ~found; Nsegs = 7; end
if numel(Nsegs)==1
    Nsegs = repmat(Nsegs,1,length(PC));
end
if numel(Nsegs) ~= numel(PC)
    error('mismatch between lengths of NSEGS and PC');
end

[found, erode_flag, varargin] = argParse(varargin, 'erode');
if ~found; erode_flag = 0; end

[found, invert_flag, varargin] = argParse(varargin, 'invert');
if ~found; invert_flag = false; end

[found, param, varargin] = argParse(varargin, 'param');
if ~found; param = 'unknown_parameter'; end

[found, units, varargin] = argParse(varargin, 'units');
if ~found; units = 'unknown_units'; end

[found_alternative, AlternativeAxes, varargin] = argParse(varargin, 'apply_alternative_axes');


[found_m, max_change, varargin] = argParse(varargin, 'max_change');
if found_m
    if size(max_change,1)==1
        max_change = repmat(max_change,NROIs,1);
    end
    if size(max_change,1)~=NROIs
        error('MAX_CHANGE should have a row for each ROI');
    else
        msg = 'using user''s input priors for PCs directionality sign';
        disp(msg)
    end
end

[found, BL_normalize, varargin] = argParse(varargin, 'BL_normalize');
if ~found; BL_normalize = false; end

[found, isfigs, varargin] = argParse(varargin, 'figures');
if ~found; isfigs = 0; end

[found, outfile, varargin] = argParse(varargin, 'outfile');


mrgrad_defs.ROI = ROI;
mrgrad_defs.roi_names = roi_names;
mrgrad_defs.segmentingMethod = segmentingMethod;
mrgrad_defs.stat = stat;
mrgrad_defs.PC = PC;
mrgrad_defs.Nsegs = Nsegs;
mrgrad_defs.erode_flag = erode_flag;
mrgrad_defs.invert_flag = invert_flag;
mrgrad_defs.param = param;
mrgrad_defs.units = units;
mrgrad_defs.max_change = max_change;
mrgrad_defs.BL_normalize = BL_normalize;
mrgrad_defs.isfigs = isfigs;
mrgrad_defs.outfile = outfile;

if ~isempty(AlternativeAxes)
    if ~isa(AlternativeAxes.seg_list{1},'cell')
        AlternativeAxes.seg_list = {AlternativeAxes.seg_list};
    end
    if numel(AlternativeAxes.ROI)==1
        AlternativeAxes.ROI = repmat(AlternativeAxes.ROI,size(ROI));
    end
    mrgrad_defs.Alternative_seg_list = AlternativeAxes.seg_list;
    mrgrad_defs.Alternative_ROI = AlternativeAxes.ROI;
end
