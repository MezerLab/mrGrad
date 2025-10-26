function mrgrad_defs = setGlobalmrgrad(varargin)
mrgrad_defs = [];
[found, ROI, varargin] = argParse(varargin, 'ROI');
if ~found
    error('please specify ROI label');
end

[found, roi_names, varargin] = argParse(varargin, 'roi_names');
if ~found
    fprintf(' roi_names not provided; mrGrad will use Freesurfer''s look-up table for roi_names.\n')
    roi_names = ROI_name(ROI);
end
NROIs = numel(ROI);

[found, segmentingMethod, varargin] = argParse(varargin, 'segmentingMethod');
if ~found; segmentingMethod = 'equidistance'; end
segmentingMethod = lower(segmentingMethod);


% map old vesions argument naming into current naming
old_args = {'nsegs','outfile', 'segmentingmethod'};
new_args = {'n_segments','output_name', 'segmenting_method'};

idx_ischar = cellfun(@(x) isa(x,'char'),varargin);
string_args = varargin(idx_ischar);
for jj = 1:length(old_args)
    string_args = regexprep(string_args,['(?i)^',old_args{jj},'$'],new_args{jj});
end
varargin(idx_ischar) = string_args;

if strcmpi(segmentingMethod,'spacing')
    segmentingMethod = 'equidistance';
elseif strcmpi(segmentingMethod,'VoxN')
    segmentingMethod = 'equivolume';
end

[found, Parallel, varargin] = argParse(varargin, 'Parallel');
if ~found; Parallel = false; end

[found, stat, varargin] = argParse(varargin, 'stat');
if ~found; stat = 'median'; end


[found, PC, varargin] = argParse(varargin, 'PC');
if ~found; PC = 1:3; end

[found, n_segments, varargin] = argParse(varargin, 'n_segments');
if ~found; n_segments = 7; end
if numel(n_segments)==1
    n_segments = repmat(n_segments,1,length(PC));
end
if numel(n_segments) ~= numel(PC)
    error('mismatch between lengths of n_segments and PC');
end

[found, erode_flag, varargin] = argParse(varargin, 'erode');
if ~found; erode_flag = 0; end

[found, invert_flag, varargin] = argParse(varargin, 'invert');
if ~found; invert_flag = false; end

[found, param, varargin] = argParse(varargin, 'param');


[found, units, varargin] = argParse(varargin, 'units');


[found, ignore_missing, varargin] = argParse(varargin, 'ignore_missing');
if ~found; ignore_missing = false; end

[found_alternative, AlternativeAxes, varargin] = argParse(varargin, 'apply_alternative_axes');


[found_m, max_change, varargin] = argParse(varargin, 'max_change');
if found_m
    if size(max_change,1)==1
        max_change = repmat(max_change,NROIs,1);
    end
    if size(max_change,1)~=NROIs
        error('MAX_CHANGE should have a row for each ROI');
    else
        msg = 'Using user''s input priors for PCs directionality sign';
        disp(msg)
    end
end

[found, BL_normalize, varargin] = argParse(varargin, 'BL_normalize');
if ~found; BL_normalize = false; end

[found, isfigs, varargin] = argParse(varargin, 'figures');
if ~found; isfigs = 0; end

[found, output_mode, varargin] = argParse(varargin, 'output_mode');
if ~found; output_mode = 'default'; end

% Output Name
[found, output_name, varargin] = argParse(varargin, 'output_name');
if ~found
    output_name = sprintf("mrGrad_%s.mat",string(param));
end
output_name = string(output_name);
[a,b,c] = fileparts(output_name);
output_name_path = a;
output_name = b+c;
if ~endsWith(output_name,".mat")
    output_name = output_name + ".mat";
end
if ~startsWith(output_name,"mrGrad")
    output_name = "mrGrad_" + output_name;
end
output_name = strrep(output_name,"_.",".");

% Output directory: (1) output_dir input or (2) path from output_name input
% or (3) current directory
[found, output_dir, varargin] = argParse(varargin, 'output_dir');
if isempty(output_dir)
    output_dir = output_name_path;
end
if isempty(char(output_dir))
    output_dir = fullfile(cd(),"mrGrad_analysis");
    tmp_output_dir = output_dir;
    num = 0;
    while exist(tmp_output_dir,"dir")
        num = num+1;
        tmp_output_dir = strjoin([output_dir,string(num)],"_");
    end
    output_dir = tmp_output_dir;
end

mrgrad_defs.ROI = ROI;
mrgrad_defs.roi_names = roi_names;
mrgrad_defs.segmentingMethod = segmentingMethod;
mrgrad_defs.stat = stat;
mrgrad_defs.PC = PC;
mrgrad_defs.n_segments = n_segments;
mrgrad_defs.erode_flag = erode_flag;
mrgrad_defs.invert_flag = invert_flag;
mrgrad_defs.param = param;
mrgrad_defs.units = units;
mrgrad_defs.max_change = max_change;
mrgrad_defs.BL_normalize = BL_normalize;
mrgrad_defs.isfigs = isfigs;
mrgrad_defs.ignore_missing = ignore_missing;
mrgrad_defs.output_dir = output_dir;
mrgrad_defs.output_name = output_name;
mrgrad_defs.output_mode = output_mode;
mrgrad_defs.parallel = Parallel;

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
