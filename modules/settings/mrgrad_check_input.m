function [Data, mrgrad_defs] = mrgrad_check_input(Data, mrgrad_defs)
% Validates mrGrad v2.0 main input.

if ~isa(Data,'struct')
    error('DATA input should be a struct. Please refer to mrGrad() for help.');
end

% make sure obligatory input lists exist
obfields = {'map_list','seg_list'};
if any(~isfield(Data,obfields))
    error('DATA fields ''map_list'',''seg_list'' are obligatory');
end

% make sure input lists are cell arrays
Data.seg_list = cellstr(Data.seg_list);
Data.map_list = cellstr(Data.map_list);

N = length(Data.seg_list);
if ~N
    error('Please provide at least one segmentation/mask path');
end
    
% Verify matching lengths of input lists
if size(Data.map_list,1)~=N
    error('height of map_list should match height of seg_list (number of subject records)');
end

% verify matching number of map_lists and number of parameter names
n_map_lists = size(Data.map_list,2);
if isempty(mrgrad_defs.parameter_names)
    mrgrad_defs.parameter_names = cellstr("unknown_parameter" + (1:n_map_lists));
elseif ~isequal(length(mrgrad_defs.parameter_names), n_map_lists)
    error('length of PARAMETER_NAMES should match the width of DATA.map_list');
end
if isempty(mrgrad_defs.units)
    mrgrad_defs.units = repmat({'unknown_units'},1,n_map_lists);
elseif ~isequal(length(mrgrad_defs.units), n_map_lists)
    error('number of UNITS strings should match the width of DATA.map_list');
end

% SUBJECT IDS
if ~isfield(Data,'subject_ids')
    maxDigits = max([floor(log10(N))+1,3]);
    fmt = ['sub-%0' num2str(maxDigits) 'd'];
    Data.subject_ids = arrayfun(@(num) sprintf(fmt,num),(1:N)','un',0);
end

warn_prefx = "\n ";

% Verify existence of input files (segmentations)
idx = cellfun(@(x) ~exist(x,'file'),Data.seg_list);
if all(idx)
    fprintf(warn_prefx); warn_prefx = " ";
    error('None of the input segmentation files exist. Please check your file paths.');
elseif any(idx)
    if mrgrad_defs.allow_missing
        fprintf(warn_prefx); warn_prefx = " ";
        warning('Some input segmentation files are missing. Proceeding with available data (allow_missing = true).');
    else
        fprintf(warn_prefx); warn_prefx = " ";
        error('One or more input segmentation files are missing. To continue despite missing files, set allow_missing = true.');
    end
end

% Verify existence of input files (maps)
idx = cellfun(@(x) ~exist(x,'file'),Data.map_list);
if all(idx)
    fprintf(warn_prefx); warn_prefx = " ";
    error('None of the input image files exist. Please check your file paths.');
elseif any(idx)
    if mrgrad_defs.allow_missing
        fprintf(warn_prefx); warn_prefx = " ";
        warning('Some input image files are missing. Proceeding with available data (allow_missing = true).');
    else
        fprintf(warn_prefx); warn_prefx = " ";
        error('One or more input image files are missing. To continue despite missing files, set allow_missing = true.');
    end
end
