function Data = mrgrad_check_input(Data, mrgrad_defs)

if isa(Data,'struct')
    Data = {Data};
end
Ngroups = numel(Data);

warn_prefx = "\n ";

for gg = 1:Ngroups
    
    % make sure obligatory input lists exist
    obfields = {'map_list','seg_list'};
    if any(~isfield(Data{gg},obfields))
        error('DATA fields ''map_list'',''seg_list'' are obligatory');
    end
    if isempty(Data{gg}.map_list)
        error('some map lists are empty');
    end

    % make sure input lists are cell arrays
    if isa(Data{gg}.map_list,"string")
        Data{gg}.map_list = cellstr(Data{gg}.map_list);
    elseif isa(Data{gg}.map_list,"char")
        Data{gg}.map_list = {Data{gg}.map_list};
    end
    if isa(Data{gg}.seg_list,"string")
        Data{gg}.seg_list = cellstr(Data{gg}.seg_list);
    elseif isa(Data{gg}.seg_list,"char")
        Data{gg}.seg_list = {Data{gg}.seg_list};
    end

    % make sure lengths of input lists match
    if ~isequal(length(Data{gg}.map_list),length(Data{gg}.seg_list))
        error('number of maps and segmentations don''t match');
    end
    
    % GROUP NAME
    if ~isfield(Data{gg},'group_name')
        Data{gg}.group_name = sprintf('SubjectGroup_%d',gg);
    end

    % SUBJECT NAMES
    Nsubs = length(Data{gg}.map_list);
    if ~isfield(Data{gg},'subject_names')
        maxDigits = max([floor(log10(Nsubs))+1,3]);
        fmt = ['sub-%0' num2str(maxDigits) 'd'];
        Data{gg}.subject_names = arrayfun(@(num) sprintf(fmt,num),(1:Nsubs)','un',0);
    end

    % check for non-existing inputs
    idx = cellfun(@(x) ~exist(x,'file'),Data{gg}.map_list);
    if all(idx)
        fprintf(warn_prefx); warn_prefx = " ";
        error('None of the input image files exist for Group: %s. Please check your file paths.',Data{gg}.group_name);
    elseif any(idx) && ~mrgrad_defs.ignore_missing
        fprintf(warn_prefx); warn_prefx = " ";
        error('One or more input image files are missing for Group: %s. To continue despite missing files, set ignore_missing = true.',Data{gg}.group_name);
    elseif mrgrad_defs.ignore_missing
        fprintf(warn_prefx); warn_prefx = " ";
        warning('Some input image files are missing for Group: %s. Proceeding with available data (ignore_missing = true).',Data{gg}.group_name);
    end

    idx = cellfun(@(x) ~exist(x,'file'),Data{gg}.seg_list);
    if all(idx)
        fprintf(warn_prefx); warn_prefx = " ";
        error('None of the input segmentation files exist for Group: %s. Please check your file paths.',Data{gg}.group_name);
    elseif any(idx) && ~mrgrad_defs.ignore_missing
        fprintf(warn_prefx); warn_prefx = " ";
        error('One or more input segmentation files are missing for Group: %s. To continue despite missing files, set ignore_missing = true.',Data{gg}.group_name);
    elseif mrgrad_defs.ignore_missing
        fprintf(warn_prefx); warn_prefx = " ";
        warning('Some input segmentation files are missing for Group: %s. Proceeding with available data (ignore_missing = true).',Data{gg}.group_name);
    end
    
end
