function Data = mrgrad_check_input(Data)

if isa(Data,'struct')
    Data = {Data};
end
Ngroups = numel(Data);

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
    
    % check for non-existing inputs
    idx = cellfun(@(x) ~exist(x,'file'),Data{gg}.map_list);
    if any(idx)
        error('one or all input image files not exist.')
    end
    idx = cellfun(@(x) ~exist(x,'file'),Data{gg}.seg_list);
    if any(idx)
        error('one or all input segmentation files not exist.')
    end
end
