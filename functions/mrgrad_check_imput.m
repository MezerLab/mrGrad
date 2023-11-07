function Data = mrgrad_check_imput(Data)

if isa(Data,'struct')
    Data = {Data};
end
Ngroups = numel(Data);

% make sure obligatory input and files exist
for gg = 1:Ngroups
    obfields = {'map_list','seg_list'};
    if any(~isfield(Data{gg},obfields))
        error('DATA fields ''map_list'',''seg_list'' are obligatory');
    end
    if isempty(Data{gg}.map_list)
        error('some map lists are empty');
    end
    if ~isequal(length(Data{gg}.map_list),length(Data{gg}.seg_list))
        error('number of maps and segmentations don''t match');
    end
    
    idx = cellfun(@(x) ~exist(x,'file'),Data{gg}.map_list);
    if any(idx)
        error('one or all input image files not exist.')
    end
    idx = cellfun(@(x) ~exist(x,'file'),Data{gg}.seg_list);
    if any(idx)
        error('one or all input segmentation files not exist.')
    end
end
