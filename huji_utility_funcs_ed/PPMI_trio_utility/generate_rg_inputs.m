function [groups_data, group_names] = generate_rg_inputs(varargin)

addpath(genpath('/ems/elsc-labs/mezer-a/code/elior/Datasets_utilities'));

[found, dataset, varargin] = argParse(varargin, 'dataset');
if ~found; dataset = []; end
dataset = upper(dataset);

[found, param, varargin] = argParse(varargin, 'param');
if ~found; param = []; end

[found, noedge, varargin] = argParse(varargin, 'noEdge');
if ~found; noedge = true; end



[found,~, varargin] = argParse(varargin, 'bsatlas');
if found
    bsatlas_val = 'bsatlas';
else
    bsatlas_val = ' ';
end


switch dataset

    case 'HUJI'
        group_names = {'young','old'};
        groups_data = cell(size(group_names));
        for gg = 1:2
            group = group_names{gg};
            [subjects,age,sex] = HUJI_subjects(group,param,bsatlas_val);
            Nsubs = length(subjects);

            n = Nsubs; % chosen subjects
            subjects = subjects(1:n);
            age = age(1:n);
            sex = sex(1:n);

            [map_list, seg_list] = huji_MapAndSeg(subjects, param,'noedge',noedge,bsatlas_val,varargin{:});
    %         tv_list = huji_MapAndSeg(subjects, 'MTV');

            groups_data{gg}.group_name = group;
            groups_data{gg}.subjects = subjects;
            groups_data{gg}.age = age;
            groups_data{gg}.sex = sex;
            groups_data{gg}.map_list = map_list;
            groups_data{gg}.seg_list = seg_list;
    %         groups_data{gg}.tv_list = tv_list;
        end
        
    case 'PD_SZ'
        [groups_data, group_names] = generate_rg_inputs_PD_SZ(param);
        
    case {'HCP','HCPY'}
    [groups_data, group_names] = generate_rg_inputs_HCP(param);

    case 'HCP_1MM'
        [groups_data, group_names] = generate_rg_inputs_HCP(param,'1mm');    

    case {'STANFORD','WH'}
        [groups_data, group_names] = generate_rg_inputs_Stanford(param,noedge);    

    case 'AHEAD'
        [groups_data, group_names] = generate_rg_inputs_AHEAD(param);

    case {'PPMI','PPMI_TRIO'}
        [groups_data, group_names] = generate_rg_inputs_PPMI(param);

    case {'PPMI_PET','PPMI_ESPREE'} % PPMI Espree (1.5T) with PET data (only PD, no control)
        [groups_data, group_names] = generate_rg_inputs_PPMI_pet(param);
    
    otherwise
        error('Dataset Unknown')
end



