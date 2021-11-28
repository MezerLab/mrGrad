function [groups_data, group_names] = generate_rg_inputs(varargin)

[found, dataset, varargin] = parseparam(varargin, 'dataset');
if ~found; dataset = []; end
dataset = upper(dataset);

[found, param, varargin] = parseparam(varargin, 'param');
if ~found; param = []; end

[found, noedge, varargin] = parseparam(varargin, 'noEdge');
if ~found; noedge = true; end


% HUJI
if strcmpi(dataset,'HUJI')
    group_names = {'young','old'};
    groups_data = cell(size(group_names));
    for gg = 1:2
        group = group_names{gg};
        [subjects,age,sex] = HUJI_subjects(group,param,'bsatlas');
%         [subjects,age,sex] = HUJI_subjects(group,'all'); % With H046
%         [subjects,age,sex] = HUJI_subjects(group,'R2s',param); warning('only subjects with R2*'); % only subjects  with R2*
        Nsubs = length(subjects);
        
        n = Nsubs; % chosen subjects
        subjects = subjects(1:n);
        age = age(1:n);
        sex = sex(1:n);
        
        [map_list, seg_list] = huji_MapAndSeg(subjects, param,'noedge',noedge,'bsatlas');
        tv_list = huji_MapAndSeg(subjects, 'MTV');

        groups_data{gg}.group_name = group;
        groups_data{gg}.subjects = subjects;
        groups_data{gg}.age = age;
        groups_data{gg}.sex = sex;
        groups_data{gg}.map_list = map_list;
        groups_data{gg}.seg_list = seg_list;
        groups_data{gg}.tv_list = tv_list;    
    end

% HUJI PD_SZ
elseif strcmpi(dataset,'PD_SZ')
    [groups_data, group_names] = generate_rg_inputs_PD_SZ(param);

% % Hadassah PD T2w/T1w 2mm
% elseif strcmpi(dataset,'PD')
%     [groups_data, group_names] = generate_rg_inputs_PD(param);

% HCP Young Adults
elseif ismember(dataset,{'HCP','HCPY'})
    [groups_data, group_names] = generate_rg_inputs_HCP(param);
    
elseif strcmpi(dataset,'HCP_1MM')
    [groups_data, group_names] = generate_rg_inputs_HCP(param,'1mm');    
    
% Stanford
elseif ismember(dataset,{'STANFORD','WH'})
    [groups_data, group_names] = generate_rg_inputs_Stanford(param,noedge);    
    
% AHEAD
elseif strcmpi(dataset,'AHEAD')
    [groups_data, group_names] = generate_rg_inputs_AHEAD(param);
    
% Parkinson PPMI
elseif ismember(dataset,{'PPMI','PPMI_TRIO'})
    [groups_data, group_names] = generate_rg_inputs_PPMI(param);
    
% Parkinson PPMI Espree (1.5T) with PET data (only PD, no control)
elseif ismember(dataset,{'PPMI_PET','PPMI_ESPREE'})
    [groups_data, group_names] = generate_rg_inputs_PPMI_pet(param);
end



