function [groups_data, group_names] = generate_rg_inputs_PD_SZ(param,varargin)

%     rawDir = '/ems/elsc-labs/mezer-a/Mezer-Lab/rawData/HUJI/Parkinson_SZ';
%     metadata = fullfile(rawDir,'metadata_do_not_edit.csv');
%     metadata = readtable(metadata);
%     group_names = unique(metadata.ClinicalGroup);
    
    group_names = {'PD','CTL','DYT','GBA','MSA'};
    groups_data = cell(size(group_names));
    for gg = 1:length(group_names)
        group = group_names{gg};
        [subjects,age,sex] = PD_SZ_subjects(group,param);
        Nsubs = length(subjects);
        
        n = Nsubs; % chosen subjects
        subjects = subjects(1:n);
        age = age(1:n);
        sex = sex(1:n);
        
        [map_list, seg_list] = PD_SZ_MapAndSeg(subjects,param,varargin{:});
%         tv_list = PD_SZ_MapAndSeg(subjects, 'MTV');
        
        groups_data{gg}.group_name = group;
        groups_data{gg}.subject_names = subjects;
        groups_data{gg}.age = age;
        groups_data{gg}.sex = sex;
        groups_data{gg}.map_list = map_list;
        groups_data{gg}.seg_list = seg_list;
%         groups_data{gg}.tv_list = tv_list;
    end

%     %% add HUJI Old as Control (mrQ_fixbias)
%     gg = gg+1;
%     [subjects,age,sex] = HUJI_subjects('old',param);    
%     [map_list, seg_list] = huji_MapAndSeg(subjects,param,'noedge',noEdge,'mrQver',2);
%     groups_data{gg}.group_name = 'Old HC';
%     groups_data{gg}.subjects = subjects;
%     groups_data{gg}.age = age;
%     groups_data{gg}.sex = sex;
%     groups_data{gg}.map_list = map_list;
%     groups_data{gg}.seg_list = seg_list;
%     
%     %% add HUJI Young as Control (mrQ_fixbias)
%     gg = gg+1;
%     [subjects,age,sex] = HUJI_subjects('young',param);    
%     [map_list, seg_list] = huji_MapAndSeg(subjects,param,'noedge',noEdge,'mrQver',2);
%     groups_data{gg}.group_name = 'Young HC';
%     groups_data{gg}.subjects = subjects;
%     groups_data{gg}.age = age;
%     groups_data{gg}.sex = sex;
%     groups_data{gg}.map_list = map_list;
%     groups_data{gg}.seg_list = seg_list;
% 
    %% add HUJI Old as Control (New)
    [subjects,age,sex,analysisDir] = HUJI_subjects('old',param);
    
    switch upper(param)
        case 'R1'
            map_sfx = 'mrQ_tests_ed/mrQ_2022_03_09/OutPutFiles_1/BrainMaps/R1_map.nii.gz';
        case {'R2*','R2S','R2STAR'}
            map_sfx = 'multiecho_flash_R2s/R2_mean_2TV.nii.gz';
        case {'TV','MTV'}
            map_sfx = 'mrQ_tests_ed/mrQ_2022_03_09/OutPutFiles_1/BrainMaps/TV_map.nii.gz';
        case upper('T1deg20overT1deg4')
            map_sfx = 'mrQ_tests_ed/mrQ_2022_03_09/semi-quantitative/T1deg20overT1deg4.nii.gz';            
        otherwise
            return
    end
    
    map_list = cellfun(@(x) fullfile(analysisDir,x,map_sfx),subjects,'un',0);
    seg_list = cellfun(@(x) fullfile(analysisDir,x,'FSLFIRST_v2/first_all_none_firstseg_noEdge.nii.gz'), subjects,'un', 0);
    
    
    JoinControlGroups = 1;
    
    if ~JoinControlGroups
        gg = gg+1;
        group_names{gg} = 'Old HC';
        groups_data{gg}.group_name = group_names{gg};
        groups_data{gg}.subject_names = subjects;
        groups_data{gg}.age = age;
        groups_data{gg}.sex = sex;
        groups_data{gg}.map_list = map_list;
        groups_data{gg}.seg_list = seg_list;
    else % join the control groups
        warning('Joining HUJI old subjects with new CTL subjects');
        gg = find(ismember(group_names,'CTL'));
        group_names{gg} = 'HC';
        groups_data{gg}.group_name = group_names{gg};
        groups_data{gg}.subject_names = [groups_data{gg}.subject_names; subjects];
        groups_data{gg}.age = [groups_data{gg}.age; age];
        groups_data{gg}.sex = [groups_data{gg}.sex; sex];
        groups_data{gg}.map_list = [groups_data{gg}.map_list; map_list];
        groups_data{gg}.seg_list = [groups_data{gg}.seg_list; seg_list];
    end
    
    
    
%     %% add HUJI Young as Control (New)
%     gg = gg+1;
%     [subjects,age,sex,analysisDir] = HUJI_subjects('young',param);
%     map_list = cellfun(@(x) fullfile(analysisDir,x,'mrQ_tests_ed/mrQ_2022_03_09/OutPutFiles_1/BrainMaps/R1_map.nii.gz'),subjects,'un',0);
%     if noEdge
%         seg_list = cellfun(@(x) fullfile(analysisDir,x,'FSLFIRST_v2/first_all_none_firstseg_noEdge.nii.gz'), subjects,'un', 0);
%     else
%         seg_list = cellfun(@(x) fullfile(analysisDir,x,'FSLFIRST_v2/first_all_none_firstseg.nii.gz'), subjects,'un', 0);
%     end
%     
%     groups_data{gg}.group_name = 'Young HC new';
%     groups_data{gg}.subject_names = subjects;
%     groups_data{gg}.age = age;
%     groups_data{gg}.sex = sex;
%     groups_data{gg}.map_list = map_list;
%     groups_data{gg}.seg_list = seg_list;
%     
%     
%     group_names = [group_names,'Old HC','Young HC','Old HC new','Young HC new'];
    
end
    