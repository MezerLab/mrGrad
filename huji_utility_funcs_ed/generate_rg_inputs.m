function [groups_data, group_names] = generate_rg_inputs(varargin)

warning('off','all'); % supress warnings

addpath(genpath('/ems/elsc-labs/mezer-a/code/elior/Datasets_utilities'));

[found, dataset, varargin] = argParse(varargin, 'dataset');
if ~found; error('Dataset not specified.'); end
dataset = upper(dataset);

[found, param, varargin] = argParse(varargin, 'param');
if ~found; error('MRI parameter not specified.'); end

switch dataset

    case 'HUJI'
        [groups_data, group_names] = generate_rg_inputs_HUJI(param,varargin{:});
        
    case {'PD_SZ','HUJI_PD'}
        [groups_data, group_names] = generate_rg_inputs_PD_SZ(param,varargin{:});
        
    case {'HCP','HCPY'}
    [groups_data, group_names] = generate_rg_inputs_HCP(param);

    case 'HCP_1MM'
        [groups_data, group_names] = generate_rg_inputs_HCP(param,'1mm');    

    case {'STANFORD','WH'}
        [groups_data, group_names] = generate_rg_inputs_Stanford(param,varargin{:});    

    case 'AHEAD'
        [groups_data, group_names] = generate_rg_inputs_AHEAD(param,varargin{:});

    case {'PPMI','PPMI_TRIO'}
        [groups_data, group_names] = generate_rg_inputs_PPMI(param);

    case {'PPMI_PET','PPMI_ESPREE'} % PPMI Espree (1.5T) with PET data (only PD, no control)
        [groups_data, group_names] = generate_rg_inputs_PPMI_pet(param);
    
    otherwise
        error('Dataset Unknown')
end

% PRINT GROUPS N
fprintf('\nInput data groups:\n');
for jj = 1:length(groups_data)
    fprintf('"%s": N=%d\n',groups_data{jj}.group_name,length(groups_data{jj}.subject_names));
end
fprintf('\n');
warning('on','all');
