function [subjects,age,sex,analysisDir] = HUJI_subjects(varargin)
%--------------------------------------------------------------------------
% OPTIONAL INPUT:
%--------------------------------------------------------------------------
%     'young' / 'old':      specify age group
%
%     'M' / 'F':            specify sex group
%
%     'R2*':                only subjects with R2* map
%
%     'BSatlas':            only subjects with brainstem atlas segmentation
%
%     'all':                do not remove corrupt subjects
%--------------------------------------------------------------------------
% OUTPUT:
%--------------------------------------------------------------------------
%     The function returns a list of HUJI subjects who have MDM maps, and
%     R2* maps (optional), that belong to a specified age group (optional)
%     and sex group (optional). 
%
%     The function uses Shir's function gen_Subjects.m
%     Elior
%--------------------------------------------------------------------------
varargin = varargin(cellfun(@(x) isa(x,'str')||  isa(x,'char'), varargin)); % remove non-string arguments
analysisDir='/ems/elsc-labs/mezer-a/Mezer-Lab/analysis/HUJI/Calibration/Human/'; % qMRI data path
analysis_d='/ems/elsc-labs/mezer-a/Mezer-Lab/projects/analysis/MTV_Vs_qMRI'; % path to save MDM analysis
T = gen_Subjects(analysisDir,analysis_d);

if any(strcmpi(varargin,'ALL'))
    all_flag = true;
else
    all_flag = false;
end

% MDM index
mdm_idx = cellfun(@(s) logical(exist(fullfile(analysisDir,s,'MDMmap'),'dir')),{T.sub_path}');

% R2* index
if any(ismember(upper(varargin),{'R2*','R2S','R2STAR'}))
    r2s_flag = true;
    r2s_idx = cellfun(@(s) logical(exist(fullfile(analysisDir,s,'multiecho_flash_R2s/R2_mean_2TV.nii.gz'),'file')),{T.sub_path}');
else
    r2s_flag = false;
    r2s_idx = true(length(T),1);
end

% BSatlas index
if any(ismember(upper(varargin),{'BSATLAS'}))
    bs_flag = true;
    bs_idx = logical(cellfun(@(x) exist(fullfile(analysisDir,x,'BSAtlasSpace/MASSP/massp_massp-label_2sub.nii.gz'),'file'), {T.sub_path}'));
else
    bs_flag = false;
    bs_idx = true(length(T),1);
end

% age index
if any(startsWith(varargin,'YOU','IgnoreCase',true))
    age_idx = [T.age]' <= 31;
elseif any(startsWith(varargin,'OLD','IgnoreCase',true))
    age_idx = [T.age]' >= 57;
else
    age_idx = true(length(T),1);
end

% sex index
if any(strcmpi(varargin,'M'))
    sex_idx = strcmp({T.sex}','M');
elseif any(strcmpi(varargin,'F'))
    sex_idx = strcmp({T.sex}','F');
else
    sex_idx = true(length(T),1);
end

%% intersection
% idx = mdm_idx & r2s_idx & age_idx & sex_idx; warning('using N=20 young, N=17 old');
idx = r2s_idx & age_idx & sex_idx & bs_idx;
T(~idx) = [];

%% remove specific subjects
% remove old female subject H046_NB due to right-caudate deterioration
if ~all_flag
    
    % outlier to remove always
    art_idx = strcmp('H046_NB/2018_04_10',{T.sub_path});
    if nnz(art_idx)>0
        T(art_idx)=[];
        warning('Subject H046_NB removed due to right-caudate deterioration');
    end

%     % outlier to remove for R2*
%     if r2s_flag
%         art_idx = strcmp('H039_JH/2018_03_21',{T.sub_path}); % R2* outlier
%         if nnz(art_idx)>0
%             T(art_idx)=[];
%             warning('Subject H039_JH removed due to R2* outlier map');
%         end
%     end
    warning('did not remove R2* outlier H039_JH')
end
%% Output
subjects = {T.sub_path}';
age = [T.age]';
sex = {T.sex}';







