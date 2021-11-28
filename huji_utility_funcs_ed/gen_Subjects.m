function [subjects,run_on]=gen_Subjects(analysis,analysis_d)

% This function will generate a list of all huji subjects used for the MDM paper 
% (the rest of the subject might not have the appropriate data). 
% To use this function:
% analysis='/ems/elsc-labs/mezer-a/Mezer-Lab/analysis/HUJI/Calibration/Human/'; % qMRI data path
% analysis_d='/ems/elsc-labs/mezer-a/Mezer-Lab/projects/analysis/MTV_Vs_qMRI'; % path to save 
% [subjects,~]=gen_Subjects(analysis,analysis_d); % generate list of subjects

%% subjects list
analysisDir = analysis;
subjects_dir = dir(analysisDir);
subjects_name = {subjects_dir.name};
subjects_name = subjects_name(~cellfun(@isempty,regexp(subjects_name, 'H\d\d\d_[A-Z][A-Z]')));

% Take only subject 009 and on
subjects_path = subjects_name([1,find(strcmp(subjects_name, 'H009_IR')):52]);
subjects_path(strcmp(subjects_path,'H010_AG')) = []; % Ignore this subjects, because his mrQ is not in ACPC
subjects_path(strcmp(subjects_path,'H011_GP')) = []; % Ignore this subjects, because he doesn't have mrQ data
subjects_path(strcmp(subjects_path,'H014_ZW')) = []; % Ignore this subjects, because he doesn't have mrQ data
subjects_path(strcmp(subjects_path,'H029_ON')) = []; % Ignore this subjects, because he doesn't have data
age=[26 30	27	27	27	27	26	57	69	31	26	27	61	63	75	57	27	23	25	68	70	71	65	63	73	24	31	73	25	65	27	29	26	31	77	27	26	24 67 75 69 ];
sex={'F' 'M'	'M'	'F'	'M'	'F'	'F'	'M'	'M'	'M'	'M'	'M'	'F'	'M'	'M'	'M'	'M'	'F'	'F'	'M'	'M'	'F'	'M'	'M'	'M'	'M'	'M'	'M'	'F'	'F'	'F'	'F'	'F'	'M'	'F'	'M'	'M'	'F'	'M'	'M'	'F'};

for ii=1:length(subjects_path)
    subjects(ii).Subject_D=fullfile(analysis_d,subjects_path{ii});
    ind=find(strcmp({subjects_dir.name},subjects_path{ii}));
    path=[analysisDir '/' subjects_dir(ind).name];
    date=dir(path);
    date={date.name};
    if any(strcmp(date,'readme'))
        fid = fopen([path,'/readme']);
        A= fscanf(fid, '%s');
        subjects(ii).sub_path=[subjects_path{ii},'/',A];
    else
        date = date(~cellfun(@isempty,regexp(date, '\d')));
        subjects(ii).sub_path=[subjects_path{ii},'/',date{1}];
    end
    subjects(ii).age=age(ii);
    subjects(ii).sex=sex{ii};
    
    sus_path=fullfile('/ems/elsc-labs/mezer-a/Mezer-Lab/rawData/HUJI/Calibration/Human', subjects(ii).sub_path,'/nifti_dicm2nii/susceptibility');
    if exist(sus_path)
        subjects(ii).T2s=1;
    else
        subjects(ii).T2s=0;
    end
    ME_path=fullfile('/ems/elsc-labs/mezer-a/Mezer-Lab/analysis/HUJI/Calibration/Human/',subjects(ii).sub_path,'multiecho_flash_R2s');
    if exist(ME_path)
        subjects(ii).T2f=1;
    else
        subjects(ii).T2f=0;
    end
        T2_path=fullfile('/ems/elsc-labs/mezer-a/Mezer-Lab/analysis/HUJI/Calibration/Human/',subjects(ii).sub_path,'/T2/T2map.nii.gz');
    if exist(T2_path)
        subjects(ii).T2=1;
    else
        subjects(ii).T2=0;
    end
    MD_path=fullfile('/ems/elsc-labs/mezer-a/Mezer-Lab/analysis/HUJI/Calibration/Human/',subjects(ii).sub_path,'DTIoutput_1.5mm');
    if exist(MD_path)
        subjects(ii).MD=1;
    else
        subjects(ii).MD=0;
    end
        MD_path=fullfile('/ems/elsc-labs/mezer-a/Mezer-Lab/analysis/HUJI/Calibration/Human/',subjects(ii).sub_path,'Dif_fsl_preprocessed/eddy/aligned2T1/dtifit');
    if exist(MD_path)
        subjects(ii).MDtopup=1;
    else
        subjects(ii).MDtopup=0;
    end
    mprage_path=fullfile('/ems/elsc-labs/mezer-a/Mezer-Lab/rawData/HUJI/Calibration/Human/',subjects(ii).sub_path,'/nifti_dicm2nii/MPRAGE');
    if exist(mprage_path)
        subjects(ii).mprage=1;
    else
        subjects(ii).mprage=0;
    end
    
   if subjects(ii).age>55
        subjects(ii).group=['old ' num2str(subjects(ii).age)];
    else
        subjects(ii).group='young';
    end
%     if any(strcmp(old,subjects_path{ii}))
%         ind=find(strcmp(old,subjects_path{ii}));
%         subjects(ii).group=old_age{ind};
%     else
%         subjects(ii).group='young';
%     end
end

subjects(5).T2f=1;
run_on=1:length(subjects);
end