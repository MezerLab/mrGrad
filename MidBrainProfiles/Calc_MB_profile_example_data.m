% MidBrainProfile RUN EXAMPLE
%
%   SOFTWARE REQUIREMENTS: 
%
%            * MATLAB          - http://www.mathworks.com/products/matlab/
%               * Statistics and Machine Learning Toolbox
%               * Signal Processing Toolbox
%               * Symbolic Math Toolbox
%               * image processing toolbox
%            * Vistasoft       - https://github.com/vistalab/vistasoft   
%            * boundedline-pkg - https://github.com/kakearney/boundedline-pkg (recommended)
%
% (C) Mezer lab, the Hebrew University of Jerusalem, Israel, Copyright 2021


% addpath(genpath('..\vistasoft-master'))
% addpath(genpath('..\MBprofiles'))
 

%% first create a segmentation of the right and left sides of the midbrain
data_dir =  fullfile(fileparts(matlab.desktop.editor.getActiveFilename),'example_data');
brainstem_mask =  fullfile(data_dir, 'brainstem.nii.gz'); % this is the freesurfer outptut
output_dir = data_dir;       % this is where I would keep the resulting Midbrain 
MBfile = SplitMB(brainstem_mask,output_dir);

%% Now, we can create the profiles along those mask
%% Data input
r1_map_path = fullfile(data_dir,'R1_map.nii.gz');

Data{1}.group_name = 'Younger Adults';
Data{1}.map_list = {r1_map_path};
Data{1}.seg_list = {MBfile};
Data{1}.subject_names = {'sub1'};

% % % % % % % % % for nalyzing several subjects groups the following input
% % % % % % % % % can be provided to the function
 
%     Data{1}.group_name = 'Younger Adults';
%     Data{1}.map_list = {'path1','path2',...}; % to R1 for examples
%     Data{1}.seg_list = {'path1','path2',...}; % to the MB masks
%     Data{1}.subject_names = {'sub1','sub2',...};
% 
%     Data{2}.group_name = 'Older Adults';
%     Data{2}.map_list = {'path1','path2',...};
%     Data{2}.seg_list = {'path1','path2',...};
%     Data{2}.subject_names = {'sub3','sub4',...};

%% Function call
 
% set analysis flags and options
rois = [1,2];               % freesurfer labels of bilateral caudate and putamen (see RG/functions/freesurfer_LUT.txt)
segmentation_method = 'spacing';    % equally-spaced segments
Axes = 1:3;                         % analysis axes
Ns = [7,7,7];                       % number of segments along each axis
stat = 'mean';                    % statistic to obtain. % for sampling masks (and not maps) this should be mean (and not median)

% provide additional information
parameter = 'R1';                   % MRI parameter
parameter_units = '1/sec';          % units

% run tool
RG = mrGrad(Data,'ROI',rois,'Nsegs',Ns,'roi_names',{'MB_L','MB_R'},'segmentingMethod',segmentation_method,...
    'stat',stat,'PC',Axes,'erode',1,'invert',0,'param',parameter,'units',parameter_units,'max_change', [2 3 1]);

% % SAVE OUTPUT
% save(FileNameOut,'RG');
%% Visualize result spatial functions
% Display group-averaged gradients (and possibly individual subjects
% gradients) of 1 or more ROIs, in 1 or more subject-groups

group_ids = 1; % subject groups (rows of RG)
roi_ids = 1:2; % ROIs (columns of RG)
rg = RG(group_ids,roi_ids);
PC = 1:3;
ind = [0 0 0 0]; % if more than one subject, indicate for each group/roi whether to display individual gradients

y_lim = [];%[.71 .88];
fig1 = showRG(rg,'errorType',2,'ind',ind,'PC',PC,...
    'markershapes',1,'legend',1,'labels',1,'ylim',y_lim);

fig1.WindowState = 'maximized';
%% Visualize example subjects segmentation along 1 axis
% Display the segmentation of 1 or more subjects along 1 axis of ROI(s)
% Colors indicate different segments

Axis = 1; % region axis
group_id = 1; % group number (rows of RG)
roi_ids = 1:2; % ROIs (columns of RG)
sub_ind = 1;%1:length(RG{group_id,1}.individual_data); % show all subjects

fig2=region_axis_visualize(RG,group_id,roi_ids,sub_ind,Axis,...
    'pc_lines',1,'planes',0,'alpha',0.1);
fig2.WindowState = 'maximized';
%--------------------------------------------------------------------------