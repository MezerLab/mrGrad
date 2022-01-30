% MRGRAD RUN EXAMPLE
%
%   SOFTWARE REQUIREMENTS: 
%
%            * MATLAB          - http://www.mathworks.com/products/matlab/
%               * Statistics and Machine Learning Toolbox
%               * Signal Processing Toolbox
%               * Symbolic Math Toolbox
%            * Vistasoft       - https://github.com/vistalab/vistasoft   
%            * boundedline-pkg - https://github.com/kakearney/boundedline-pkg (recommended)
%
% (C) Mezer lab, the Hebrew University of Jerusalem, Israel, Copyright 2021
%%   EXAMPLE DATA INPUT:
%   (see mrGrad.m documentation for more information)

%     Data{1}.group_name = 'Younger Adults';
%     Data{1}.map_list = {'path1','path2',...};
%     Data{1}.seg_list = {'path1','path2',...};
%     Data{1}.subject_names = {'sub1','sub2',...};
% 
%     Data{2}.group_name = 'Older Adults';
%     Data{2}.map_list = {'path1','path2',...};
%     Data{2}.seg_list = {'path1','path2',...};
%     Data{2}.subject_names = {'sub3','sub4',...};

example_data_dir = fullfile(fileparts(matlab.desktop.editor.getActiveFilename),'example_data');
r1_map_path = fullfile(example_data_dir,'R1_map.nii.gz');
segmentation_path = fullfile(example_data_dir,'segmentation.nii.gz');


Data{1}.group_name = 'Younger Adults';
Data{1}.map_list = {r1_map_path};
Data{1}.seg_list = {segmentation_path};
Data{1}.subject_names = {'sub1'};

%%   EXAMPLE FUNCTION CALL:

% set analysis flags and options
rois = [11,50,12,51];               % freesurfer labels of bilateral caudate and putamen (see RG/functions/freesurfer_LUT.txt)
segmentation_method = 'spacing';    % equally-spaced segments
Axes = 1:3;                         % analysis axes
Ns = [7,7,7];                       % number of segments along each axis
stat = 'median';                    % statistic to obtain

% provide additional information
parameter = 'R1';                   % MRI parameter
parameter_units = '1/sec';          % units

% run tool
RG = mrGrad(Data,'ROI',rois,'Nsegs',Ns,'segmentingMethod',segmentation_method,...
    'stat',stat,'PC',Axes,'erode',0,'invert',0,'param',parameter,'units',parameter_units);

% % SAVE OUTPUT
% save(FileNameOut,'RG');
%% Visualize result spatial functions
% Display group-averaged gradients (and possibly individual subjects
% gradients) of 1 or more ROIs, in 1 or more subject-groups

group_ids = 1; % subject groups (rows of RG)
roi_ids = 1:4; % ROIs (columns of RG)
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
roi_ids = 1:4; % ROIs (columns of RG)
sub_ind = 1%1:length(RG{group_id,1}.individual_data); % show all subjects

fig2=region_axis_visualize(RG,group_id,roi_ids,sub_ind,Axis,...
    'pc_lines',1,'planes',0,'alpha',0.1);
fig2.WindowState = 'maximized';
%--------------------------------------------------------------------------