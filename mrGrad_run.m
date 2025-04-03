% MRGRAD RUN EXAMPLE
%
%   SOFTWARE REQUIREMENTS: 
%
%            * MATLAB          - http://www.mathworks.com/products/matlab/
%
%   Recommended:
%
%            * boundedline-pkg - https://github.com/kakearney/boundedline-pkg
%
% (C) Mezer lab, the Hebrew University of Jerusalem, Israel, Copyright 2021
%%   EXAMPLE DATA INPUT:
%   (see mrGrad.m documentation for more information)
%
%     Data{1}.group_name = 'Younger Adults';
%     Data{1}.map_list = {'path1','path2',...};
%     Data{1}.seg_list = {'path1','path2',...};
%     Data{1}.subject_names = {'sub1','sub2',...};
% 
%     Data{2}.group_name = 'Older Adults';
%     Data{2}.map_list = {'path1','path2',...};
%     Data{2}.seg_list = {'path1','path2',...};
%     Data{2}.subject_names = {'sub3','sub4',...};

% Alternatively, you can put your entire cohort together:
%     Data.group_name = 'my_datatset_name';
%     Data.map_list = {'path1','path2',...};
%     Data.seg_list = {'path1','path2',...};
%     Data.subject_names = {'sub1','sub2','sub3','sub4',...};

example_data_dir = fullfile(fileparts(matlab.desktop.editor.getActiveFilename),'example_data');
r1_map_path = fullfile(example_data_dir,'R1_map.nii.gz');
segmentation_path = fullfile(example_data_dir,'segmentation.nii.gz');


Data{1}.group_name = 'Younger Adults';
Data{1}.map_list = {r1_map_path};
Data{1}.seg_list = {segmentation_path};
Data{1}.subject_names = {'sub1'};

%%   EXAMPLE FUNCTION CALL:
% set analysis flags and options
roi_labels = [11,50,12,51];               % freesurfer labels of bilateral caudate and putamen (see RG/functions/freesurfer_LUT.txt)
segmenting_method = 'equidistance'; % equally-spaced segments
Axes = 1:3;                         % analysis axes
n_segments = [7,7,7];               % number of segments along each axis
stat = 'median';                    % statistic to obtain

% provide additional information
parameter = 'R1';                   % MRI parameter
parameter_units = '1/sec';          % units

% output folder
output_dir = fullfile(example_data_dir,'ExampleResults');
output_mode = 'extended'; % extended mode for saving result segmentation niftis

% run tool
RG = mrGrad(Data,'ROI',roi_labels,'n_segments',n_segments,'segmenting_method',segmenting_method,...
    'stat',stat,'PC',Axes,'erode',0,'invert',0,'param',parameter,'units',parameter_units,...
    'output_dir',output_dir,'output_mode',output_mode);

%% Visualize result spatial functions
% Display group-averaged gradients (and possibly individual subjects
% gradients) of 1 or more ROIs, in 1 or more subject-groups

group_ids = 1; % subject groups (rows of RG)
roi_ids = 1:4; % ROIs (columns of RG)
rg = RG(group_ids,roi_ids);
PC = 1:3;
ind = [0 0 0 0]; % if more than one subject, indicate for each group/roi whether to display individual gradients

y_lim = [.71 .88];
fig1 = showRG(rg,'errorType',2,'ind',ind,'PC',PC,...
    'markershapes',1,'legend',1,'labels',1,'ylim',y_lim);

fig1.WindowState = 'maximized';
%% Visualize example subjects segmentation along 1 axis
% Display the segmentation of 1 or more subjects along 1 axis of ROI(s)
% Colors indicate different segments

Axis = 1;      % region axis
group_id = 1;  % group number (rows of RG)
roi_ids = 1:4; % ROIs (columns of RG)
sub_ind = 1;   % subjects indices. To show all subjects: sub_ind = 1:length(RG{group_id,1}.individual_data);

fig2=region_axis_visualize(RG,group_id,roi_ids,sub_ind,Axis,...
    'pc_lines',1,'planes',0,'alpha',0.2);
fig2.WindowState = 'maximized';
%--------------------------------------------------------------------------