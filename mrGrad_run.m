% MRGRAD V2.0 RUN EXAMPLE
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
%   Given N subject sessions and M MRI parameters (t1, t2,...), the input
%   should be a DATA struct:
%
%     Data.seg_list (N x 1) = {'seg_path1'; 'seg_path2'; ...};
%     Data.map_list  (N x M) = {'t1_path_1', 't2_path_1'; 't1_path_2', 't2_path_2',...};
%     Data.subject_ids (N x 1) = {'sub1-ses1'; 'sub1-ses2'; 'sub2'; ...};
clear;
example_data_dir = fullfile(fileparts(matlab.desktop.editor.getActiveFilename),'example_data');

subject_ids = {'sub-1'; 'sub-2'; 'sub-3'};
segmentation_path = fullfile(example_data_dir,subject_ids,'segmentation.nii.gz');
r1_path = fullfile(example_data_dir,subject_ids,'R1_map.nii.gz');
wf_path = fullfile(example_data_dir,subject_ids,'WF_map.nii.gz');
mtsat_path = fullfile(example_data_dir,subject_ids,'MTsat_map.nii.gz');

Data.seg_list = segmentation_path;
Data.map_list = [r1_path, wf_path, mtsat_path];
Data.subject_ids = subject_ids;

%%   EXAMPLE FUNCTION CALL:
% set analysis flags and options
roi_labels = [11,50,12,51];         % freesurfer labels of bilateral caudate and putamen (see RG/functions/freesurfer_LUT.txt)
segmenting_method = 'equidistance'; % equally-spaced segments
Axes = 1:3;                         % analysis axes
n_segments = [7,7,7];               % number of segments along each axis
stat = 'median';                    % statistic to obtain

% provide additional information
parameter_names = {'R1', 'WaterFraction', 'MTsaturation'}; % MRI parameter
parameter_units = {'1/s', 'fraction', 'p.u.'}; % units

% output folder
output_dir = fullfile(example_data_dir,'ExampleResults');
output_mode = 'default'; % extended mode for saving result segmentation niftis

% run mrGrad
RG = mrGrad(Data,'ROI',roi_labels,'n_segments',n_segments,'segmenting_method',segmenting_method,...
    'stat',stat,'PC',Axes,'erode',0,'parameter_names',parameter_names,'units',parameter_units,...
    'output_dir',output_dir,'output_mode',output_mode, 'Parallel', false);

%% Visualize result spatial functions
% Display group-averaged gradients (and possibly individual subjects
% gradients) of 1 or more ROIs, in 1 or more subject-groups

mrgrad_show_gradients(RG,'error_name','SEM');

%% Alternatively, show only left putamen, MTsat, group by sex:
rg = rmfield(RG,{'Left_Caudate', 'Right_Caudate', 'Right_Putamen'});

group_labels = categorical(["M", "F", "M"]);
mrgrad_show_gradients(rg,'error_name','STD', 'group_labels', group_labels,...
    'parameter_names','MTsaturation');

%% Alternatively, show mean(left,right) ROIs
RG_avg = mrGrad_average_LR(RG);
mrgrad_show_gradients(RG_avg,'error_name','SEM');

%% Altrnatively, remove subjects
subIdx = [1 0 1];
RG_sub = mrGrad_subset(RG, subIdx);
mrgrad_show_gradients(RG_sub,'error_name','SEM');

%% Visualize example subjects segmentation along 1 axis
% Display the segmentation of 1 or more subjects along 1 axis of ROI(s)
% Colors indicate different segments

axis_name = 'axis1'; % region axis to display
subIdx = true(1,3);   % subjects indices to show (up to 30)

fig2=mrgrad_axis_visualize(RG, axis_name, subIdx,...
    'pc_lines',1,'planes',0,'alpha',0.2,'ax_view',[-45, 45],'cmap','autumn');
%--------------------------------------------------------------------------