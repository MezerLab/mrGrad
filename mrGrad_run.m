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

example_data_dir = fullfile(fileparts(matlab.desktop.editor.getActiveFilename),'example_data');
r1_map_path = fullfile(example_data_dir,'R1_map.nii.gz');
segmentation_path = fullfile(example_data_dir,'segmentation.nii.gz');

Data.seg_list = {segmentation_path};
Data.map_list = {r1_map_path};
Data.subject_ids = {'sub-1'};

%%   EXAMPLE FUNCTION CALL:
% set analysis flags and options
roi_labels = [11,50,12,51];         % freesurfer labels of bilateral caudate and putamen (see RG/functions/freesurfer_LUT.txt)
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

% run mrGrad
RG = mrGrad(Data,'ROI',roi_labels,'n_segments',n_segments,'segmenting_method',segmenting_method,...
    'stat',stat,'PC',Axes,'erode',0,'parameter_names',parameter,'units',parameter_units,...
    'output_dir',output_dir,'output_mode',output_mode, 'Parallel', false);

%% Visualize result spatial functions
% Display group-averaged gradients (and possibly individual subjects
% gradients) of 1 or more ROIs, in 1 or more subject-groups


mrgrad_show_gradients(RG,'error_name','SEM');

% Alternatively, select regions to show:
rg = rmfield(RG,{'Left_Caudate', 'Right_Caudate', 'Right_Putamen'});
mrgrad_show_gradients(rg,'error_name','STD');


%% Visualize example subjects segmentation along 1 axis
% Display the segmentation of 1 or more subjects along 1 axis of ROI(s)
% Colors indicate different segments

axis_name = 'axis1'; % region axis to display

subject_idx = 1;   % subjects indices. To show all subjects: sub_ind = 1:length(RG{group_id,1}.individual_data);

fig2=mrgrad_axis_visualize(RG, axis_name, subject_idx,...
    'pc_lines',1,'planes',0,'alpha',0.2,'ax_view',[-45, 45],'cmap','autumn');
%--------------------------------------------------------------------------