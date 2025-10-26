% MidBrainProfile RUN DEMO
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
example_data_dir =  fullfile(fileparts(matlab.desktop.editor.getActiveFilename),'example_data');
brainstem_mask =  fullfile(example_data_dir, 'brainstem.nii.gz'); % this is the freesurfer outptut
output_dir = fullfile(example_data_dir,'ExampleResults'); % this is where I would keep the resulting Midbrain 
MBfile = SplitMB(brainstem_mask,output_dir);

%% Now, we can create the profiles along those mask
%% Data input
r1_map_path = fullfile(example_data_dir,'R1_map.nii.gz');

Data.seg_list = {MBfile};
Data.map_list = {r1_map_path};
Data.subject_names = {'sub1'};

%% Function call
 
% set analysis flags and options
rois = [1,2];               % freesurfer labels of bilateral caudate and putamen (see RG/functions/freesurfer_LUT.txt)
segmentation_method = 'spacing';    % equally-spaced segments
Axes = 1:3;                         % analysis axes
n_segments = [7,7,7];                       % number of segments along each axis
stat = 'mean';                    % statistic to obtain. % for sampling masks (and not maps) this should be mean (and not median)

% provide additional information
parameter = 'R1';                   % MRI parameter
parameter_units = '1/sec';          % units

% run tool
RG = mrGrad(Data,'ROI',rois,'n_segments',n_segments,'roi_names',{'MB_L','MB_R'},'segmenting_method',segmentation_method,...
    'stat',stat,'PC',Axes,'erode',1,'parameter_names',parameter,'units',parameter_units,'max_change', [2 3 1],...
    'output_dir',output_dir);

%% Visualize result spatial functions
mrgrad_show_gradients(RG,'error_name','SEM');

%% Visualize example subjects segmentation along 1 axis
% Display the segmentation of 1 or more subjects along 1 axis of ROI(s)
% Colors indicate different segments

axis_name = 'axis1'; % region axis to display
mrgrad_axis_visualize(RG, axis_name, [],...
    'pc_lines',1,'planes',0,'alpha',0.2,'ax_view',[-30, 75],'cmap','autumn');
%--------------------------------------------------------------------------