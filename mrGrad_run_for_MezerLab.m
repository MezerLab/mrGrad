% MRGRAD RUN HUJI
clearvars -except gcp;
dataset = 'HUJI'; % See generate_rg_inputs.m for details
param = 'R1'; % See huji_MapAndSeg.m for details
units = '';%mri_units(param);
Data = generate_rg_inputs('dataset',dataset,'param',param);

rois = [11,50,12,51];
segmentation_method = 'spacing';
Axes = 1:3;
Ns = [7,7,7];
stat = 'median';

gcp();

RG = mrGrad_parallel(Data,'ROI',rois,'Nsegs',Ns,'segmentingMethod',segmentation_method,...
    'stat',stat,'PC',Axes,'erode',0,'invert',0,'param',param,'units',units);

% SAVE OUTPUT
outDir = '/tmp';
outFile = fullfile(outDir,sprintf('RG_%s_%s.mat',dataset,param));
save(outFile,'RG');
%% Visualize result spatial functions

rg_dir = outDir;
rg_path = fullfile(rg_dir,sprintf('RG_%s_%s.mat',dataset,param));

RG = load(rg_path);
RG = RG.RG;

group_ids = 1:2; % subject groups (rows of RG)
roi_ids = 3; % ROIs (columns of RG)
rg = RG(group_ids,roi_ids);
PC = 1:3;
ind = [0 0 0 0];

y_lim = [];%[.71 .88];
fig1 = showRG(rg,'errorType',2,'ind',ind,'PC',PC,...
    'markershapes',1,'legend',1,'labels',1,'ylim',y_lim,'gradnumber',0);

fig1.WindowState = 'maximized';
%% Visualize example subjects segmentation along 1 axis
cls;
Axis = 3;
group_id = 1;
roi_ids = 1:4;
sub_ind = 1;%1:length(RG{group_id,1}.individual_data); % show all subjects

fig2=region_axis_visualize(RG,group_id,roi_ids,sub_ind,Axis,...
    'pc_lines',1,'planes',0,'centroids',1,'alpha',0.1);
fig2.WindowState = 'maximized';
% view([ 0 90]) % xy
% view([-90 0]) % yz
% view([  0 0]) % xz
%--------------------------------------------------------------------------