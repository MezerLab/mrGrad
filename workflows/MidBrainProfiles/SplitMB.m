function [MBfile] = SplitMB(brainstem_mask,output_dir)
% this function is a semi-automatic aid for finding the midsagittal point
% of the midbrain.

% brainstem_mask is the output file of th efreesyrfer brainstem
% segmentation. 
% output_dir is where the new files will be saved
%% test the midbrain seg is there
if exist(brainstem_mask,'file')
    brainstem = niftiread(brainstem_mask);
else
    error('brain stem file does not exist or the path is incorrect');
end

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end
%% 1. save the MB from the brainstem
MBval = 173;
MBmask = single(brainstem == MBval);
info = niftiinfo(brainstem_mask);
info.Datatype = 'single';
MBmask_file = fullfile(output_dir,'MBmask.nii.gz');
niftiwrite(MBmask,regexprep(MBmask_file,'.gz$',''),info,Compressed=true)
    
%% 2. find mid sagittal x coordinate

MBfile = fullfile(output_dir,'MB.nii.gz');

% find the general middle to ssee the MB well
mbinds = find(MBmask);
[x,y,z] = ind2sub(size(MBmask),mbinds);
zslice = median(z); 
    
% open figure and choose the middle
fig=figure(WindowState="maximized");
imshow(MBmask(:,:,zslice));
warning('if you are using a keyboard and not a mouse, and the window doesnot close after pressing, change "button" to 2')
    
button=2;
while button==2
    figure(fig),title('USER ACTION REQUIRED: Choose midsagittal line:',Color='r',FontSize=20);
    [xx, yy, button] = ginput(1);
end

xxMidSag = cast(floor(yy), 'int16');  % size(MBold.data,1)-
close(fig)

% save the right and left MB
MB = single(zeros(size(MBmask)));
   
leftInds = find(x>=xxMidSag);
MB(mbinds(leftInds)) = 2;
rightInds = find(x<xxMidSag);
MB(mbinds(rightInds)) = 1;

niftiwrite(MB, regexprep(MBfile,'.gz$',''), info, Compressed=true);

fprintf('Midbrain (left, right) segmentation saved in: %s\n',MBfile);
