function [MBfile] = SplitMB(brainstem_mask,output_dir)
% this function is a semi-automatic aid for finding the midsagittal point
% of the midbrain.

% brainstem_mask is the output file of th efreesyrfer brainstem
% segmentation. 
% output_dir is where the new files will be saved
%% test the midbrain seg is there
brainstem = readFileNifti(brainstem_mask);
if isempty(brainstem.data)
    error(['brain stem file does not exist or the path is incorrect'])    
end

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end
%% 1. save the MB from the brainstem
MBmask_file = fullfile(output_dir,'MBmask.nii.gz'); 
if ~exist(MBmask_file,'file')
    
    MBval = 173;
    MBmask = single(brainstem.data == MBval);
    dtiWriteNiftiWrapper(MBmask, brainstem.qto_xyz,MBmask_file);
else
    
    MBold = readFileNifti(MBmask_file);
    MBmask = MBold.data;
    
end
    
%% 2. find mid sagittal x coordinate

MBfile = fullfile(output_dir,'MB.nii.gz');

if ~exist(MBfile,'file')
    % find the general middle to ssee the MB well
     mbinds = find(MBmask);
    [x,y,z] = ind2sub(size(MBmask),mbinds);
    zslice = median(z); 
    
    % open figure and choose the middle
    f=figure( 'position',[680,25,  1235,1050]);
    imshow(MBmask(:,:,zslice));%imagesc(MBold.data(:,:,zslice))
    warning('if you are using a keyboard and not a mouse, and the window doesnot close after pressing, change "button" to 2')
    
    button=2;
    while button==2
        figure(f),title('Choose midsagittal line')
        [xx, yy, button] = ginput(1);
    end
    
    xxMidSag = cast(floor(yy), 'int16');  % size(MBold.data,1)-
    close(f)

    % save the right and left MB
    
    MB = zeros(size(MBmask));
       
    leftInds = find(x>=xxMidSag);
    MB(mbinds(leftInds ))=2;
    rightInds = find(x<xxMidSag);
    MB(mbinds(rightInds))=1;
    
    dtiWriteNiftiWrapper(single(MB),brainstem.qto_xyz,MBfile);
     
end
