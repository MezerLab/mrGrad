function singlsb_rgs = mrgrad_inner_loop(ii, rr, seg_list, map_list, mrgrad_defs, maxchange_roi)
%             fprintf('%d\n',ii); % uncomment for debugging
            %----------------------------------------------------------------------
            % load subject's qMRI data
            %----------------------------------------------------------------------
            roi = mrgrad_defs.ROI(rr);

            stat = mrgrad_defs.stat;
            PC = mrgrad_defs.PC;
            n_segments = mrgrad_defs.n_segments;
            sampling_method = mrgrad_defs.segmentingMethod;
            isfigs = mrgrad_defs.isfigs;

            if ~exist(seg_list{ii},'file')
                singlsb_rgs = [];
                return
            end

            images_exist = logical(cellfun(@(x) exist(x,'file'), map_list(ii,:)));
            if all(~images_exist)
                singlsb_rgs = [];
                return
            end

            % load segmentation and keep nifti info
            mask = ROImask(seg_list{ii},roi,mrgrad_defs.erode_flag);
            image_info = niftiinfo(seg_list{ii});
            [strides,im_dims] = keep_strides(image_info);

            images = cellfun(@(im) double(niftiread(im)), map_list(ii,images_exist),'un',0);
            imsize = cellfun(@size,images,'un',0);

            sz = size(mask);

            % Make sure mask and image have the same dimensions
            if any(cellfun(@(x) ~isequal(sz,x), imsize))
                error('Input image and mask/segmentation'' dimensions must agree.');
            end
            %----------------------------------------------------------------------
            % Outliers removal
            %----------------------------------------------------------------------
            % warn about outliers
            remove_vox_outliers = false;
            if remove_vox_outliers
%                 out_idx = cellfun(@(im) isoutlier(im(mask>0)), images, 'un',0);
%                 out_idx = cat(2,out_idx{:});
% 
%                 mask(mask>0) = ~Outlier;
%                 warning('%d outliers removed.',nnz(Outlier));
            end
            %----------------------------------------------------------------------
            % mask the image
            %----------------------------------------------------------------------
            images = cellfun(@(im) (mask.* im),images,'un',0);

            %----------------------------------------------------------------------
            % if an alternative ROI was given for Axes calculation, make mask
            %----------------------------------------------------------------------
            alternative_mask = [];
            if isfield(mrgrad_defs,'Alternative_ROI')
                alternative_roi = mrgrad_defs.Alternative_ROI(rr);
                alternative_seg = mrgrad_defs.Alternative_seg_list{ii};
                alternative_mask = ROImask(alternative_seg,alternative_roi,0);
            end
            %----------------------------------------------------------------------
            % make sure data is in positive strides (L>R P>A I>S)
            %----------------------------------------------------------------------

            % change image strides order to [1,2,3]
            [~, im_perm] = sort(im_dims);
            images = cellfun(@(im) permute(im,im_perm), images,'un',0);
            mask = permute(mask,im_perm);
            alternative_mask = permute(alternative_mask,im_perm);

            % flip negative strides to achieve [+1,+2,+3]
            dimsflip = im_dims(strides < 0);
            for d = dimsflip
                images = cellfun(@(im) flip(im,d), images,'un',0);
                mask = flip(mask,d);
                alternative_mask = flip(alternative_mask,d);
            end
            % Warn once about strides' change
            if ~isequal(strides,[1,2,3]) && ...
                    (ii==1 || ~mrgrad_defs.parallel)
                warning('mrGrad:Strides','\nImages of some/all subjects are flipped to match positive strides.')
                warning('off','mrGrad:Strides');
            end

            %----------------------------------------------------------------------
            % MAIN FUNCTION CALL mrgrad_per_sub.m
            %----------------------------------------------------------------------
            % single subject mrgrads in (up to) 3 PCs
            singlsb_rgs = arrayfun(@(x,y)...
                mrgrad_per_sub(images,mask,'PC',x,'n_segments',y,'sampling_method',sampling_method,...
                'stat',stat,'maxchange',maxchange_roi,...
                'subID',ii,'isfigs',isfigs,'apply_alternative_axes',alternative_mask),...
                PC,n_segments);

            % add nan values for parameters that do not exist
            for pc = 1:length(mrgrad_defs.PC)
                function_data = double(nan(mrgrad_defs.n_segments(pc),length(map_list(ii,:))));
                function_data(:,images_exist) = singlsb_rgs(pc).function;
                singlsb_rgs(pc).function = function_data;
            end

            if mrgrad_defs.output_mode=="minimal"
                singlsb_rgs = rmfield(singlsb_rgs, setdiff(fieldnames(singlsb_rgs), 'function'));
            else
                % Keep original stride info for generating segmentation files
                [singlsb_rgs.analysis_image_size] = deal(size(mask));
                [singlsb_rgs.original_image_size] = deal(sz);
                singlsb_rgs(1).original_nifti_info = image_info;
                
                [singlsb_rgs.analysis_strides] = deal([1,2,3]);
                [singlsb_rgs.original_strides] = deal(strides);
                [singlsb_rgs.parameter_names] = deal(mrgrad_defs.parameter_names(:)');
            end
            %-------------------------------------
        end