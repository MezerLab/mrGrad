function RG_lr_func = mrGrad_operation_LR(RG,FunctionHandle,operation_name,hemisphere_notations)
% This is a generalized function of "mrGrad_average_LR.m", which compute
% any given function on L, R mrGrad matrices

% FunctionHandle    A function handle with two mtrices arguments l, r
%                   Example: FunctionHandle = @(l,r) (l+r)/2;

if ~exist('operation_name','var') || isempty(operation_name)
    operation_name = 'computation';
end

if ~exist('hemisphere_notations','var')
    hemisphere_notations = {'Left-','Right-'};
end
roi_labels = cellfun(@(x) x.ROI_label,RG(:),'un',0);
region_names = unique(cellfun(@(x) erase(x,hemisphere_notations),roi_labels,'un',0),'stable')';

ngroups = size(RG,1);
nregions = length(region_names);

RG_lr_func = cell(ngroups,nregions);
for gg = 1:ngroups
    for rr = 1:nregions
        rg = RG(gg,:);
        rg(cellfun(@(x) ~contains(x.ROI_label,region_names{rr}),rg)) = [];
        if numel(rg)~=2
            error('number of hemispheres per regions are different than 2')
        end
        rg_l = rg{cellfun(@(x) contains(x.ROI_label,hemisphere_notations{1}),rg)};
        rg_r = rg{cellfun(@(x) contains(x.ROI_label,hemisphere_notations{2}),rg)};

        RG_lr_func{gg,rr}.Y = cellfun(@(a,b) FunctionHandle(a,b), rg_l.Y, rg_r.Y,'un',0);
        RG_lr_func{gg,rr}.Y_mean = cellfun(@(x) mean(x,2,"omitnan"),RG_lr_func{gg,rr}.Y,'un',0);

        RG_lr_func{gg,rr}.Y_std = cellfun(@(x) std(x,0,2,"omitnan"),RG_lr_func{gg,rr}.Y,'un',0);
        RG_lr_func{gg,rr}.Y_SEM = cellfun(@(x) std(x,0,2,"omitnan")/sqrt(size(x,2)),RG_lr_func{gg,rr}.Y,'un',0);

        additional_fields = fieldnames(rg{1});
        additional_fields(ismember(additional_fields,fieldnames(RG_lr_func{gg,rr}))) = [];
        for ff = 1:length(additional_fields)
            RG_lr_func{gg,rr}.(additional_fields{ff}) = rg{1}.(additional_fields{ff});
        end
        RG_lr_func{gg,rr}.ROI_label = sprintf('LR-%s-%s',operation_name,region_names{rr});
        
        if isfield(RG_lr_func{gg,rr},'individual_data')
            RG_lr_func{gg,rr}.WARNING_NOTE = sprintf('"individual_data" field is of %s side only',rg{1}.ROI_label);
        end
    end
end


