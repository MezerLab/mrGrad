function [maxchange_roi,msg] = get_roi_priors(roi)

msg = ['using priors for PCs directionality sign in ROI ',num2str(roi)];
roi_name = ROI_name(roi);
if contains(roi_name,'Caudate')
    maxchange_roi = [2 3 1]; % [y z x]: The longest projection of the 1st PC is on the Y-axis, etc.
elseif contains(roi_name,'Putamen')
    maxchange_roi = [2 3 1]; % [y z x]
elseif contains(roi_name,'Pallidum')
    maxchange_roi = [2 3 1];
elseif contains(roi_name,'Thalamus')
    maxchange_roi = [2 3 3];
else
    maxchange_roi = [2 3 2];
    msg = ['no default directionslity specs. for ROI ',num2str(roi),'. Agreement between subjects might be compromised'];
end