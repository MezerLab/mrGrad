function [Out, LUT] = ROI_name(r)
%--------------------------------------------------------------------------
% INPUT
%         ROI label(s)    e.g. 10
%           OR
%         ROI name(s)     e.g. 'Left-Thalamus'
% 
%         (In case no argument provided - function outputs all LUT)
%
% OUTPUT
%
%         Out:    ROI name(s)     e.g. 'Left-Thalamus'
%                       OR
%                 ROI label(s)    e.g. 10
% 
%         LUT:    A shortened freesurfer lookup table
% 
% (C) Elior Drori, Mezer lab, the Hebrew University of Jerusalem, Israel, Copyright 2021
%--------------------------------------------------------------------------

% load lookup table
path = mfilename('fullpath');
a=fileparts(path);
LUT = readtable(fullfile(a,'freesurfer_LUT.txt'));

All_labels = table2array(LUT(:,1));
All_names = table2cell(LUT(:,2));
%--------------------------------------------------------------------------
% CASE: input is ROI label; RETURN: ROI name
if exist('r','var') && isnumeric(r)
    name = cell(length(r),1);
    for j = 1:length(r)
        ind = ismember(All_labels,r(j));
        if any(ind)
            name{j} = All_names{ind};
        else
            name{j} = 'Unknown';
        end
    end
    Out = name;
%--------------------------------------------------------------------------    
% CASE: input is ROI name; RETURN: ROI label
elseif exist('r','var') && (isa(r,'cell') || isa(r,'char'))
    if isa(r,'char')
        r = {r};
    end
    label = zeros(length(r),1);
    for j = 1:length(r)
        ind = ismember(All_names,r(j));
        if any(ind)
            label(j) = All_labels(ind);
        else
            label(j) = NaN;
        end            
    end
    Out = label;
%--------------------------------------------------------------------------    
% CASE: no argument provided; RETURN: whole list labels + names
elseif ~exist('r','var')
    Out = LUT;
else
    Out = [];
end
end