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
Path = mfilename('fullpath');
lut_path=fileparts(Path);
LUT = readtable(fullfile(lut_path,'freesurfer_LUT.txt'));

%--------------------------------------------------------------------------
if ~exist('r','var') || isempty(r)
    Out = LUT;
    return
elseif isnumeric(r)
    % CASE: input is ROI label; RETURN: ROI name
    name = arrayfun(@(x) LUT.Name(LUT.Label==x),r,'un',0);
    unknown = cellfun(@isempty,name);
    name(unknown) = {'Unknown'};
    name = cat(1,name{:});
    Out = name;
else
    % CASE: input is ROI name; RETURN: ROI label
    r = string(r);
    label = arrayfun(@(x) LUT.Label(LUT.Name==x),r,'un',0);
    unknown = cellfun(@isempty,label);
    label(unknown) = {nan};
    label = cat(1,label{:});
    Out = label;
%--------------------------------------------------------------------------    
end
