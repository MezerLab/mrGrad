function RG_avg = mrGrad_average_LR(RG, hemisphere_notations)
% mrGrad_average_LR
%
% DESCRIPTION:
%   Computes the left–right average of ROI-based mrGrad structures.
%   Each hemisphere-specific ROI (e.g. 'Left_Caudate', 'Right_Caudate')
%   is averaged to create a new ROI entry in RG_avg with the mean of 
%   corresponding numeric table values in the field `.Results.(subfield).y`.
%
% SYNTAX:
%   RG_avg = mrGrad_average_LR(RG)
%   RG_avg = mrGrad_average_LR(RG, hemisphere_notations)
%
% INPUTS:
%   RG : struct
%       mrGrad-style struct containing hemisphere-specific ROI fields.
%       Each RG.(ROI).Results.(subfield).y must be a table.
%
%   hemisphere_notations : cell array of char (optional)
%       Prefixes defining hemisphere naming, e.g. {'Left','Right'}, {'lh','rh'}, or {'_L', '_R'}.
%       Defaults to {'Left','Right'}.
%
% OUTPUT:
%   RG_avg : struct
%       Struct with averaged ROIs (e.g. 'Caudate', 'Putamen'), each containing
%       averaged `.Results.(subfield).y` tables and same metadata.
%
% NOTES:
%   - Assumes left/right tables are identical in size, row names, and variable names.
%   - Removes the field 'individual_data' (cannot be averaged).
%   - Skips ROIs missing a hemisphere counterpart.

% -------------------------------------------------------------------------
% Input Validation
% -------------------------------------------------------------------------
if ~isa(RG,'struct')
    error('Input RG must be a struct.');
end
rg_names = fieldnames(RG)';

if ~exist('hemisphere_notations','var') || isempty(hemisphere_notations)
    hemisphere_notations = {'Left','Right'};
else
    hemisphere_notations = cellstr(hemisphere_notations);
    if numel(hemisphere_notations)~=2
        error('hemisphere_notations must contain exactly two string elements, e.g. {''Left'',''Right''}.');
    end
end

warn = true;

% -------------------------------------------------------------------------
% Derive base ROI names (remove hemisphere notations)
% -------------------------------------------------------------------------
base_names = cellfun(@(f) erase(f, hemisphere_notations), rg_names, 'un', 0);
base_names = regexprep(unique(base_names,'stable'), '^[-_]+|[-_]+$', ''); % Remove any leading or trailing underscores or dashes from each string 

RG_avg = struct;

% -------------------------------------------------------------------------
% Loop over each base ROI
% -------------------------------------------------------------------------
for ii = 1:numel(base_names)
    base = base_names{ii};

    % Identify hemisphere-specific fieldnames
    left_field  = rg_names{contains(rg_names, base) & contains(rg_names, hemisphere_notations{1})};
    right_field  = rg_names{contains(rg_names, base) & contains(rg_names, hemisphere_notations{2})};

    % Skip if either hemisphere missing
    if isempty(left_field) || isempty(right_field)
        warning('Skipping %s: missing hemisphere field.', base);
        continue;
    end

    % Initialize averaged struct from left (copy metadata)
    RG_avg.(base) = RG.(left_field);
    RG_avg.(base).roi_name = ['LR-average-',base];

    % Remove individual-level axes data (cannot be averaged)
    if isfield(RG_avg.(base),'individual_data')
        RG_avg.(base) = rmfield(RG_avg.(base), 'individual_data');
        if warn
            warning('Removed ''individual_data'' field (non-averagable).');
            warn = false;
        end
    end

    % ---------------------------------------------------------------------
    % Average Results subfields
    % ---------------------------------------------------------------------
    subfields = fieldnames(RG.(left_field).Results);

    for j = 1:numel(subfields)
        subf = subfields{j};

        yL = RG.(left_field).Results.(subf).y;
        yR = RG.(right_field).Results.(subf).y;

        % Validate table compatibility
        if ~isequal(size(yL), size(yR)) || ~isequal(yL.Properties.VariableNames, yR.Properties.VariableNames)
            warning('Skipping %s.%s: table mismatch between hemispheres.', base, subf);
            continue;
        end

        % Compute LR average
        RG_avg.(base).Results.(subf).y{:,:} = (yL{:,:} + yR{:,:}) / 2;
        
        % Re-compute mean, std, sem
        RG_avg.(base).Results.(subf).y_mean = mean(RG_avg.(base).Results.(subf).y{:,:},1,"omitnan");
        RG_avg.(base).Results.(subf).y_std = std(RG_avg.(base).Results.(subf).y{:,:},[],1,"omitnan");
        n_actual = sum(~isnan(RG_avg.(base).Results.(subf).y{:,:}),1);
        RG_avg.(base).Results.(subf).y_sem = RG_avg.(base).Results.(subf).y_std ./ sqrt(n_actual);

    end
end