function RG_sub = mrGrad_subset(RG, subIdx)
% mrGrad_subset
%
% DESCRIPTION:
%   Creates a subset of an mrGrad structure by selecting subjects
%   according to a logical index vector (subIdx).
%   All .Results.(subfield).y tables, .subject_ids, and .individual_data
%   fields are subset to match subIdx.
%
% SYNTAX:
%   RG_sub = mrGrad_subset(RG, subIdx)
%
% INPUTS:
%   RG : struct
%       mrGrad-style struct
%
%   subIdx : logical vector
%       Index of subjects to retain. Must have length equal to the
%       number of subjects in RG.
%
% OUTPUT:
%   RG_sub : struct
%       Same as RG, but only with selected subjects retained.
%
% NOTES:
%   - Validates consistency across all ROIs and subfields.
% -------------------------------------------------------------------------
% Input validation
% -------------------------------------------------------------------------
if ~isstruct(RG)
    error('Input RG must be a struct.');
end

if nargin < 2 || isempty(subIdx)
    error('subIdx must be provided (logical vector).');
end

subIdx = logical(subIdx);
RG_sub = RG;
% -------------------------------------------------------------------------
% Sanity checks
% -------------------------------------------------------------------------
if all(subIdx)
    warning('mrGrad_subset:AllSelected', ...
        'All subjects are selected — nothing to subset. Returning input RG unchanged.');
    return;
end

if ~any(subIdx)
    error('mrGrad_subset:NoSubjectsSelected', ...
        'subIdx selects no subjects — resulting dataset would be empty.');
end

warn = true;
% -------------------------------------------------------------------------
% Validate and subset
% -------------------------------------------------------------------------
roi_names = fieldnames(RG);
for ii = 1:numel(roi_names)
    roi = roi_names{ii};
    subfields = fieldnames(RG.(roi).Results);

    % --- Subset all y tables ---
    for j = 1:numel(subfields)
        subf = subfields{j};
        y = RG.(roi).Results.(subf).y;
        if ~isequal(height(y), numel(subIdx))
            error('ROI "%s", field "%s": subIdx length (%d) does not match table height (%d).', ...
                roi, subf, numel(subIdx), height(y));
        end
        y(~subIdx, :) = [];  % remove excluded rows
        RG_sub.(roi).Results.(subf).y = y;

        % Re-compute mean, std, sem
        RG_sub.(roi).Results.(subf).y_mean = mean(RG_sub.(roi).Results.(subf).y{:,:},1,"omitnan");
        RG_sub.(roi).Results.(subf).y_std = std(RG_sub.(roi).Results.(subf).y{:,:},[],1,"omitnan");
        n_actual = sum(~isnan(RG_sub.(roi).Results.(subf).y{:,:}),1);
        RG_sub.(roi).Results.(subf).y_sem = RG_sub.(roi).Results.(subf).y_std ./ sqrt(n_actual);

    end

    % --- Subset subject_ids ---
    if isfield(RG.(roi), 'subject_ids')
        if ~isequal(numel(RG.(roi).subject_ids), numel(subIdx))
            error('ROI "%s": subIdx length (%d) does not match subject_ids length (%d).', ...
                roi, numel(subIdx), numel(RG.(roi).subject_ids));
        end
        RG_sub.(roi).subject_ids(~subIdx) = [];
    end

    % --- Subset individual_data ---
    if isfield(RG.(roi), 'individual_data')
        if ~isequal(numel(RG.(roi).individual_data), numel(subIdx))
            error('ROI "%s": subIdx length (%d) does not match individual_data length (%d).', ...
                roi, numel(subIdx), numel(RG.(roi).individual_data));
        end
        RG_sub.(roi).individual_data(~subIdx) = [];
    end

    % --- Warn about user_input_fields ---
    if isfield(RG.(roi), 'user_input_fields')
        if warn
            warning('''user_input_fields'' is not modified by mrGrad_subset.');
            warn = false;
        end
    end
end
end
