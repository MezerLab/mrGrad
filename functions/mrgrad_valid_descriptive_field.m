function is_valid = mrgrad_valid_descriptive_field(var)
    try
        string(var);
        is_valid = true;
    catch
        % If conversion fails, return false
        is_valid = false;
    end
end
