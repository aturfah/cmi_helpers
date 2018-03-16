function invalid = invalid_C(C)
%INVALID_C Check if C has been set
%   Since C is a cell array, this can be done
%   with isempty().
invalid = isempty(C);
if invalid
    debug_print("You must set C.");
end

end

