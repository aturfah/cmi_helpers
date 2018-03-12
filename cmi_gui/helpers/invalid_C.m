function invalid = invalid_C(C)
%INVALID_C Summary of this function goes here
%   Detailed explanation goes here
invalid = isempty(C);
if invalid
    debug_print("You must select some C.");
end

end

