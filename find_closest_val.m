function [ closest_val, index ] = find_closest_val( x, val )
% Finds the closest value to val in array x, and returns that value as well
% as the index

[temp, index] = min(abs(x - val));
closest_val = x(index);


end

