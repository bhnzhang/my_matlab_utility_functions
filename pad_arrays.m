function [ array1, array2 ] = pad_arrays( array1, array2 )
% Zero pads two arrays of different length so they are the same length 
% so I can do FFTs properly

% get longer length of the two
longer_len = max( [ length(array1), length(array2) ] );

% pad array 1
pad_len     = longer_len - length(array1);
pad_low     = ceil( pad_len/2 );
pad_high    = floor( pad_len/2 );
array1      = [ zeros(1,pad_low), array1, zeros(1,pad_high) ];

% pad array 2
pad_len     = longer_len - length(array2);
pad_low     = ceil( pad_len/2 );
pad_high    = floor( pad_len/2 );
array2      = [ zeros(1,pad_low), array2, zeros(1,pad_high) ];

end

