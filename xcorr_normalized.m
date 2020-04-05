function [ overlap_12 ] = xcorr_normalized( field1, field2 )
% very simple function that just computes normalized cross correlation of two vectors, 
% field1 and field2

% transpose fields so they are nx1
field1 = reshape( field1, [ length(field1), 1 ] );
field2 = reshape( field2, [ length(field2), 1 ] );

% first normalize both fields
norm_factor_1   = field1' * field1;
norm_factor_2   = field2' * field2;

% put them on same grid size
maxlen = max(length(field1), length(field2));
field1(end+1:maxlen) = 0;
field2(end+1:maxlen) = 0;

% calculate field1 and field2 xcorr
overlap_12 = ifftshift( ifft( fft( fftshift( field1 ) ) .* conj( fft( fftshift( field2 ) ) ) ) );
overlap_12 = overlap_12./sqrt( ( norm_factor_1 .* norm_factor_2 ) );
            
end     % end function xcorr_normalized()