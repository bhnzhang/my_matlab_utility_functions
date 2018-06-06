function [ E_XZ, FZ, E_fxfz ] = RS_prop( x, fx, e_x0, z, lambda0, n, DEBUG )
% Written by Bohan Zhang
% 
% TODO:
%     clean up this code sometime
%
% Propagates the incident, monochromatic E field at a boundary
% to some distance z
% To propagate to multiple distances z, make sure the following inputs are appropriate (meshgridded)
%     matrices: x, fx, z
% 
% inputs:
%     x     - space array/matrix, where x = 0 is in the middle of the array
%                 if matrix, then each row is the x array, repeated for the # of rows
%     fx    - frequency array/matrix, same deal as x
%     e_x0  - electric field on boundary given by x
%     z     - distance to propagate to, can be array (a row array). 
%                If x and fx are arrays, z can be a scaler. If x and fx are matrices, then we are
%                propagating to multiple z's, and z is an array
%     lambda0
%         : scalar
%           free-space wavelength
%     n
%         : scalar
%           index of refraction of homogenous material
%     DEBUG   - true if in debug mode (enable extra plotting)
%     
% 
% outputs:
%       e_xz  - electric field at distance z, dim (z, x)
%         fz  - transfer function exponent
%         E_fxfz - E spectrum
        
if nargin < 7
    DEBUG = false;
end

[X, Z]  = meshgrid(x, z);
[FX, Z] = meshgrid(fx, z);

% Create the transfer function, with origin at center of array
FZ = ( (n/lambda0).^2 - FX.^2 ).^(1/2); % dims ( # z, freq )
H_fz = exp(1i * 2 * pi .* Z .* FZ); % dims ( # z, freq )

% return origin to index 1
x = ifftshift(x);
e_x0 = ifftshift(e_x0);

% take the first FFT
E_fx0 = fft( e_x0 );

% DEBUG
% plot the first FFT, aka k space of E at z = 0
if DEBUG
    figure;
    plot( fx, fftshift(abs(E_fx0)) );
    xlabel('fx'); ylabel('Absolute value E_{fx}');
    title('E(fx) at z = 0');
end

% DEBUG
% plot TF overlayed with FFT
if DEBUG
    fz = FZ(1,:);
    fz_norm = fz./max(abs(fz));
    figure;
    plot( fx, real(fz)./max(abs(real(fz))) ); hold on;
    plot( fx, imag(fz)./max(abs(imag(fz))) );
    plot( fx, fftshift(abs(E_fx0)./max(abs(E_fx0))) );
    xlabel('f_x'); ylabel('f_z');
    legend('Real TF', 'Imag TF', 'Abs(E(fx))');
    title('Spectrum and transfer function overplot');
end

% repeat depending on the number of z's to propagate to
E_fx0 = repmat(E_fx0, size(Z,1), 1);

% multiply by transfer function H
E_fxfz = E_fx0.*(ifftshift(H_fz));

% take inverse FT to get the actual field
E_XZ = ifft(E_fxfz, size(E_fxfz,2), 2); % take ifft's along rows

% shift it before returning for convenience
E_XZ = fftshift(E_XZ);
E_fxfz = fftshift(E_fxfz,2);

end