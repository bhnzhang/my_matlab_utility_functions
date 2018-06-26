function [ E_XZ, FZ, E_fxfz ] = f_free_space_propagation( x, z, e_x0, lambda0, n, DEBUG )
% Written by Bohan Zhang
% 
% Propagates the incident, monochromatic E field at a boundary through free
% space
%
% NOTE: this is an updated ver. of RS_prop() which I was using for fourier
% optics class
% 
% Inputs:
%   x
%       type: double, array
%       desc: transverse coordinates, units arbitrary as long as
%               consistent?
%   z
%       type: double, array
%       desc: propagation coordinates, units arbitrary as long as
%               consistent?
%   e_x0
%       type: double, array
%       desc: transverse electric field to propagate
%   lambda0
%       type: double, scalar
%       desc: free space wavelength, units arbitrary as long as
%               consistent?
%   n
%       type: double, scalar
%       desc: index of refraction of homogenous material
%   DEBUG
%       type: boolean
%       desc: OPTIONAL, set to true to enable debug mode
%     
% outputs:
%   E_XZ
%       type: double, matrix
%       desc: propagated field, dimensions x vs z
%   FZ
%       type: double, matrix
%       desc: propagation constants, dimensions fx vs fz
%   E_fxfz
%       type: double, matrix
%       desc: k space propagated field, dimensions fx vs fz


% use debug mode?
if nargin < 7
    DEBUG = false;
end

% create freq. domain vector
Nx = length(x);
dx = x(2) - x(1);
fx = ( -Nx/2:((Nx/2)-1) ) .* (1/(Nx*dx));                                   % frequency vector (centered around 0) (1/m)

% coordinates, dimensions z vs x
[X, Z]  = meshgrid(x, z);
[FX, Z] = meshgrid(fx, z);

% Create the transfer function, with origin at center of array
FZ      = ( (n/lambda0).^2 - FX.^2 ).^(1/2);                                % dims ( # z, freq )
H_fz    = exp(1i * 2 * pi .* Z .* FZ);                                      % dims ( # z, freq )

% return origin to index 1
x       = ifftshift(x);
e_x0    = ifftshift(e_x0);

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
E_XZ = ifft(E_fxfz, size(E_fxfz,2), 2);                                     % take ifft's along rows

% shift and transpose it before returning for convenience
E_XZ    = fftshift(E_XZ).';
E_fxfz  = fftshift(E_fxfz,2).';
FZ      = FZ.';

end