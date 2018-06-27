function [Ez, Hx] = f_fiberModeGaussian( w0, lambda0, xvec, theta, d0, nclad )
% somewhat adapted from Cale's code
%
% Generate Gaussian-beam mode profile at a plane through y = d0 at angle of
% theta (in degrees)
%
% NOTE that the angle theta is with respect to the positive y axis. When
% overlapping downwards, a "positive theta" in our typical definition (from
% the downwards normal) is actually a negative theta with respect to the
% positive y axis.
%
% currently solving for TE mode only
%
% (using H.A. Haus, Waves & Fields, Chapter 5)
% CMG November 21st, 2014
%
%
% Inputs:   
%   w0  
%       type: double, scalar
%       desc: 1/e beam RADIUS at waist, in units 'units'
%           
%   lambda0
%       type: double, scalar
%       desc: free space wavelength, in units 'units'
%
%   xvec
%       type: double, array
%       desc: coordinates along direction of propagation of
%             grating, in units 'units'
%
%   unit_scale
%       type: double, scalar
%       desc: scaling factor, such that 'units' * unit_scale = meters
%               p sure this isn't needed
%           
%   theta 
%       type: double, scalar
%       desc: angle from normal in degrees
%
%   d0  
%       type: double, scalar
%       desc: distance from beam waist to slice
%
%   nclad 
%       type: double, scalar
%       desc: cladding index
%
%
% Outputs: 
%   u
%       type: double, array
%       desc: returned slice of gaussian beam, normalized to total
%             power


% Constants, in units of meters
% lambda0 = lambda0 * unit_scale;                 % m
% c       = 3e8;                                  % m/s
% mu0     = 4*pi * 1e-7;                          % H/m
% omega0  = 2*pi*c/lambda0;                       % rad/s
% lambda  = lambda0 / nclad;                      % wavelength in cladding, units m
% k       = 2*pi/lambda;                          % 1/m
% w0      = w0 * unit_scale;                      % [meters] radius
% d0      = d0 * unit_scale;                      % [meters] offset
c       = 3e8;                                  % m/s
mu0     = 4*pi * 1e-7;                          % H/m
omega0  = 2*pi*c/lambda0;                       % rad/s
lambda  = lambda0 / nclad;                      % wavelength in cladding
k       = 2*pi/lambda;                          

% Convert to radians
theta = (pi/180)*theta;

% Scale coordinates
% xvec = xvec * unit_scale;                                              % units m

% coordinates in fiber frame
xprime = xvec.*cos(theta) - d0*sin(theta);
zprime = xvec.*sin(theta) + d0*cos(theta);

% b (confocal parameters) is used instead of z0 so that z0 = -1j.*b removes the singularity of the solution on the real z axis (see Haus pg 109)
b = k*w0^2/2;                                                                                   

% Equation (5.2) in Haus [1/meters]
u00 =   1j .* sqrt(k*b/pi) .* ( 1./(zprime + 1j.*b) ).*...
    exp( -1j.*k.*( xprime.^2 )./( 2*(zprime + 1j.*b) ) );     

% % normalize the mode to intensity, makes things nicer
% u00   = u00/sqrt( sum( abs( u00 ).^2 ) );


% calculate fields
E_yprime = -1i * omega0 * u00 .* exp( -1i * k * zprime );
H_xprime = ( -1i * k/mu0 ) * u00 .* exp( -1i * k * zprime );
H_zprime = ( -xprime./( zprime + 1i*b ) ) .* H_xprime; 

% project fields back into original coordinates
Ez = E_yprime;                                                              % this one is used (TE pol grating)
Hx = H_xprime * cos(theta) + H_zprime * sin(theta) ;



end     % end fiberModeGaussian()












































