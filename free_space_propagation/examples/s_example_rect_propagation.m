% Fourier Optics, Homework #2 - Bohan Zhang
% Part 1
clear; clc; close all;

% add path to upper folder
addpath('..');

% % save img path (end it with a '\')
% % make sure the folder already exists beforehand
% figs_path   = 'C:\Users\beezy_000\Google Drive\CU Boulder\coursework\2016 Fall\Fourier Optics\HW\hw2\latex\figs\part2_e\';
% save_off    = true;

% -------------------------------------------------------------------------
% Parameters
% -------------------------------------------------------------------------

% assuming an incident E field of a plane wave
lambda0 = 1e-6;                 % 1 um wavelength
n0      = 1;                    % medium index of refraction
lambda  = lambda0/n0;           % material wavelength
k       = 2*pi/lambda;          % material wavenumber

% Defining the space for proper DFT
L   = 2*lambda;               % length of slit
W   = 100*L;
dx  = (lambda/4);
Nx  = W/dx;
x   = ( -Nx/2:((Nx/2)-1) ) .* dx;               % space vector (centered around 0) (m)
fx  = ( -Nx/2:((Nx/2)-1) ) .* (1/(Nx*dx));      % frequency vector (centered around 0) (1/m)

% propagation distance(s)
n_z             = 100;                          % number of Z's i wanna make
ff_factor       = 1;                            % how many far fields i wanna go to
dz              = 0.1e-6;
z_final         = 10e-6;
z               = 0:dz:z_final; % z's to propagate to

% -------------------------------------------------------------------------
% Propagate E field
% -------------------------------------------------------------------------

% E field is a rect
e_x0 = double( abs(x-L/4) <= L/2);

% Propagate field
[ E_XZ_l, FZ_l, E_fxfz_l ] = f_free_space_propagation( x, z, e_x0, lambda0, n0 );

% Plot the field in 2D space, amp
figure;
imagesc( z, x, abs( E_XZ_l ) );
title('Propagated field, absolute value');
xlabel('z (m)'); ylabel('x (m)');
colorbar
ylim([ -10*L, 10*L ]);
set(gca,'YDir','normal')

% Plot the field in 2D space, real
figure;
imagesc( z, x, real( E_XZ_l ) );
title('Propagated field, real component');
xlabel('z (m)'); ylabel('x (m)');
colorbar
ylim([ -10*L, 10*L ]);
set(gca,'YDir','normal')







