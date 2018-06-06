% Fourier Optics, Homework #3 - Bohan Zhang
% Part c
clear; clc; close all;

% save img path (end it with a '\')
% make sure the folder already exists beforehand
figs_path = 'C:\Users\beezy_000\Google Drive\CU Boulder\coursework\2016 Fall\Fourier Optics\HW\hw2\latex\figs\part2_e\';
save_off = true;
DEBUG = true;

%% Parameters

% assuming an incident E field of a plane wave
lambda0 = 1e-6;         % 1 um wavelength
n0 = 1;                 % medium index of refraction
lambda = lambda0/n0;    % material wavelength
k = 2*pi/lambda;        % material wavenumber

% Talbot grating parameters
m = 8;              % taken from book
% qmax = 6;
L = 10*lambda;      % period of the grating
% fo_grat = 10/L;   

% Defining the space for proper DFT
% z_ff = ((L/2)^2)/lambda;   % propagate into fraunhoffer
W = 10*L;
dx = (lambda/8);
Nx = W/dx;
x = ( -Nx/2:((Nx/2)-1) ) .* dx;               % space vector (centered around 0) (m)
fx = ( -Nx/2:((Nx/2)-1) ) .* (1/(Nx*dx));     % frequency vector (centered around 0) (1/m)

% Lens parameters
% f = z_ff;                       % focal point = object distance
% D = [50] .* 2*(lambda*z_ff)/L;

% propagation distance(s)
zmax = 2*(L^2)/lambda;
n_z = 100;                           % number of Z's i wanna make
ff_factor = 1;                      % how many far fields i wanna go to
z = linspace(0,ff_factor,n_z).*zmax; % z's to propagate to
% z = z_ff*ff_factor;


%% Generate incident field

e_x0 = (1/2)*(1 + m * cos(2 * pi .*x/L));
% e_x0 = exp( 1i * (m/2) .* sin( 2*pi*fo_grat.*x) );

% plot the incident field
% figure;
% plot(x, angle(e_x0));
% xlabel('x (m)'); ylabel('Incident E field Phase (radians)'); title('Phase E(x,0)');
% save_fig_multiformat( gcf, figs_path, ['Ex_0'], save_off);
% 
% % plot the incident field
% figure;
% plot(x, angle(e_x0));
% xlabel('x (m)'); ylabel('Incident E field Phase (radians)'); title('Phase E(x,0)');
% xlim([ -L/2, L/2 ]);
% save_fig_multiformat( gcf, figs_path, ['Ex_0'], save_off);

% % Compute the R.S. Propagation to the lens
% [ E_XZ_l, FZ_l, E_fxfz_l ] = RS_prop( x, fx, e_x0, z, lambda0, n0 );
% 
% % take the e field at the lens
% e_xz_l = E_XZ_l(end,:);
% fz_l = FZ_l(end,:);
% 
% % Multiply by lens factor
% t_l = exp(-1i* (k/(2*f)) .* (x.^2) );
% pupil_l = double(abs(x) <= D/2);
% e_xz_l_plus = e_xz_l .* t_l .* pupil_l;

% % Compute the R.S. Propagation to the image
% [ E_XZ, FZ, E_fxfz ] = RS_prop( x, fx, e_xz_l_plus, z, lambda0, n0 );

% Compute the R.S. Propagation to the image
[ E_XZ, FZ, E_fxfz ] = RS_prop( x, fx, e_x0, z, lambda0, n0, DEBUG );

% take the e field at zmax
e_xz = E_XZ(end,:);
fz = FZ(end,:);

fprintf('Finished RS Prop\n\n');

%% Plotting field


% plot the talbot field
figure;
plot(x, (e_xz)); hold on;
plot(x, E_XZ(end/2,:));
plot(x, E_XZ(end/4,:));
xlabel('x (m)'); ylabel('E Field');
title('E(x,z) of Talbot images');
legend('Talbot', 'Phase Reversed Talbot', 'Talbot subimage');
save_fig_multiformat( gcf, figs_path, ['Ex_zmax'], save_off);
makeFigureNice();

% plot the inverted talbot field
figure;
plot(x, (e_xz));
xlabel('x (m)'); ylabel('E Field');
title('E(x,z) at Phase reversed Talbot image');
save_fig_multiformat( gcf, figs_path, ['Ex_zmax'], save_off);

% plot the talbot field
figure;
plot(x, (e_xz));
xlabel('x (m)'); ylabel('E Field');
title('E(x,z) at Talbot image');
save_fig_multiformat( gcf, figs_path, ['Ex_zmax'], save_off);

% Plot the field in 2D space, amp
figure;
imagesc( z, x, (abs([E_XZ])).' );
title('Amplitude of propagated field');
xlabel('z (m)'); ylabel('x (m)');
colorbar
% ylim([ -10*L, 10*L ]);
set(gca,'YDir','normal')
save_fig_multiformat( gcf, figs_path, ['Ex_z_all'], save_off);

% Plot the field in 2D space, real
figure;
imagesc( z, x, real([E_XZ].') );
title('Propagated Field');
xlabel('z (m)'); ylabel('x (m)');
colorbar
% ylim([ -10*L, 10*L ]);
set(gca,'YDir','normal')
save_fig_multiformat( gcf, figs_path, ['Ex_z_all'], save_off);

% % Plot the field in 2D space, amp
% figure;
% imagesc( [z, z+z(end)], x, (abs([E_XZ_l; E_XZ])).' );
% title('Amplitude of propagated field');
% xlabel('z (m)'); ylabel('x (m)');
% colorbar
% % ylim([ -10*L, 10*L ]);
% set(gca,'YDir','normal')
% save_fig_multiformat( gcf, figs_path, ['Ex_z_all'], save_off);

% %% Analytical solution
% zmax = max(z);
% % x = x.*(1e4);
% U = (1/(1i*lambda*zmax)) .* exp(1i*k*zmax) .* exp(1i * (k/(2*zmax)) .* (x.^2) );
% q = -qmax:1:qmax;
% s_sincs = zeros(size(U));
% for ii = 1:length(q)
%     s_sincs = s_sincs + besselj(q(ii), m/2).*sinc( (L/(lambda*zmax)) .* (x - q(ii) * fo_grat * zmax * lambda) );
% end
% 
% U = U .* s_sincs;
% 
% % plot analytical
% figure;
% plot( x, abs(U) );
% xlabel('x (m)'); ylabel('U (Abs)');
% title('Theoretical field');
% % makeFigureNice();
% 
% % plot analytical and numerical
% figure;
% plot( x, abs(e_xz)./max(abs(e_xz)), x, abs(U)./max(abs(U)) );
% xlabel('x (m)'); ylabel('E field (abs)');
% title('Numerical versus theoretical field');
% legend('Numerical', 'Analytical');
% makeFigureNice();







