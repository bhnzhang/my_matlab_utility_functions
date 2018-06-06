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

% Lens parameters, uhhh what are these for?
L           = 100*lambda;               % length of slit
z_l         = 1*((L/2)^2)/lambda;       % distance to lens (m)
f           = z_l/2;
theta_parax = atan( lambda/L ); % paraxial angle
D           = [20] .* 2*(lambda*z_l)/L;
NA          = n0*D/(2*f);

% Defining the space for proper DFT
W   = 100*L;
dx  = (lambda/4);
Nx  = W/dx;
x   = ( -Nx/2:((Nx/2)-1) ) .* dx;               % space vector (centered around 0) (m)
fx  = ( -Nx/2:((Nx/2)-1) ) .* (1/(Nx*dx));      % frequency vector (centered around 0) (1/m)

% propagation distance(s)
n_z             = 100;                          % number of Z's i wanna make
ff_factor       = 1;                            % how many far fields i wanna go to
z               = linspace(0,1,n_z).*z_l;       % z's to propagate to

% -------------------------------------------------------------------------
% Generate incident field
% -------------------------------------------------------------------------

e_x0 = double( abs(x-L/4) <= L/2);

% % plot the incident field
% figure;
% plot(x, e_x0);
% xlabel('x (m)'); ylabel('Incident E field'); title('E(x,0)');
% save_fig_multiformat( gcf, figs_path, ['Ex_0'], save_off);

% Compute the R.S. Propagation to the lens
% [ E_XZ_l, FZ_l, E_fxfz_l ] = RS_prop( x, fx, e_x0, z, lambda0, n0 );
[ E_XZ_l, FZ_l, E_fxfz_l ] = f_free_space_propagation( x, e_x0, z, lambda0, n0 );

% take the e field at the lens
e_xz_l  = E_XZ_l(end,:);
fz_l    = FZ_l(end,:);

e_xz_image_all = [];
for ii = 1:length(D)
    
    fprintf('running RS prop for D = %f\n', D(ii));

    % Multiply by lens factor
    t_l = exp(-1i* (k/(2*f)) .* (x.^2) );
    pupil_l = double(abs(x) <= D(ii)/2);
    e_xz_l_plus = e_xz_l .* t_l .* pupil_l;

    % Compute the R.S. Propagation to the lens
    [ E_XZ_image, FZ_image, E_fxfz_image ] = RS_prop( x, fx, e_xz_l_plus, z, lambda0, n0 );
    
    % take the e field at the maximum z to plot
    e_xz_image_all = [e_xz_image_all; E_XZ_image(end,:)];
    
end

fprintf('Finished RS Prop\n\n');

%% Plotting fields at lens

% % plot the propagated field
% figure;
% plot(x, abs(e_xz_l));
% xlabel('x (m)'); ylabel('Propagated E Field');
% title('E(x,z) (abs)');
% save_fig_multiformat( gcf, figs_path, ['Ex_zmax'], save_off);
% 
% % Plot the field in 2D space, amp
% figure;
% imagesc( z, x, (abs(E_XZ_l))');
% title('Amplitude of propagated field up to lens');
% xlabel('z (m)'); ylabel('x (m)');
% colorbar
% set(gca,'YDir','normal')
% save_fig_multiformat( gcf, figs_path, ['Ex_z_all'], save_off);
% 
% % Plot the field in 2D space, amp
% figure;
% imagesc( z, x, log10((abs(E_XZ_l))'));
% title('Amplitude of propagated field up to lensm, log scale');
% xlabel('z (m)'); ylabel('x (m)');
% colorbar
% set(gca,'YDir','normal')
% save_fig_multiformat( gcf, figs_path, ['Ex_z_all'], save_off);

%% Plotting image

% take the e field at the maximum z to plot
e_xz_image = E_XZ_image(end,:);
fz_image = FZ_image(end,:);

% plot the propagated field
figure;
plot(x, abs(e_xz_image));
xlabel('x (m)'); ylabel('Propagated E Field');
title('E(x,z) (abs) at image');
save_fig_multiformat( gcf, figs_path, ['Ex_zmax'], save_off);

% Plot the field in 2D space, amp
z_all = [z, z+z_l];
figure;
imagesc( z_all, x, (abs([E_XZ_l; E_XZ_image])).' );
title('Amplitude of propagated field from object to image');
xlabel('z (m)'); ylabel('x (m)');
colorbar
ylim([ -10*L, 10*L ]);
set(gca,'YDir','normal')
save_fig_multiformat( gcf, figs_path, ['Ex_z_all'], save_off);

% plot the propagated field versus incident field
figure;
plot(x, abs(e_xz_image), 'b-', 'LineWidth', 1.5); hold on;
plot(x, abs(e_x0), 'r-', 'LineWidth', 1.5);
xlabel('x (m)'); ylabel('E Field');
legend('Field at image', 'Incident Field');
title('E(x,z) (abs) at image versus incident field');
xlim([min(x)/10.0, max(x)/10.0]);
save_fig_multiformat( gcf, figs_path, ['Ex_zmax'], save_off);

% plot image versus NA
figure;
plot(x, abs(e_xz_image_all), x, e_x0, '--');
xlabel('x (m)'); ylabel('E field');
title('Comparison of field at image versus NA');
xlim( [ -L, L  ] );
legendstr = {};
for ii = 1:length(D)
   legendstr{end+1} = ['NA = ', num2str(NA(ii))];
end
legendstr{end+1} = 'Object';
legend(legendstr);
makeFigureNice( [800, 500] );

% 
% % plot the propagated field
% % figure;
% % plot(x, (angle(e_xz)));
% % xlabel('x (m)'); ylabel(['Propagated E Field @ Z = ' + num2str(z(end))]);
% % title('E(x,z) (phase)');
% 
% % Plot the field in 2D space, amp
% figure;
% imagesc( z, x, (abs(E_XZ))');
% title('Amplitude of propagated field');
% xlabel('z (m)'); ylabel('x (m)');
% colorbar
% set(gca,'YDir','normal')
% save_fig_multiformat( gcf, figs_path, ['Ex_z_all'], save_off);
% 
% % Plot the field in 2D space, phase
% figure;
% imagesc( z, x, (angle(E_XZ))');
% title('Phase of propagated field');
% xlabel('z (m)'); ylabel('x (m)');
% colorbar
% set(gca,'YDir','normal')
% 
% % Plot the field in 2D space, real
% figure;
% imagesc( z, x, (abs(real(E_XZ)))');
% title('Real propagated field');
% xlabel('z (m)'); ylabel('x (m)');
% colorbar
% set(gca,'YDir','normal')
% 
% % Plot the field in 2D space, imag
% figure;
% imagesc( z, x, (imag(E_XZ))');
% title('Imag propagated field');
% xlabel('z (m)'); ylabel('x (m)');
% colorbar
% set(gca,'YDir','normal')
% 
% % Plot power
% S = abs(E_XZ).^2; % intensity
% P = sum( S, 2); % integral of intensity
% figure;
% plot(z, P);
% xlabel('z (m)'); ylabel('Total power'); title('Power vs. Z');

%% The analytical Solution

% % plot the fraunhoffer limit (far field)
% % max_z = z(end);
% z = z( z >= ff );
% [X,Z] = meshgrid(x,z);
% e_xz_fraun = (exp(1i*Z*k)./(1i*lambda*Z)) .* (exp(1i*k*(X.^2)./(2.*Z))) .* sinc(L.*X./(lambda*Z));
% E_analytical = e_xz_fraun;
% 
% % Plot the field in 2D space
% figure;
% imagesc( z, x, (abs(E_analytical))');
% title('Amplitude of propagated field, analytical');
% xlabel('z (m)'); ylabel('x (m)');
% colorbar
% set(gca,'YDir','normal')
% save_fig_multiformat( gcf, figs_path, ['Ex_z_all_fraun'], save_off);
% 
% % Plot the field in 2D space
% % figure;
% % imagesc( z, x, (angle(E_analytical))');
% % title('Phase of propagated field, analytical');
% % xlabel('z (m)'); ylabel('x (m)');
% % colorbar
% % set(gca,'YDir','normal')
% 
% % Plot the field in 2D space
% % figure;
% % imagesc( z, x, (real(E_analytical))');
% % title(['E(x,z) versus z, real, Analytical solution, \theta_i = ', num2str(theta_i*180/pi) ]);
% % xlabel('z (m)'); ylabel('x (m)');
% % colorbar
% % set(gca,'YDir','normal')
% % save_fig_multiformat( gcf, figs_path, ['Exz_vs_z_analytical']);
% % 
% % % Plot the field in 2D space
% % figure;
% % imagesc( z, x, (imag(E_analytical))');
% % title('Imag propagated field, analytical');
% % xlabel('z (m)'); ylabel('x (m)');
% % colorbar
% % set(gca,'YDir','normal')
% 
% % e_xz_fraun = e_xz_fraun./length(x);
% % e_xz_fraun = L .* sinc(L.*x./(lambda*z));
% 
% % Plot fraunhoffer vs. angle at max z
% theta = 180*atan(x/z(end))/pi;
% figure;
% plot(theta, abs(E_analytical(end,:) ) );
% xlabel('Angle (degrees)'); ylabel('E');
% title(['E(x,z) Fraunhoffer (amplitude) vs. angle at z = ', num2str(z(end)), ' m']);
% save_fig_multiformat( gcf, figs_path, ['Ex_zmax_fraun'], save_off );
% 
% % Plot comparison between fraunhoffer and actual
% figure;
% plot(theta, abs(e_xz)./max(abs(e_xz)), 'LineWidth', 2 ); hold on;
% plot(theta, abs(E_analytical(end,:))./max(abs(E_analytical(end,:))), 'LineWidth', 2, 'LineStyle', '--' ); 
% xlabel('Angle (degrees)'); ylabel('E');
% title('Comparison of Rayleigh Sommerfield and Fraunhoffer');
% legend( 'Numerical', 'Fraunhoffer');
% save_fig_multiformat( gcf, figs_path, ['exz_num_vs_fraun_cmp'], save_off );
% 
% % plot transfer function fz vs fx
% figure;
% plot( fx, real(FZ(1,:)), 'b-', 'LineWidth', 2 ); hold on;
% plot( fx, imag(FZ(1,:)), 'r-', 'LineWidth', 2);
% xlabel('f_x (m^{-1})'); ylabel('f_z (m^{-1})');
% legend('Real', 'Imaginary');
% title('f_z versus f_x');
% 
% % plot transfer function fz vs n
% figure;
% plot( real(FZ(1,:)), 'b-', 'LineWidth', 2 ); hold on;
% plot( imag(FZ(1,:)), 'r-', 'LineWidth', 2);
% xlabel('index'); ylabel('f_z (m^{-1})');
% legend('Real', 'Imaginary');
% title('f_z versus f_x');

% %% Plotting power vs spatial freq
% 
% % indexes of the f_x's i wanna choose
% n_fx1 = 2300;
% n_fx2 = 3000;
% n_fx3 = 3300;
% n_fx4 = 4000;
% n_fx =[2300 3000 3300 4000].*(400/4000);
% 
% % spectrums at those f_x's
% Efx1 = ifftshift(E_fxfz(:,n_fx(1)));
% Efx2 = ifftshift(E_fxfz(:,n_fx(2)));
% Efx3 = ifftshift(E_fxfz(:,n_fx(3)));
% Efx4 = ifftshift(E_fxfz(:,n_fx(4)));
% 
% % powers
% Pfx1 = abs(Efx1).^2;
% Pfx2 = abs(Efx2).^2;
% Pfx3 = abs(Efx3).^2;
% Pfx4 = abs(Efx4).^2;
% 
% % plotting
% figure;
% plot( z, Pfx1, 'b-' ); hold on
% plot( z, Pfx2, 'r-' );
% plot( z, Pfx3, 'b--' );
% plot( z, Pfx4, 'r--' );
% legendstrs = {};
% for ii = 1:4
%     legendstrs{end+1} = ['f_x = ', num2str(fx(n_fx(ii))), ' m^{-1}'];
% end
% xlabel('z (m)'); ylabel('Power');
% legend(legendstrs);
% title('Power of E(f_x) for different f_x vs. z');
% 
% % theory
% fz_chosen = FZ(:,n_fx);
% Pdecay = exp(1i * 2 * 2 * pi * fz_chosen .* repmat(z',1,4));
% 
% % plotting
% figure;
% plot( z, abs(Pdecay(:,1)), 'r-', 'LineWidth', 2 ); hold on;
% plot( z, Pfx1./max(Pfx1), 'b--', 'LineWidth', 2  );
% legend('Analytical', 'Numerical');
% xlabel('z (m)'); ylabel('Power');
% title(['Power vs. z for f_x = ', num2str( fx(n_fx(1))), ' m^{-1}' ]);
% save_fig_multiformat( gcf, figs_path, ['p_vs_z_fx1'], save_off);
% 
% figure;
% plot( z, abs(Pdecay(:,3)), 'r-', 'LineWidth', 2 ); hold on;
% plot( z, Pfx3./max(Pfx3), 'b--', 'LineWidth', 2  );
% legend('Analytical', 'Numerical');
% xlabel('z (m)'); ylabel('Power');
% title(['Power vs. z for f_x = ', num2str( fx(n_fx(3))), ' m^{-1}' ]);
% save_fig_multiformat( gcf, figs_path, ['p_vs_z_fx3'], save_off);
% 
% figure;
% plot( z, abs(Pdecay(:,2)), 'r-', 'LineWidth', 2 ); hold on;
% plot( z, Pfx2./max(Pfx2), 'b--', 'LineWidth', 2  );
% legend('Analytical', 'Numerical');
% xlabel('z (m)'); ylabel('Power');
% title(['Power vs. z for f_x = ', num2str( fx(n_fx(2))), ' m^{-1}' ]);
% save_fig_multiformat( gcf, figs_path, ['p_vs_z_fx2'], save_off);
% 
% figure;
% plot( z, abs(Pdecay(:,4)), 'r-', 'LineWidth', 2  ); hold on;
% plot( z, Pfx4./max(Pfx4), 'b--', 'LineWidth', 2  );
% legend('Analytical', 'Numerical');
% xlabel('z (m)'); ylabel('Power');
% title(['Power vs. z for f_x = ', num2str( fx(n_fx(4))), ' m^{-1}' ]);
% save_fig_multiformat( gcf, figs_path, ['p_vs_z_fx4'], save_off);

% figure;
% plot(x, angle(e_xz_fraun));
% xlabel('x (m)'); ylabel('E');
% title('E(x,z) Fraunhoffer (phase)');
% 
% figure;
% plot(x, abs(e_xz_fraun)./max(abs(e_xz_fraun)) ); hold on;
% plot(x, abs(e_xz)./max(abs(e_xz)) );
% xlabel('x (m)'); ylabel('E');
% title('Comparison of Rayleigh Sommerfield and Fraunhoffer');






