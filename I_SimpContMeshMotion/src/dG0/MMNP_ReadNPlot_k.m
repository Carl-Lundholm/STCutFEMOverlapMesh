clc
clear all
close all

% Storage info
% dirname = 'MMNP_ECC_datafiles';
% filename = 'MMNP_ECC_k.txt';
dirname = 'MMNP_ECC_datafiles/MMNP_ECC_k';
filename_mu_vec = 'mu_vec.txt';
filename_kn_vec = 'kn_vec_mu0p%i.txt';    % OBS! One for each mu
filename_error_vec = 'error_vec_mu0p%i.txt';  % OBS! One for each mu

% Read data
% DATA = readmatrix([dirname, '/', filename]);

% mu_vec = DATA(1, :);
%mu_vec = readmatrix([dirname, '/', filename_mu_vec]);

mu_vec = [0.0, 0.1, 0.2, 0.4, 0.6];
% mu_vec = [0.6];

mu_instance = 0;
for mu_iter = mu_vec

    % mu instance and value of mu
    mu_instance = mu_instance + 1;
    disp(['mu = ', num2str(mu_iter)])
    
    % Read kn_vec and error_vec 
%     kn_vec = DATA(2*mu_instance, :);
%     error_vec = DATA(2*mu_instance + 1, :);

    kn_vec = readmatrix([dirname, '/', sprintf(filename_kn_vec, 10*mu_iter)]);
    error_vec = readmatrix([dirname, '/', sprintf(filename_error_vec, 10*mu_iter)]);
    
    kn_veclog10 = log10(kn_vec);
    error_veclog10 = log10(error_vec);

    % Plot error convergence w.r.t. kn 
    cfun_error = polyfit(kn_veclog10(4:12), error_veclog10(4:12), 1)
    disp(['Slope of LLS = ', num2str(cfun_error(1))])
    
    figure(1)

    loglog(kn_vec, error_vec, 'bo', 'Linewidth', 2)
    hold on
    loglog(kn_vec(4:12), 10.^(polyval(cfun_error, kn_veclog10(4:12))),'r-', 'Linewidth', 4)
    loglog(kn_vec, kn_vec.^0.5*(error_vec(1)/kn_vec(1)^0.5), 'k-.')
    loglog(kn_vec, kn_vec.^1.0*(error_vec(1)/kn_vec(1)^1.0), 'k-')
    loglog(kn_vec, kn_vec.^1.5*(error_vec(1)/kn_vec(1)^1.5), 'k--')
    loglog(kn_vec, error_vec, 'bo', 'Linewidth', 2)

    error_str = '$\left|\mkern-1.5mu\left|\mkern-1.5mu\left| e \right|\mkern-1.5mu\right|\mkern-1.5mu\right|_{X}$';
    legend(error_str,'LLS','$k^{0.5}$','$k^{1.0}$','$k^{1.5}$', ...
        'interpreter','latex','Fontsize', 24, 'Location', 'SouthEast');
    xlabel('$k$','interpreter','latex','Fontsize', 24);
    ylabel(error_str,'interpreter','latex','Fontsize', 24);
    grid
    
    %pause    % Press any key to continue 
    input('') % ENTER-key needs to be pressed to continue
    
    close all
    
end