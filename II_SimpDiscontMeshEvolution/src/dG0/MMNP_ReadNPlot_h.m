clc
clear all
close all

% Storage info
% dirname = 'MMNP_ECC_datafiles';
% filename = 'MMNP_ECC_h.txt';
dirname = 'MMNP_ECC_datafiles/MMNP_ECC_h';
filename_mu_vec = 'mu_vec.txt';
filename_h_vec = 'h_vec_mu0p%i.txt';    % OBS! One for each mu
filename_error_vec = 'error_vec_mu0p%i.txt';  % OBS! One for each mu

% Read in data
%DATA = readmatrix([dirname, '/', filename]);

%mu_vec = DATA(1, :);
% mu_vec = readmatrix([dirname, '/', filename_mu_vec]);

mu_vec = [0.0, 0.1, 0.2, 0.4, 0.6];
% mu_vec = [0.6];

mu_instance = 0;
for mu_iter = mu_vec

    % mu instance and value of mu
    mu_instance = mu_instance + 1;
    disp(['mu = ', num2str(mu_iter)])
    
    % Read in h_vec and error_vec 
%     h_vec = DATA(2*mu_instance, :);
%     error_vec = DATA(2*mu_instance + 1, :);

    h_vec = readmatrix([dirname, '/', sprintf(filename_h_vec, 10*mu_iter)]);
    error_vec = readmatrix([dirname, '/', sprintf(filename_error_vec, 10*mu_iter)]);

    h_veclog10 = log10(h_vec);
    error_veclog10 = log10(error_vec);

    % Plot error convergence w.r.t. h 
    cfun_error = polyfit(h_veclog10(1:11), error_veclog10(1:11),1)
    disp(['Slope of LLS = ', num2str(cfun_error(1))])
    
    figure(1)

    loglog(h_vec, error_vec, 'bo', 'Linewidth', 2)
    hold on
    loglog(h_vec(1:11), 10.^(polyval(cfun_error, h_veclog10(1:11))),'r-', 'Linewidth', 4)
    loglog(h_vec, h_vec.^1*(error_vec(1)/h_vec(1)^1), 'k-.')
    loglog(h_vec, h_vec.^2*(error_vec(1)/h_vec(1)^2), 'k-')
    loglog(h_vec, h_vec.^3*(error_vec(1)/h_vec(1)^3), 'k--')
    loglog(h_vec, error_vec, 'bo', 'Linewidth', 2)

    error_str = '$\| e(T) \|_{L^2(\Omega_0)}$';
    legend(error_str,'LLS','$h^1$','$h^2$','$h^3$', ...
        'interpreter','latex','Fontsize', 24, 'Location', 'SouthEast');
    xlabel('$h$','interpreter','latex','Fontsize', 24);
    ylabel(error_str,'interpreter','latex','Fontsize', 24);
    grid
    
    %pause    % Press any key to continue 
    input('') % ENTER-key needs to be pressed to continue
    
    close all
    
end