clc
clear all
close all

global kn_min kn_max

% Storage info
% dirname = 'MMNP_ECC_datafiles';
% filename = 'MMNP_ECC_k.txt';
dirname = 'MMNP_ECC_datafiles/MMNP_ECC_k';
filename_mu_vec = 'mu_vec.txt';
filename_kn_vec = 'kn_vec_mu0p%i.txt';     % OBS! One for each mu
filename_error_vec = 'error_vec_mu0p%i.txt';  % OBS! One for each mu

%mu_vec = [0.0, 0.1, 0.2, 0.4, 0.6];
mu_vec = [0.0];

% lekn_vec = 15;
lekn_vec = 3;
kn_vec_init = zeros(lekn_vec,1);
kn_vec = zeros(lekn_vec,1);
error_vec = zeros(lekn_vec,1); 

% kn_min = 3e-2;
% kn_max = 1e0;
kn_min = 5e-3;
kn_max = 5e-1;
kn_vec_init(1) = kn_max;

stepfac = (log10(kn_max) - log10(kn_min))/(lekn_vec - 1);

for k = 1:lekn_vec - 1
    kn_vec_init(k+1) = kn_vec_init(k)*10^-stepfac;
end

% Write mu_vec to file
writematrix(mu_vec, [dirname, '/', filename_mu_vec])

for mu_iter = mu_vec

    % Compute error and kn for given mu
    iter = 0;
    for kn_init = kn_vec_init'
        iter = iter + 1;
        disp(['Iteration = ', num2str(iter)])

        [error, kn, h0] = MMNP_main_ECC(kn_init, mu_iter);
        disp(' ')
        
        kn_vec(iter) = kn;
        error_vec(iter) = error;

    end
    
    % Save kn_vec and error_vec for given mu. 1st kn, 2nd error
%     writematrix(kn_vec', [dirname, '/', filename], 'WriteMode', 'append')
%     writematrix(error_vec', [dirname, '/', filename], 'WriteMode', 'append')
    
    writematrix(kn_vec', [dirname, '/', sprintf(filename_kn_vec, 10*mu_iter)])
    writematrix(error_vec', [dirname, '/', sprintf(filename_error_vec, 10*mu_iter)])

end

disp('Iteration over different mu-values complete.')

% kn_veclog10 = log10(kn_vec);
% error_veclog10 = log10(error_vec);
% 
% % Plot error convergence w.r.t. kn 
% cfun_error = polyfit(kn_veclog10, error_veclog10,1)
% 
% figure(1)
% 
% loglog(kn_vec, error_vec, 'bo', 'Linewidth', 2)
% hold on
% loglog(kn_vec, 10.^(polyval(cfun_error, kn_veclog10)),'r-', 'Linewidth', 2)
% loglog(kn_vec, kn_vec.^0.5*(error_vec(1)/kn_vec(1)^0.5), 'k-.')
% loglog(kn_vec, kn_vec.^1.0*(error_vec(1)/kn_vec(1)^1.0), 'k-')
% loglog(kn_vec, kn_vec.^1.5*(error_vec(1)/kn_vec(1)^1.5), 'k--')
% loglog(kn_vec, error_vec, 'bo', 'Linewidth', 2)
% 
% error_str = '$\left|\mkern-1.5mu\left|\mkern-1.5mu\left| e \right|\mkern-1.5mu\right|\mkern-1.5mu\right|_{X}$';
% legend(error_str,'LLS','$k^{0.5}$','$k^{1.0}$','$k^{1.5}$', ...
%     'interpreter','latex','Fontsize', 24, 'Location', 'SouthEast');
% xlabel('$k$','interpreter','latex','Fontsize', 24);
% ylabel(error_str,'interpreter','latex','Fontsize', 24);
% grid