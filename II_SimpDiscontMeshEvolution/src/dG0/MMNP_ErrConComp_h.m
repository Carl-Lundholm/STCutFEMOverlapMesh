clc
clear all
close all

% Storage info
% dirname = 'MMNP_ECC_datafiles';
% filename = 'MMNP_ECC_h.txt';
dirname = 'MMNP_ECC_datafiles/MMNP_ECC_h';
filename_mu_vec = 'mu_vec.txt';
filename_h_vec = 'h_vec_mu0p%i.txt';     % OBS! One for each mu
filename_error_vec = 'error_vec_mu0p%i.txt';  % OBS! One for each mu

mu_vec = [0.1, 0.2, 0.4];
% mu_vec = [0.0, 0.6];

leh_vec = 15;
% leh_vec = 3;
h_vec_init = zeros(leh_vec,1);
h_vec = zeros(leh_vec,1);
error_vec = zeros(leh_vec,1); 

h_min = 1e-3;
h_max = 1e-1;
h_vec_init(1) = h_max;

stepfac = (log10(h_max) - log10(h_min))/(leh_vec - 1);

for k = 1:leh_vec - 1
    h_vec_init(k+1) = h_vec_init(k)*10^-stepfac;
end

% Write mu_vec to file
writematrix(mu_vec, [dirname, '/', filename_mu_vec])

for mu_iter = mu_vec

    % Compute error and h for given mu
    iter = 0;
    for h_init = h_vec_init'
        iter = iter + 1;
        disp(['Iteration = ', num2str(iter)])

        [error, kn, h0] = MMNP_main_ECC(h_init, mu_iter);
        disp(' ')
        
        h_vec(iter) = h0;
        error_vec(iter) = error;

    end
    
    % Save h_vec and error_vec for given mu. 1st h, 2nd error
%     writematrix(h_vec', [dirname, '/', filename], 'WriteMode', 'append')
%     writematrix(error_vec', [dirname, '/', filename], 'WriteMode', 'append')

    writematrix(h_vec', [dirname, '/', sprintf(filename_h_vec, 10*mu_iter)])
    writematrix(error_vec', [dirname, '/', sprintf(filename_error_vec, 10*mu_iter)])

end

disp('Iteration over different mu-values complete.')

% h_veclog10 = log10(h_vec);
% error_veclog10 = log10(error_vec);
% 
% % Plot error convergence w.r.t. h 
% cfun_error = polyfit(h_veclog10, error_veclog10,1)
% 
% figure(1)
% 
% loglog(h_vec, error_vec, 'bo', 'Linewidth', 2)
% hold on
% loglog(h_vec, 10.^(polyval(cfun_error, h_veclog10)),'r-', 'Linewidth', 2)
% loglog(h_vec, h_vec.^1*(error_vec(1)/h_vec(1)^1), 'k-.')
% loglog(h_vec, h_vec.^2*(error_vec(1)/h_vec(1)^2), 'k-')
% loglog(h_vec, h_vec.^3*(error_vec(1)/h_vec(1)^3), 'k--')
% loglog(h_vec, error_vec, 'bo', 'Linewidth', 2)
% 
% error_str = '$\| e(T) \|_{L^2(\Omega_0)}$';
% legend(error_str,'LLS','$h^1$','$h^2$','$h^3$', ...
%     'interpreter','latex','Fontsize', 24, 'Location', 'SouthEast');
% xlabel('$h$','interpreter','latex','Fontsize', 24);
% ylabel(error_str,'interpreter','latex','Fontsize', 24);
% grid