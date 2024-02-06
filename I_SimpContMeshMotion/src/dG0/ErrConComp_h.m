clc
clear all
% close all

leh_vec = 15;
% leh_vec = 5;
h_vec_init = zeros(leh_vec,1);
h_vec = zeros(leh_vec,1);
error_vec = zeros(leh_vec,1); 
h_min = 1e-3;
h_max = 1e-1;
% h_min = 5e-2;
% h_max = 5e-1;
h_vec_init(1) = h_max;

stepfac = (log10(h_max) - log10(h_min))/(leh_vec - 1);
for k = 1:leh_vec - 1
    h_vec_init(k+1) = h_vec_init(k)*10^-stepfac;
end

lekfix_vec = 3; % If "= 1", then it is as in the non-MANY standard case. 
kfix_min = 1e-4; % kfix in orig fig
kfix_max = 1e-2;
kfix_vec = zeros(lekfix_vec,1);
kfix_vec(1) = kfix_min;
stepfac = (log10(kfix_max) - log10(kfix_min))/(lekfix_vec - 1);
for kfix_iter = 1:lekfix_vec-1
    kfix_vec(kfix_iter+1) = kfix_vec(kfix_iter)*10^stepfac;
end

kfix_vec(2) = 1e-3;

kfix_symbvec = ['o', 'x', 'd'];
for kfix_iter = 2:lekfix_vec
    
    kfix = kfix_vec(kfix_iter);
    iter = 0;
    for h_init = h_vec_init'
        kfix
        iter = iter + 1
        
        [error, kn, h0] = main_ECC(h_init, kfix); 
        h_vec(iter) = h0;
        error_vec(iter) = error;
    end

% Plotting of more kfix cases for kfix MANY. Orig fig should be open! 
figure(1)
hold on
loglog(h_vec, error_vec, ['b',kfix_symbvec(kfix_iter)], 'Linewidth', 2)
error_str = '$\left|\mkern-1.5mu\left|\mkern-1.5mu\left| e \right|\mkern-1.5mu\right|\mkern-1.5mu\right|_{X}$';
legend(error_str,'LLS','$h^{0.5}$','$h^{1.0}$','$h^{1.5}$', ...
    'interpreter','latex', 'Fontsize', 24, 'Location', 'SouthEast');

end

kfix_vec

% % Plotting of original figure, to be commented out for kfix MANY
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
% loglog(h_vec, h_vec.^0.5*(error_vec(1)/h_vec(1)^0.5), 'k-.')
% loglog(h_vec, h_vec.^1.0*(error_vec(1)/h_vec(1)^1.0), 'k-')
% loglog(h_vec, h_vec.^1.5*(error_vec(1)/h_vec(1)^1.5), 'k--')
% loglog(h_vec, error_vec, 'bo', 'Linewidth', 2)
% 
% error_str = '$\left|\mkern-1.5mu\left|\mkern-1.5mu\left| e \right|\mkern-1.5mu\right|\mkern-1.5mu\right|_{X}$';
% legend(error_str,'LLS','$h^{0.5}$','$h^{1.0}$','$h^{1.5}$', ...
%     'interpreter','latex', 'Fontsize', 24, 'Location', 'SouthEast');
% xlabel('$h$','interpreter','latex', 'Fontsize', 24);
% ylabel(error_str,'interpreter','latex', 'Fontsize', 24);
% grid