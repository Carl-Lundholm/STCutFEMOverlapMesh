clc
clear all
%close all

global kn_min kn_max

lekn_vec = 15;
%lekn_vec = 5;
kn_vec_init = zeros(lekn_vec,1);
kn_vec = zeros(lekn_vec,1);
error_vec = zeros(lekn_vec,1);

%kn_min = 1e-4;
% kn_max = 1e-0;
% kn_min = 5e-5;
kn_max = 5e-1;
kn_min = 5e-3;
% kn_min = 1e-1;
% kn_min = 1e-2;
% kn_max = 5e-1;
kn_vec_init(1) = kn_max;

stepfac = (log10(kn_max) - log10(kn_min))/(lekn_vec - 1);
for k = 1:lekn_vec - 1
    kn_vec_init(k+1) = kn_vec_init(k)*10^-stepfac;
end

lehfix_vec = 3; % If "= 1", then it is as in the non-MANY standard case. 
hfix_min = 1e-3; % hfix in orig fig
% hfix_max = 1e-1; 
hfix_max = 1.5e-1;
hfix_vec = zeros(lehfix_vec,1);
hfix_vec(1) = hfix_min;
stepfac = (log10(hfix_max) - log10(hfix_min))/(lehfix_vec - 1);
for hfix_iter = 1:lehfix_vec-1
    hfix_vec(hfix_iter+1) = hfix_vec(hfix_iter)*10^stepfac;
end

hfix_vec(2) = 7e-2;
%hfix_vec(2) = 1e-2; % Looks a bit weird due to earlier saturation bump

hfix_symbvec = ['o', 'x', 'd'];
for hfix_iter = 2:lehfix_vec
    
    hfix = hfix_vec(hfix_iter);
    iter = 0;
    for kn_init = kn_vec_init'
        hfix
        iter = iter + 1
    
        [error, kn, h0] = main_ECC(kn_init, hfix);
        kn_vec(iter) = kn;
        error_vec(iter) = error;
    end

% Plotting of more hfix cases for hfix MANY. Orig fig should be open! 
figure(1)
hold on
loglog(kn_vec, error_vec, ['b',hfix_symbvec(hfix_iter)], 'Linewidth', 2)
error_str = '$\left|\mkern-1.5mu\left|\mkern-1.5mu\left| e \right|\mkern-1.5mu\right|\mkern-1.5mu\right|_{X}$';
legend(error_str,'LLS','$k^{0.5}$','$k^{1.0}$','$k^{1.5}$', ...
    'interpreter','latex','Fontsize', 24, 'Location', 'SouthEast');

end

hfix_vec

% % Plotting of original figure, to be commented out for hfix MANY
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