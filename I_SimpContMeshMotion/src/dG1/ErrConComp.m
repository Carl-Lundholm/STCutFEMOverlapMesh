clc
clear all
close all

global kn_min kn_max

lekn_vec = 15;
kn_vec = zeros(lekn_vec,1);

kn_min = 3e-2;
kn_max = 5e-1;
kn_vec(1) = kn_max;

stepfac = (log10(kn_max) - log10(kn_min))/(lekn_vec - 1);

for k = 1:lekn_vec - 1
    kn_vec(k+1) = kn_vec(k)*10^-stepfac;
end

h_vec = zeros(lekn_vec,1);

errdiff_vec = zeros(lekn_vec,1);
errq2p_vec = zeros(lekn_vec,1);
errexa_vec = zeros(lekn_vec,1); 
delxh0_vec = zeros(lekn_vec,1);

err1L_vec = zeros(lekn_vec,1);
err1R_vec = zeros(lekn_vec,1);
err2_vec = zeros(lekn_vec,1);

iter = 0;
for k_ECC = kn_vec';
    iter = iter + 1
    
    [errexa, errq2p, delxh0, err1L_L2Norm, err1R_L2Norm, err2_L2Norm, ...
        h0] = main_ECC(k_ECC, kn_max);

    h_vec(iter) = h0;
    
    errexa_vec(iter) = errexa;
    errq2p_vec(iter) = errq2p;
    delxh0_vec(iter) = delxh0;
    
    err1L_vec(iter) = sqrt(err1L_L2Norm);
    err1R_vec(iter) = sqrt(err1R_L2Norm);
    err2_vec(iter) = sqrt(err2_L2Norm);
    
    errdiff_vec(iter) = abs(errexa - errq2p);

end

h_veclog10 = log10(h_vec);

errq2p_veclog10 = log10(errq2p_vec);

err1L_veclog10 = log10(err1L_vec);
err1R_veclog10 = log10(err1R_vec);
err2_veclog10 = log10(err2_vec);

%cfun_errq2p = fit(kn_veclog10,errq2p_veclog10,'poly1')

% cfun_err1L = fit(kn_veclog10,err1L_veclog10,'poly1')
% cfun_err1R = fit(kn_veclog10,err1R_veclog10,'poly1')
% cfun_err2 = fit(kn_veclog10,err2_veclog10,'poly1')

kn_veclog10 = log10(kn_vec);
errexa_veclog10 = log10(errexa_vec);

% % Error convergence w.r.t. kn 
% cfun_errexa = fit(kn_veclog10,errexa_veclog10,'poly1')
% figure(1)
% loglog(kn_vec, errexa_vec, 'bo')
% hold on
% loglog(kn_vec, 10.^(cfun_errexa(kn_veclog10)),'r-')
% %loglog(kn_vec, errq2p_vec, 'rx')
% %loglog(kn_vec, 10.^(cfun_errq2p(kn_veclog10)))
% % loglog(kn_vec, kn_vec.^1, 'k-')
% % loglog(kn_vec, kn_vec.^2, 'k-')
% % loglog(kn_vec, kn_vec.^3, 'k-')
% %legend('exa err','q2p err','interpreter','latex'); %,'k^1')
% legend('error','lls of error','interpreter','latex'); %,'k^1')
% % title('$error$ on $\Omega_1(T) \cup \Omega_2(T)$','interpreter','latex');
% xlabel('$k_n$','interpreter','latex');
% ylabel('$error$','interpreter','latex');
% grid 

% PrintError(kn_vec, errexa_vec)

% Error convergence w.r.t. h
cfun_errexa = fit(h_veclog10,errexa_veclog10,'poly1')
figure(1)
loglog(h_vec, errexa_vec, 'bo')
hold on
loglog(h_vec, 10.^(cfun_errexa(h_veclog10)),'r-')
%loglog(kn_vec, errq2p_vec, 'rx')
%loglog(kn_vec, 10.^(cfun_errq2p(kn_veclog10)))
% loglog(h_vec, h_vec.^1, 'k-')
% loglog(h_vec, h_vec.^2, 'k-')
% loglog(h_vec, h_vec.^3, 'k-')
%legend('exa err','q2p err','interpreter','latex'); %,'k^1')
% legend('error','lls of error','k^1','k^2','k^3','interpreter','latex'); %,'k^1')
% title('$error$ on $\Omega_1(T) \cup \Omega_2(T)$','interpreter','latex');
xlabel('$h$','interpreter','latex');
ylabel('$error$','interpreter','latex');
grid 

% figure(2)
% % plot(kn_vec, delxh0_vec ,'ko')
% loglog(kn_vec, delxh0_vec ,'ko')
% hold on
% %loglog([1e-3, 1], [1 1], 'k-')
% xlabel('$k_n$','interpreter','latex');
% ylabel('$\Delta x/ h_0$','interpreter','latex');
% grid 

% figure(3)
% loglog(kn_vec, err1L_vec, 'bo')
% hold on
% loglog(kn_vec, 10.^(cfun_err1L(kn_veclog10)),'r-')
% loglog(kn_vec, kn_vec.^1, 'k-')
% loglog(kn_vec, kn_vec.^2, 'k-')
% loglog(kn_vec, kn_vec.^3, 'k-')
% %legend('exa err','q2p err','interpreter','latex'); %,'k^1')
% % legend('error','lls of error','k^1','k^2','k^3','interpreter','latex'); %,'k^1')
% title('$error$ on $\Omega_1(T)$ to the left of $\Omega_2(T)$','interpreter','latex');
% xlabel('$k_n$','interpreter','latex');
% ylabel('$error$','interpreter','latex');
% grid 
% 
% figure(4)
% loglog(kn_vec, err1R_vec, 'bo')
% hold on
% loglog(kn_vec, 10.^(cfun_err1R(kn_veclog10)),'r-')
% loglog(kn_vec, kn_vec.^1, 'k-')
% loglog(kn_vec, kn_vec.^2, 'k-')
% loglog(kn_vec, kn_vec.^3, 'k-')
% %legend('exa err','q2p err','interpreter','latex'); %,'k^1')
% % legend('error','lls of error','k^1','k^2','k^3','interpreter','latex'); %,'k^1')
% title('$error$ on $\Omega_1(T)$ to the right of $\Omega_2(T)$','interpreter','latex');
% xlabel('$k_n$','interpreter','latex');
% ylabel('$error$','interpreter','latex');
% grid 
% 
% figure(5)
% loglog(kn_vec, err2_vec, 'bo')
% hold on
% loglog(kn_vec, 10.^(cfun_err2(kn_veclog10)),'r-')
% loglog(kn_vec, kn_vec.^1, 'k-')
% loglog(kn_vec, kn_vec.^2, 'k-')
% loglog(kn_vec, kn_vec.^3, 'k-')
% %legend('exa err','q2p err','interpreter','latex'); %,'k^1')
% % legend('error','lls of error','k^1','k^2','k^3','interpreter','latex'); %,'k^1')
% title('$error$ on $\Omega_2(T)$','interpreter','latex');
% xlabel('$k_n$','interpreter','latex');
% ylabel('$error$','interpreter','latex');
% grid 
%maxdiff = max(errdiff_vec)