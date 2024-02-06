clc
clear all
close all

lekn_vec = 15;
kn_vec = zeros(lekn_vec,1);

kn_min = 5*1e-3;
kn_max = 1e-0;
kn_vec(1) = kn_max;

stepfac = (log10(kn_max) - log10(kn_min))/(lekn_vec - 1);

for k = 1:lekn_vec - 1
    kn_vec(k+1) = kn_vec(k)*10^-stepfac;
end

errdiff_vec = zeros(lekn_vec,1);
errq2p_vec = zeros(lekn_vec,1);
errexa_vec = zeros(lekn_vec,1); 
delxh0_vec = zeros(lekn_vec,1); 

iter = 0;
for k_ECC = kn_vec';
    iter = iter + 1
    
    [errexa, errq2p, delxh0] = main_ECC(k_ECC);

    errexa_vec(iter) = errexa;
    errq2p_vec(iter) = errq2p;
    delxh0_vec(iter) = delxh0;
    
    errdiff_vec(iter) = abs(errexa - errq2p);

end

kn_veclog10 = log10(kn_vec);
errexa_veclog10 = log10(errexa_vec);
errq2p_veclog10 = log10(errq2p_vec);

cfun_errexa = fit(kn_veclog10,errexa_veclog10,'poly1')
%cfun_errq2p = fit(kn_veclog10,errq2p_veclog10,'poly1')

figure(1)
loglog(kn_vec, errexa_vec, 'bo')
hold on
loglog(kn_vec, 10.^(cfun_errexa(kn_veclog10)),'r-')
%loglog(kn_vec, errq2p_vec, 'rx')
%loglog(kn_vec, 10.^(cfun_errq2p(kn_veclog10)))
%loglog(kn_vec, kn_vec, 'ko-')
%legend('exa err','q2p err','interpreter','latex'); %,'k^1')
legend('error','lls of error','interpreter','latex'); %,'k^1')
xlabel('$k_n$','interpreter','latex');
ylabel('$error$','interpreter','latex');
grid 

% figure(2)
% plot(kn_vec, delxh0_vec ,'ko')
% %loglog(kn_vec, delxh0_vec ,'ko')
% hold on
% %loglog([1e-3, 1], [1 1], 'k-')
% xlabel('$k_n$','interpreter','latex');
% ylabel('$\Delta x/ h_0$','interpreter','latex');
% grid 

maxdiff = max(errdiff_vec)