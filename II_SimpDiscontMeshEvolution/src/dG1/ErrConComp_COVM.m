clc
clear all
% close fig

global kn_min kn_max mu_COVM

% Create a logarithmically equidistant mu array 
lemu_vec = 30;
mu_vec = zeros(lemu_vec,1);
mu_min = 1e-7;
mu_max = 1.3;
mu_vec(1) = mu_max;
stepfac_mu = (log10(mu_max) - log10(mu_min))/(lemu_vec - 1);
for m = 1:lemu_vec - 1
    mu_vec(m+1) = mu_vec(m)*10^-stepfac_mu;
end

% Create a logarithmically equidistant kn array 
lekn_vec = 12;            %12
kn_vec = zeros(lekn_vec,1);
kn_min = 5e-2;
kn_max = 5e-1;
kn_vec(1) = kn_max;
stepfac_kn = (log10(kn_max) - log10(kn_min))/(lekn_vec - 1);
for k = 1:lekn_vec - 1
    kn_vec(k+1) = kn_vec(k)*10^-stepfac_kn;
end
kn_veclog10 = log10(kn_vec);

tic
conv_ord_vec = zeros(lemu_vec, 1);
m = 0;
for mu_COVM = mu_vec';
    m = m + 1
    
    errexa_vec = zeros(lekn_vec,1);
    k = 0;
    for k_ECC = kn_vec';
        k = k + 1;
    
        [errexa, errq2p, delxh0, err1L_L2Norm, err1R_L2Norm, err2_L2Norm, ...
            h0] = main_ECC(k_ECC, kn_max);
 
        errexa_vec(k) = errexa; 

    end

    errexa_veclog10 = log10(errexa_vec);
    cfun_errexa = fit(kn_veclog10, errexa_veclog10, 'poly1');

    conv_ord_vec(m) = cfun_errexa(1) - cfun_errexa(0);
end
toc

mu_crit = 0.002;     % h/kn_max

figure(1)
semilogx(mu_vec, conv_ord_vec, 'bo')
hold on
semilogx(mu_crit*[1,1], [1.5, 3], 'k--')
% semilogx(kn_min*[1,1], [0, 3], 'r-')
% semilogx(kn_max*[1,1], [0, 3], 'r-')
xplg = [kn_min, kn_min, kn_max, kn_max];
yplg = [1.5, 3, 3, 1.5];
col = [0.9, 0, 1];
fill(xplg, yplg, col, 'FaceAlpha', 0.4)
legend('order')
xlabel('$\mu$','interpreter','latex')
ylabel('$order \,\, of \,\, convergence$','interpreter','latex')
text(6e-4, 1.5 - 0.05, '$\mu_{sweep}$','interpreter','latex')
% text(6e-4, -0.14, '$k_{n,min}$','interpreter','latex')
% text(6e-4, -0.14, '$k_{n,max}$','interpreter','latex')
text(1e-1, 1.5 -0.05, '$k_n$','interpreter','latex')

PrintCOVM(mu_vec, conv_ord_vec)