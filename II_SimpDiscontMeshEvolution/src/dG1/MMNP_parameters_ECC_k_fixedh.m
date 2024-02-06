% Fixed step size
% h_fix = 5e-4; % For kn = 1e-3 this takes 1-2 h on my Mac. 
% h_fix = 1e-4; TOO SMALL, takes too much time.
h_fix = 5e-4; % OBS! This takes a long time. Only run on ozzy
% h_fix = 5e-1;

% Background mesh
x0_init = 0;
x0_fin = 1;
I0 = ceil((x0_fin - x0_init)/h_fix) + 1;
h0 = (x0_fin - x0_init)/(I0 - 1);

% Moving mesh
L = 0.25;
a0 = 0.125;
a = a0 + 0.5*h0*(1 - ceil(ceil(a0/h0) - a0/h0));
b = a + L;
IG = ceil(L/h_fix) + 1;
hG = L/(IG - 1);

% Initial velocity
mutr = mu_iter;   

% Time
t0 = 0;
T = 1;
N = ceil((T-t0)/k_ECC) + 1;
kn = (T-t0)/(N-1);

% Boundary condition parameter
eta = 1e5;

% Weights, w1 + w2 = 1 !
w1 = 0.5;
w2 = 1 - w1;           

% Penalty parameters
pnltp = 10; % >= 2*hmax*4/hmin according to a primitive error analysis
% pnltp = 1.0*8*max(h0, hG)/min(h0,hG);
pnltp_nab = 1;