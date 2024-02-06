% Fixed time step
kn = 1e-4;

% Background mesh
x0_init = 0;
x0_fin = 1;
I0 = ceil((x0_fin - x0_init)/h_ECC) + 1;
h0 = (x0_fin - x0_init)/(I0 - 1);

% Moving mesh
L = 0.25;
a0 = 0.125;
a = a0 + 0.5*h0*(1 - ceil(ceil(a0/h0) - a0/h0));
b = a + L;
IG = ceil(L/h_ECC) + 1;
hG = L/(IG - 1);

% Initial velocity
mutr = mu_iter;

% Time
t0 = 0;
T = 1;
N = (T-t0)/kn + 1;

% Boundary condition parameter
eta = 1e5;

% Weights, w1 + w2 = 1 !
w1 = 0.5;
w2 = 1 - w1;           

% Penalty parameters
pnltp = 10; % >= 2*hmax*4/hmin according to a primitive error analysis
% pnltp = 1.0*8*max(h0, hG)/min(h0,hG);
pnltp_nab = 1;