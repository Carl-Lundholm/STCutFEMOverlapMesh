global kn kn_min kn_max mu_COVM

% Background mesh
x0_init = 0;
x0_fin = 1;
h0 = 0.001; 
I0 = (x0_fin - x0_init)/h0 + 1;

% Moving mesh
L = 0.3;
a = x0_init + 0.5*h0;
b = a + L;
hG = 0.001; 
IG = L/hG + 1;

% Initial velocity
mutr = mu_COVM;    %1*(x0_fin - x0_init)/(T - t0);

% Time
t0 = 0;
% T = 1e-0; 
T = 1;

% Boundary condition parameter
eta = 1e5;

% Weights, w1 + w2 = 1 !
w1 = 0.5;
w2 = 1 - w1;           

% Penalty parameters
pnltp = 1*1e1; % >= 2*hmax*4/hmin according to a primitive error analysis
% pnltp = 1.0*8*max(h0, hG)/min(h0,hG);
pnltp_nab = 1*1e1;