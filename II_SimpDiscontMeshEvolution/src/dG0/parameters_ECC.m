% Background mesh
x0_init = 0;
x0_fin = 1;
I0 = 1e0*200;
h0 = (x0_fin - x0_init)/(I0 - 1);

% Moving mesh
L = 0.3;
a = x0_init + 10.5*h0;
b = a + L;
IG = 1e0*120;
hG = L/(IG - 1);

% Initial velocity
mutr = 0.0;    %1*(x0_fin - x0_init)/(T - t0);

% Time
t0 = 0;
T = 3;

% Boundary condition parameter
eta = 1e5;

% Weights, w1 + w2 = 1 !
w1 = 0.5;
w2 = 1 - w1;           

% Penalty parameters
pnltp = 1*1e0; % >= 2*hmax*4/hmin according to a primitive error analysis
% pnltp = 1.0*8*max(h0, hG)/min(h0,hG);
pnltp_nab = 1*1e0;