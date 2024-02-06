% Background mesh
x0_init = 0;
x0_fin = 1;
I0 = 22;
h0 = (x0_fin - x0_init)/(I0 - 1);

% Moving mesh
L = 0.25;
a = x0_init + 0.1;
b = a + L;
IG = 7;
hG = L/(IG - 1);

% Initial velocity
mutr = 1*0.1    

% Time
t0 = 0;
T = 3;
kn = 3e-1; %*min(h0, hG)/mutr;   % < h0/mutr ?
pause_time = 3*kn;

% Boundary condition parameter
eta = 1e5;

% Weights, w1 + w2 = 1 !
w1 = 0.5;
w2 = 1 - w1;         

% Penalty parameters
pnltp = 10; % >= 2*hmax*4/hmin according to a primitive error analysis
% pnltp = 1.0*8*max(h0, hG)/min(h0,hG);
pnltp_nab = 1;