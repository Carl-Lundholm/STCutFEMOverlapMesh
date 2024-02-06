% Background mesh
x0_init = 0;
x0_fin = 1;
I0 = 10; %44;
h0 = (x0_fin - x0_init)/(I0 - 1);

% Moving mesh
L = 0.25;
a = x0_init + 0.1;
b = a + L;
IG = 5; %14;
hG = L/(IG - 1);

% Initial velocity
mutr = 5e-1; %1e-1;   

% Time
t0 = 0;
T = 8;  %2.215;
kn = 0.5*min(h0,hG)/mutr;   % < h0/mutr ?
pause_time = kn; % 0.5;

% Boundary condition parameter
eta = 1e5;

% Weights, w1 + w2 = 1 !
w1 = 0.5;
w2 = 1 - w1;           

% Penalty parameters 
% pnltp = 1*1e-3; % >= 2*hmax*4/hmin according to a primitive error analysis
pnltp = 8*max(h0, hG)/min(h0,hG);
pnltp_nab = 1*1e0;


