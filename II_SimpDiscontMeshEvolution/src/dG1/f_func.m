function f = f_func(x,t)

% Heat equation
ut = -0.5*(sin(pi*x)).^2.*exp(-0.5*t);
uxx = 2*pi^2*cos(2*pi*x).*exp(-0.5*t);

f = ut - uxx;

% DEBUG_errorkslope2p6
% f = x.*(1.0-x).*(-3*t^2 -2*t - 1.0) + 2*(1.0 - t^3 - t^2 - t);
%f = -exp(-t) + (x-x);
%f = -0.5*x.*(1.0-x).*exp(-0.5*t) + 2*exp(-0.5*t);
% ut = -(sin(pi*x)).^2.*exp(-t);
% uxx = 2*pi^2*cos(2*pi*x).*exp(-t);
% 
% f = ut - uxx;


% CHECK h-convergence
% u = sin(pi*x).*exp(-0.5*t);
% ut = -0.5*sin(pi*x).*exp(-0.5*t);
% uxx = -pi^2*sin(pi*x).*exp(-0.5*t);
% 
% f = ut - uxx;

% % Heat equation + u in LHS
% ut = 0*-1*exp(-t);
% uxx = x*0;
% 
% f = ut - uxx;

% % Poisson
% f = -2*pi^2*cos(2*pi*x);


% ut = -0.5*x.*(1 - x).^2.*exp(-0.5*t);
% uxx = (-4 + 6*x).*exp(-0.5*t);
% 
% f = ut - uxx;
% f = abs((x.^2 - 2).*exp(-t) - x); 