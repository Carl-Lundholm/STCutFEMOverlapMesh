function [u, u_x, u_t] = u_func(x,t)

% Heat equation
u = (sin(pi*x)).^2.*exp(-0.5*t);

u_x = 2*pi*sin(pi*x).*cos(pi*x).*exp(-0.5*t);

u_t = -0.5*u;

% % Heat equation + u in LHS
% u = exp(-t) + x*0;

% % Poisson
% u = (sin(pi*x)).^2;

%u = x.*(1 - x).^2.*exp(-0.5*t);