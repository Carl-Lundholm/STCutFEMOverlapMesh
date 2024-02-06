function u = u_func(x,t)

% Heat equation
u = (sin(pi*x)).^2.*exp(-0.5*t);

% % Heat equation + u in LHS
% u = exp(-t) + x*0;

% % Poisson
% u = (sin(pi*x)).^2;

% u = x.*(1 - x).^2.*exp(-0.5*t);
% u = x.*(1 - x).*exp(-0.5*t);