function u = u_func(x,t)

% Heat equation
u = (sin(pi*x)).^2.*exp(-0.5*t);

% DEBUG_errorkslope2p6
% u = x.*(1.0-x).*(1.0 - t^3 - t^2 - t);
% u = x*exp(-t);
% u = x.*(1.0-x).*exp(-2.0*t);
% u = (sin(pi*x)).^2.*exp(-t);

% CHECK h-convergence
% u = sin(pi*x).*exp(-0.5*t);


% % Heat equation + u in LHS
% u = exp(-t) + x*0;

% % Poisson
% u = (sin(pi*x)).^2;

%u = x.*(1 - x).^2.*exp(-0.5*t);