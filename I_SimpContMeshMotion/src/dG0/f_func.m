function f = f_func(x,t)

% Heat equation
ut = -0.5*(sin(pi*x)).^2.*exp(-0.5*t);
uxx = 2*pi^2*cos(2*pi*x).*exp(-0.5*t);

f = ut - uxx;

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