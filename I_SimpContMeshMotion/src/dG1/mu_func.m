function mutr = mu_func(tn, mutr, a, b)

global x0_init x0_fin hG kn

% % mutr = 0.2;
% mutr = -0.215 * (1.5 - tn);
% mutr = 0.5*sin(pi*tn/1.5);

% Faked Bouncing G ------------------------------------------------------

sfac = 1.0;
absmu = abs(mutr);

if a < x0_init + sfac*absmu*kn 
    mutr = absmu;
    
elseif x0_fin - sfac*absmu*kn < b
    mutr = - absmu;
    
end

% % Faked Bouncing G with decreasing energy, G must start in the interior -
% 
% sfac = 1.0;
% enfac = 1.0;
% if x0_fin - sfac*abs(mutr)*kn < b || ...
%             a < x0_init + sfac*abs(mutr)*kn 
%     mutr = -enfac*mutr;
%
% end


end