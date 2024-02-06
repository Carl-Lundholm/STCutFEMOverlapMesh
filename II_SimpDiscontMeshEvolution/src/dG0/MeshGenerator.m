function [x1, x2, infoOmega_2] = MeshGenerator(x0, G)

global x0_init x0_fin a b h0 hG

if (x0_init <= a && b <= x0_fin)
    x1 = x0(x0 < a + h0 | b - h0 < x0);
    x2 = G;
    infoOmega_2 = 'Omega_2 helt i det inre.';
    
elseif (a < x0_init && x0_fin < b)
    x1 = zeros(1,0);
    x2 = G(x0_init - hG < G & G < x0_fin + hG);
    infoOmega_2 = 'Omega_2 = Omega_0';
    
elseif (a < x0_init && x0_init < b)
    x1 = x0(b - h0 < x0);
    x2 = G(x0_init - hG < G);
    infoOmega_2 = 'Omega_2 delvis i det inre.';
    
elseif (a < x0_fin && x0_fin < b)
    x1 = x0(x0 < a + h0);
    x2 = G(G < x0_fin + hG);
    infoOmega_2 = 'Omega_2 delvis i det inre.';
    
else
    x1 = x0;
    x2 = zeros(1,0);
    infoOmega_2 = 'Omega_2 tom';
    
end