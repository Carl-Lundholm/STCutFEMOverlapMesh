function A_pnt = PenaltyNablaTermdG1(K0, KG, x0, mutr, C, t1ga, t1gb, ...
                                amin, amax, bmin, bmax, pnltp_nab)

global  I0 IG leA M

D = 1/C;

A_pnt = zeros(leA);

% To the left of G --------------------------------------------------------
n1 = 1;
nx = n1*D;

GbL = I0 + 1;
GbR = I0 + 2;
hGb = KG(3,1);

gmin = amin;
gmax = amax;

tg = t1ga;

A_pnt = PeNaTedG1_AuxiliaryFunc(A_pnt, K0, x0, n1, nx, tg, ...
    GbL, GbR, hGb, gmin, gmax, pnltp_nab);
  
% To the right of G -------------------------------------------------------
n1 = -1;
nx = n1*D;

GbL = M - 1;
GbR = M;
hGb = KG(3, IG - 1);

gmin = bmin;
gmax = bmax;

tg = t1gb;

A_pnt = PeNaTedG1_AuxiliaryFunc(A_pnt, K0, x0, n1, nx, tg, ...
    GbL, GbR, hGb, gmin, gmax, pnltp_nab);

  
       
 
    
    
    
    
    
    
    
    
    
    
    