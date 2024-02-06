function A_pnt = PenaltyNablaTermdG0_BE(K0, KG, x0, t1ga, t1gb, ...
                                xGl, xGr, pnltp_nab)

global  I0 IG leA

A_pnt = zeros(leA);

% To the left of G --------------------------------------------------------
n1 = 1;
nx = n1;

GbL = I0 + 1;
GbR = I0 + 2;
hGb = KG(3,1);

tg = t1ga;

A_pnt = PeNaTedG0_AuxiliaryFunc(A_pnt, K0, x0, n1, nx, tg, ...
    GbL, GbR, hGb, xGl, xGl, pnltp_nab);
  
% To the right of G -------------------------------------------------------
n1 = -1;
nx = n1;

GbL = leA - 1;
GbR = leA;
hGb = KG(3, IG - 1);

tg = t1gb;

A_pnt = PeNaTedG0_AuxiliaryFunc(A_pnt, K0, x0, n1, nx, tg, ...
    GbL, GbR, hGb, xGr, xGr, pnltp_nab);

  
       
 
    
    
    
    
    
    
    
    
    
    
    