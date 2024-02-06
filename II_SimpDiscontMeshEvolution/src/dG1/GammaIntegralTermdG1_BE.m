function A_git = GammaIntegralTermdG1_BE(K0, KG, x0, t1ga, t1gb, ...
                                xGl, xGr, w1, w2, pnltp)

global  I0 IG leA M

A_git = zeros(leA);

% To the left of G --------------------------------------------------------
n1 = 1;
nt = 0;     
nx = -n1;   % The extra minus sign from the FE-form is included  

GbL = I0 + 1;
GbR = I0 + 2;
Gb = GbL;
hGb = KG(3,1);

tg = t1ga;

Kgs = K0(:, (K0(1,:) < xGl & xGl <= K0(2,:)));

A_git = GaInTedG1_AuxiliaryFunc(A_git, Kgs, x0, tg, nt, nx, ...
    GbL, GbR, Gb, hGb, xGl, xGl, w1, w2, pnltp);
  
% To the right of G -------------------------------------------------------
n1 = -1;
nt = 0;     
nx = -n1;   % The extra minus sign from the FE-form is included 

GbL = M - 1;
GbR = M;
Gb = GbR;
hGb = KG(3, IG - 1);

tg = t1gb;

Kgs = K0(:, (K0(1,:) <= xGr & xGr < K0(2,:)));

A_git = GaInTedG1_AuxiliaryFunc(A_git, Kgs, x0, tg, nt, nx, ...
    GbL, GbR, Gb, hGb, xGr, xGr, w1, w2, pnltp);
  
       
 
    
    
    
    
    
    
    
    
    
    
    