function error = ErrorXNormSlabwiseSquared_AuxiliaryGamma(Un, Unm1p, u)

global mutr K0 KG amin amax bmin bmax t1ga t1gb I0 M C

error = 0;

if mutr == 0 % amin = amax, bmin = bmax
    % Special conditions for left and right of G when mu = 0
    Kgs_left = K0(:, K0(1,:) < amin & amax <= K0(2,:));
    Kgs_right = K0(:, K0(1,:) <= bmin & bmax < K0(2,:));
else 
    % General conditions for left and right of G
    Kgs_left = K0(:, K0(1,:) < amax & amin < K0(2,:));
    Kgs_right = K0(:, K0(1,:) < bmax & bmin < K0(2,:));
end

D = 1/C;

% Absolute value of time part of space-time vector on current slab
nt_abs = abs(mutr/C);

% Absolute value of space part of space-time vector on current slab
nx_abs = D;

% Left part of Gamma
Kgs = Kgs_left;
gmin = amin;
gmax = amax;
t1gg = t1ga;
indG = I0 + 1;
indGL = I0 + 1;
indGR = I0 + 2;
hGb = KG(3, 1);

eL = ErrorXNormSlabwiseSquared_AuxiliaryGammaGen(Un, Unm1p, ...
                 u, nt_abs, nx_abs, Kgs, gmin, gmax, t1gg, ...
                 indG, indGL, indGR, hGb);
             
error = error + eL;

% Right part of Gamma
Kgs = Kgs_right;
gmin = bmin;
gmax = bmax;
t1gg = t1gb;
indG = M;
indGL = M - 1;
indGR = M;
hGb = KG(3, end); 

eR = ErrorXNormSlabwiseSquared_AuxiliaryGammaGen(Un, Unm1p, ...
                 u, nt_abs, nx_abs, Kgs, gmin, gmax, t1gg, ...
                 indG, indGL, indGR, hGb);
             
error = error + eR;  

end

