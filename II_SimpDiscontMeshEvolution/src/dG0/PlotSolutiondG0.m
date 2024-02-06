function PlotSolutiondG0(U, x0, G, Gnm1, amin, amax, bmin, bmax, ...
    tna, tnm1a, tnb, tnm1b)
                    

global I0 IG tn tnm1

U0 = U(1:I0)';
U2 = U(I0 + 1:I0 + IG)';

hold on

% Regular domains of U1 ---------------------------------------------------

% To the left of G
Xind = (x0 < amin);
Xind(1) = 1;

xind_aminL = find(Xind, 1, 'last' );
xind_aminR = find(Xind, 1, 'last' ) + 1;
xind_extr = xind_aminR;
if xind_extr <= I0
    Xind(xind_extr) = 1;
end

X = x0(Xind);
leX = length(X);

if leX > 1
    Z = U0(Xind);
    Zfin = LinearInterpolationSingVal(amin, ...
        X(leX-1), X(leX), Z(leX-1), Z(leX));
    X(leX) = amin;
    Z(leX) = Zfin;
    
    space_X = [X; X];
    
    time_Y = [tnm1; tn]*ones(1,length(X));
    
    sol_Z = [Z; Z];
    
    surf(space_X, time_Y, sol_Z);
end

% To the right of G
Xind = (bmax < x0);
Xind(I0) = 1;

xind_bmaxL = find(Xind, 1, 'first' ) - 1;
xind_bmaxR = find(Xind, 1, 'first' );

xind_extr = xind_bmaxL;
if xind_extr > 0;
    Xind(xind_extr) = 1;
end

X = x0(Xind);
leX = length(X);

if leX > 1
    Z = U0(Xind);
    Zinit = LinearInterpolationSingVal(bmax, ...
        X(1), X(2), Z(1), Z(2));
    X(1) = bmax;
    Z(1) = Zinit;
    
    space_X = [X; X];
    
    time_Y = [tnm1; tn]*ones(1,length(X));
    
    sol_Z = [Z; Z];
    
    surf(space_X, time_Y, sol_Z);
end

% Gamma domains of U1 -----------------------------------------------------

% To the left of G
Xind = (amin < x0 & x0 < amax);

xind_extr = xind_aminL;
if xind_extr > 0;
    Xind(xind_extr) = 1;
end

xind_extr = find(Xind, 1, 'last' ) + 1;
if xind_extr <= I0
    Xind(xind_extr) = 1;
end

X = x0(Xind);
leX = length(X);

if leX > 1
    Z = U0(Xind);
    Zinit = LinearInterpolationSingVal(amin, ...
        X(1), X(2), Z(1), Z(2));
    Zfin = LinearInterpolationSingVal(amax, ...
        X(leX-1), X(leX), Z(leX-1), Z(leX));
    X(1) = amin;
    X(leX) = amax;
    Z(1) = Zinit;
    Z(leX) = Zfin;
    
    space_X = [X; X];
    
    time_Y = [tnm1a; tna];
    
    sol_Z = [Z; Z];
    
    surf(space_X, time_Y, sol_Z);
end

% To the right of G
Xind = (bmin < x0 & x0 < bmax);

xind_extr = xind_bmaxR;
if xind_extr <= I0
    Xind(xind_extr) = 1;
end

xind_extr = find(Xind, 1, 'first' ) - 1;
if xind_extr > 0;
    Xind(xind_extr) = 1;
end



X = x0(Xind);
leX = length(X);

if leX > 1
    Z = U0(Xind);
    Zinit = LinearInterpolationSingVal(bmin, ...
        X(1), X(2), Z(1), Z(2));
    Zfin = LinearInterpolationSingVal(bmax, ...
        X(leX-1), X(leX), Z(leX-1), Z(leX));
    X(1) = bmin;
    X(leX) = bmax;
    Z(1) = Zinit;
    Z(leX) = Zfin;
    
    space_X = [X; X];
    
    time_Y = [tnm1b; tnb];
    
    sol_Z = [Z; Z];
    
    surf(space_X, time_Y, sol_Z);
end

% Domains of U2 -----------------------------------------------------------  
space_X = [Gnm1; G]; 

time_Y = [tnm1; tn]*ones(1,IG);

sol_Z = [U2; U2];

surf(space_X, time_Y, sol_Z);

















