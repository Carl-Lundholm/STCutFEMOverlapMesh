function PlotSolutiondG1(U, Unm1, x0, G, Gnm1, amin, amax, bmin, bmax, ...
    tna, tnm1a, tnb, tnm1b)

hold on

global I0 IG tn tnm1 M

U0n = U(1:I0)';
U2n = U(I0 + 1:M)';

U0nm1 = Unm1(1:I0)';
U2nm1 = Unm1(I0 + 1:M)';

% Regular domains of U1 ---------------------------------------------------

% To the left of G
Xind = (x0 <= amin);
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
    Znm1 = U0nm1(Xind);
    Zfinnm1 = LinearInterpolationSingVal(amin, ...
        X(leX-1), X(leX), Znm1(leX-1), Znm1(leX));
    Znm1(leX) = Zfinnm1;
    
    Zn = U0n(Xind);
    Zfinn = LinearInterpolationSingVal(amin, ...
        X(leX-1), X(leX), Zn(leX-1), Zn(leX));  
    Zn(leX) = Zfinn;
    
    X(leX) = amin;
    
    space_X = [X; X];
    
    time_Y = [tnm1; tn]*ones(1,length(X));
    
    sol_Z = [Znm1; Zn];
    
    surf(space_X, time_Y, sol_Z);
end

% To the right of G
Xind = (bmax <= x0);
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
    Znm1 = U0nm1(Xind);
    Zinitnm1 = LinearInterpolationSingVal(bmax, ...
        X(1), X(2), Znm1(1), Znm1(2));
    Znm1(1) = Zinitnm1;
    
    Zn = U0n(Xind);
    Zinitn = LinearInterpolationSingVal(bmax, ...
        X(1), X(2), Zn(1), Zn(2));
    Zn(1) = Zinitn;
    
    X(1) = bmax;
    
    space_X = [X; X];
    
    time_Y = [tnm1; tn]*ones(1,length(X));
    
    sol_Z = [Znm1; Zn];
    
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
    Znm1 = U0nm1(Xind);
    Zinitnm1 = LinearInterpolationSingVal(amin, ...
        X(1), X(2), Znm1(1), Znm1(2));
    Zfinnm1 = LinearInterpolationSingVal(amax, ...
        X(leX-1), X(leX), Znm1(leX-1), Znm1(leX));
    Znm1(1) = Zinitnm1;
    Znm1(leX) = Zfinnm1;
    
    Zn = U0n(Xind);
    Zinitn = LinearInterpolationSingVal(amin, ...
        X(1), X(2), Zn(1), Zn(2));
    Zfinn = LinearInterpolationSingVal(amax, ...
        X(leX-1), X(leX), Zn(leX-1), Zn(leX));
    Zn(1) = Zinitn;
    Zn(leX) = Zfinn;
    
    X(1) = amin;
    X(leX) = amax;
    
    if tnm1a(1) ~= tnm1a(2)
        for t = 2:leX
            Znm1(t) = LinearInterpolationSingVal(tnm1a(t), ...
            tnm1, tn, Znm1(t), Zn(t));
        end
    end
    
    if tna(1) ~= tna(2)
        for t = 2:leX
            Zn(t) = LinearInterpolationSingVal(tna(t), ...
            tnm1, tn, Znm1(t), Zn(t));  
        end
    end
    
    space_X = [X; X];
    
    time_Y = [tnm1a; tna];
    
    sol_Z = [Znm1; Zn];
    
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
    Znm1 = U0nm1(Xind);
    Zinitnm1 = LinearInterpolationSingVal(bmin, ...
        X(1), X(2), Znm1(1), Znm1(2));
    Zfinnm1 = LinearInterpolationSingVal(bmax, ...
        X(leX-1), X(leX), Znm1(leX-1), Znm1(leX));
    Znm1(1) = Zinitnm1;
    Znm1(leX) = Zfinnm1;
    
    Zn = U0n(Xind);
    Zinitn = LinearInterpolationSingVal(bmin, ...
        X(1), X(2), Zn(1), Zn(2));
    Zfinn = LinearInterpolationSingVal(bmax, ...
        X(leX-1), X(leX), Zn(leX-1), Zn(leX));
    Zn(1) = Zinitn;
    Zn(leX) = Zfinn;
    
    X(1) = bmin;
    X(leX) = bmax;
    
     if tnm1b(1) ~= tnm1b(2)
        for t = 1:leX-1
            Znm1(t) = LinearInterpolationSingVal(tnm1b(t), ...
            tnm1, tn, Znm1(t), Zn(t));
        end
    end
    
    if tnb(1) ~= tnb(2)
        for t = 1:leX-1
            Zn(t) = LinearInterpolationSingVal(tnb(t), ...
            tnm1, tn, Znm1(t), Zn(t));  
        end
    end
    
    space_X = [X; X];
    
    time_Y = [tnm1b; tnb];
    
    sol_Z = [Znm1; Zn];
    
    surf(space_X, time_Y, sol_Z);
end

% Domains of U2 -----------------------------------------------------------  
space_X = [Gnm1; G]; 

time_Y = [tnm1; tn]*ones(1,IG);

sol_Z = [U2nm1; U2n];

surf(space_X, time_Y, sol_Z);


















