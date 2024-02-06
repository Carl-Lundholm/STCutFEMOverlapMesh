function PlotSolutionAtCurrentTime(U, x0, G, a, b, color1, color2)

global I0 IG

U0 = U(1:I0)';
U2 = U(I0 + 1:I0 + IG)';

% Domains of U1 ---------------------------------------------------

% To the left of G
Xind = (x0 < a);
Xind(1) = 1;
xind_aR = find(Xind, 1, 'last' ) + 1;

xind_extr = xind_aR;
if xind_extr <= I0
    Xind(xind_extr) = 1;
end

x1L = x0(Xind);
U1L = U0(Xind);

leL = length(x1L);

U1L(leL) = LinearInterpolationSingVal(a, x1L(leL - 1), x1L(leL), ...
                                         U1L(leL - 1), U1L(leL));
x1L(leL) = a;

plot(x1L, U1L, color1);

% To the right of G
Xind = (x0 > b);
Xind(I0) = 1;
xind_bL = find(Xind, 1, 'first' ) - 1;

xind_extr = xind_bL;
if xind_extr > 0
    Xind(xind_extr) = 1;
end

x1R = x0(Xind);
U1R = U0(Xind);

U1R(1) = LinearInterpolationSingVal(b, x1R(1), x1R(2), U1R(1), U1R(2));
x1R(1) = b;

plot(x1R, U1R, color1);

% Domains of U2 -----------------------------------------------------------  

plot(G, U2, color2);

















