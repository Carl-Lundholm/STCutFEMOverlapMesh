function errL2Norm = ErrorL2NormSpaceQuadTrapezoidal(U, Ua, x0, G, a, b)

% The trapezoidal rule is used locally to approximate the error integrals
% thus yielding a quadrature error of the second order, i.e., h^2.

global I0 IG leA

U0a = Ua(1:I0)';
U2a = Ua(I0 + 1:leA)';

U0 = U(1:I0)';
U2 = U(I0 + 1:leA)';

errL2Norm = 0;

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
U1aL = U0a(Xind);

leL = length(x1L);

U1L(leL) = LinearInterpolationSingVal(a, x1L(leL - 1), x1L(leL), ...
                                         U1L(leL - 1), U1L(leL));
U1aL(leL) = LinearInterpolationSingVal(a, x1L(leL - 1), x1L(leL), ...
                                         U1aL(leL - 1), U1aL(leL));
x1L(leL) = a;

err1L = abs(U1aL - U1L).^2;

delx1L = (x1L(2:leL) - x1L(1:leL - 1));
err1L_av = 0.5*(err1L(2:leL) + err1L(1:leL - 1));

err1L_L2Norm_vec = err1L_av.*delx1L;
err1L_L2Norm = sum(err1L_L2Norm_vec);

errL2Norm = errL2Norm + err1L_L2Norm;

% plot(x1L, err1L, 'k-o');

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
U1aR = U0a(Xind);

leR = length(x1R);

U1R(1) = LinearInterpolationSingVal(b, x1R(1), x1R(2), U1R(1), U1R(2));
U1aR(1) = LinearInterpolationSingVal(b, x1R(1), x1R(2), U1aR(1), U1aR(2));
x1R(1) = b;

err1R = abs(U1aR - U1R).^2;

delx1R = (x1R(2:leR) - x1R(1:leR - 1));
err1R_av = 0.5*(err1R(2:leR) + err1R(1:leR - 1));

err1R_L2Norm_vec = err1R_av.*delx1R;
err1R_L2Norm = sum(err1R_L2Norm_vec);

errL2Norm = errL2Norm + err1R_L2Norm;

% plot(x1R, err1R, 'k-o');

% Domains of U2 -----------------------------------------------------------  

err2 = abs(U2a - U2).^2;

delx2 = (G(2:IG) - G(1:IG - 1));
err2_av = 0.5*(err2(2:IG) + err2(1:IG - 1));

err2_L2Norm_vec = err2_av.*delx2;
err2_L2Norm = sum(err2_L2Norm_vec);

errL2Norm = errL2Norm + err2_L2Norm;

errL2Norm = sqrt(errL2Norm);

% plot(G, err2, 'k-x');