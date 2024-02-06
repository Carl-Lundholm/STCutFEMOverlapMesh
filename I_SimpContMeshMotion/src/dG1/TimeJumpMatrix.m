function A_tjnm1 = TimeJumpMatrix(K0, KG, x0, G, anm1, bnm1)

global x0_init x0_fin I0 M

A_k = zeros(2);
A_tjnm1 = zeros(M);

% Regular domains of U1 ---------------------------------------------------

% To the left of G
K1a = K0(:, K0(2,:) <= anm1);
leK1a = length(K1a(1,:));
for k = 1:leK1a

    x_k = K1a(2,k);
    x_kpos = find((x0 == x_k));
    kpos = (x_kpos - 1: x_kpos);

    hk = K1a(3, k);
    A_k = hk/3*[1, 0.5; 0.5, 1];

    A_tjnm1(kpos, kpos) = A_tjnm1(kpos, kpos) + A_k;

end

% To the right of G
K1b = K0(:, K0(1,:) >= bnm1);
leK1b = length(K1b(1,:));
for k = 1:leK1b

    x_k = K1b(2,k);
    x_kpos = find((x0 == x_k));
    kpos = (x_kpos - 1: x_kpos);

    hk = K1b(3, k);
    A_k = hk/3*[1, 0.5; 0.5, 1];

    A_tjnm1(kpos, kpos) = A_tjnm1(kpos, kpos) + A_k;

end

% Boundary domains of U1 --------------------------------------------------

% To the left of G
K1ba = K0(:, K0(1,:) < anm1 & anm1 < K0(2,:));
leK1ba = length(K1ba(1,:));
if leK1ba > 0
    x_km1 = K1ba(1,1);
    x_k = K1ba(2,1);
    x_kpos = find((x0 == x_k));
    kpos = (x_kpos - 1: x_kpos);

    hk = K1ba(3, 1);
    A_k(1,1) = 1/3*(hk - (x_k - anm1)^3/hk^2);
    A_k(1,2) = 1/(6*hk^2)*(anm1 - x_km1)^2*(3*x_k - 2*anm1 - x_km1);
    A_k(2,1) = A_k(1,2);
    A_k(2,2) = (anm1 - x_km1)^3/(3*hk^2);

    A_tjnm1(kpos, kpos) = A_tjnm1(kpos, kpos) + A_k;
end


% To the right of G
K1bb = K0(:, K0(1,:) < bnm1 & bnm1 < K0(2,:));
leK1bb = length(K1bb(1,:));
if leK1bb > 0
    x_km1 = K1bb(1,1);
    x_k = K1bb(2,1);
    x_kpos = find((x0 == x_k));
    kpos = (x_kpos - 1: x_kpos);

    hk = K1bb(3, 1);
    A_k(1,1) = 1/(3*hk^2)*(x_k - bnm1)^3;
    A_k(1,2) = 1/6*(hk - (bnm1 - x_km1)^2*(3*x_k - 2*bnm1 - x_km1)/hk^2);
    A_k(2,1) = A_k(1,2);
    A_k(2,2) = 1/3*(hk - (bnm1 - x_km1)^3/hk^2);

    A_tjnm1(kpos, kpos) = A_tjnm1(kpos, kpos) + A_k;
end


% Domains of U2 -----------------------------------------------------------
K2 = KG(:, x0_init <= KG(1,:) & KG(2,:) <= x0_fin);
leK2 = length(K2(1,:));
for k = 1:leK2

    x_k = K2(2,k);
    x_kpos = find((G == x_k));
    kpos = I0 + (x_kpos - 1: x_kpos);

    hk = K2(3, k);
    A_k = hk/3*[1, 0.5; 0.5, 1];

    A_tjnm1(kpos, kpos) = A_tjnm1(kpos, kpos) + A_k;

end
% figure(1)
% test = [K1a(1:2,:), K1b(1:2,:), K1ba(1:2,:), K1bb(1:2,:)] ;
% plot(test,t*ones(2, length(test(1,:))),'go-');
% hold on
