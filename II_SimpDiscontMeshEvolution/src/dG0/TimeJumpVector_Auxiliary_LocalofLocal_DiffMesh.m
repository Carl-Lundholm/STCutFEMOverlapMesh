function b_kloc = TimeJumpVector_Auxiliary_LocalofLocal_DiffMesh(a, b, ...
                    U_lm1, U_l, Kp, Km)

% Computes the local tnm1_minus contribution from interval (a, b) on 
% simplex Km (belonging to mesh composition at time tnm1_minus) to the
% local tnm1_plus element vector corresponding to simplex Kp (belonging to
% mesh composition at time tnm1_plus).

% Note that this m-file is a generalized and somewhat more expensive 
% version of the m-file "TimeJumpVector_Auxiliary_LocalofLocal_SameMesh".
                
Pkm1lm1 = PHI_km1lm1(a, b, Kp, Km);
Pkm1l = PHI_km1l(a, b, Kp, Km);
Pklm1 = PHI_klm1(a, b, Kp, Km);
Pkl = PHI_kl(a, b, Kp, Km);

b_kloc = [U_lm1*Pkm1lm1 + U_l*Pkm1l; U_lm1*Pklm1 + U_l*Pkl];