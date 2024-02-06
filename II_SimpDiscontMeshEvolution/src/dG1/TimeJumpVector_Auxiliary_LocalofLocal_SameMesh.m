function b_kloc = TimeJumpVector_Auxiliary_LocalofLocal_SameMesh(a, b, ...
                    U_lm1, U_l, K)

% Computes the local tnm1_minus contribution from interval (a, b) on 
% simplex K to the local tnm1_plus element vector corresponding to the 
% same simplex K.

% Note that this m-file is a special case and thus a somewhat cheaper 
% version of the m-file "TimeJumpVector_Auxiliary_LocalofLocal_DiffMesh".
                                
Pkm1km1 = PHI_km1km1(a, b, K);
Pkm1k = PHI_km1k(a, b, K);
Pkk = PHI_kk(a, b, K);
     
b_kloc = [U_lm1*Pkm1km1 + U_l*Pkm1k; U_lm1*Pkm1k + U_l*Pkk];              
                
                

