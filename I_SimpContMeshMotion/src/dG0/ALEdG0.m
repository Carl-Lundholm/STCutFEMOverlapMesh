function A_ale = ALEdG0(mutr, kn)

global I0 IG 

A_ale = zeros(I0 + IG);

A_k = 0.5*mutr*kn*[1 -1; 1 -1];

for k = I0 + (2:IG) 

    A_ale(k-1:k,k-1:k) = A_ale(k-1:k,k-1:k) + A_k;
    
end
