function PrintError_kn(knvec, errvec)

len = length(knvec);

for j = 1:len
 disp(['kn_vec(', num2str(j),') = ',num2str(knvec(j),'%17.15f'), ...
     ', errexa_vec(',num2str(j),') = ', num2str(errvec(j),'%17.15f')])
end


