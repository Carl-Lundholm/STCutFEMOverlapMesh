function PrintError_h(hvec, errvec)

len = length(hvec);

for j = 1:len
 disp(['h_vec(', num2str(j),') = ',num2str(hvec(j),'%17.15f'), ...
     ', errexa_vec(',num2str(j),') = ', num2str(errvec(j),'%17.15f')])
end


