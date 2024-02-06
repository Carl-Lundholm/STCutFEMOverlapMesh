function PrintCOVM(mu_vec, conv_ord_vec)

len = length(mu_vec);

for j = 1:len
 disp(['mu_vec(', num2str(j),') = ', num2str(mu_vec(j),'%17.15f'), ...
     ', conv_ord_vec(',num2str(j),') = ', ...
     num2str(conv_ord_vec(j),'%17.15f')])
end


