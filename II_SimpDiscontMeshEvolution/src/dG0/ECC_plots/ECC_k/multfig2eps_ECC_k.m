function multfig2eps_ECC_k(pofn)

% error_str = '$\| e(T) \|_{L^2(\Omega_0)}$';
newdirname = 'epsfiles';

mkdir(newdirname);

files = dir(['*',pofn,'*']);

for file = files'

    mlfig = file.name;
    
    open(mlfig);
    
    set(gcf,'renderer','Painters')

%     xlabel('$k$', 'interpreter', 'latex', 'Fontsize', 24);
%     ylabel(error_str, 'interpreter', 'latex', 'Fontsize', 24);
%     legend(error_str,'LLS','$k^1$','$k^2$','$k^3$', ...
%         'interpreter','latex', 'Fontsize', 24, 'Location','SouthEast');
    
    newfilename = strrep(mlfig,'fig','eps');
    print('-depsc','-painters',[[newdirname,'\'],newfilename]);
    
close

end