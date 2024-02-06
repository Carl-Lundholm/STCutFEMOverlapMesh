function multfig2eps_ECC_h(pofn)

% error_str = '$\left|\mkern-1.5mu\left|\mkern-1.5mu\left| e \right|\mkern-1.5mu\right|\mkern-1.5mu\right|_{X}$';
newdirname = 'epsfiles';

mkdir(newdirname);

files = dir(['*',pofn,'*']);

for file = files'

    mlfig = file.name;
    
    open(mlfig);
    
    set(gcf,'renderer','Painters')

%     xlabel('$h$', 'interpreter', 'latex', 'Fontsize', 24);
%     ylabel(error_str, 'interpreter', 'latex', 'Fontsize', 24);
%     legend(error_str,'LLS','$h^{0.5}$','$h^{1.0}$','$h^{1.5}$', ...
%         'interpreter','latex', 'Fontsize', 24, 'Location','SouthEast');
    
    newfilename = strrep(mlfig,'fig','eps');
    print('-depsc','-painters',[[newdirname,'\'],newfilename]); 

close

end