function multfig2eps(pofn)

%newdirname = 'testepsfiles';
newdirname = 'epsfiles';

%set(gcf,'renderer','Painters')

mkdir(newdirname);

files = dir(['*',pofn,'*']);

for file = files'

    mlfig = file.name;
    
    open(mlfig);
    
    set(gcf,'renderer','Painters')

    xlabel('$x$', 'interpreter', 'latex', 'Fontsize', 24);
    ylabel('$t$', 'interpreter', 'latex', 'Fontsize', 24);
    zlabel('$u_h$', 'interpreter', 'latex', 'Fontsize', 24)
    
    newfilename = strrep(mlfig,'fig','eps');
    %exportgraphics(gcf,[[newdirname,'\'],newfilename],'ContentType','vector')
    %saveas(gcf,[[newdirname,'\'],newfilename],'epsc');
    print('-depsc','-painters',[[newdirname,'\'],newfilename]);   
    %print -depsc2 -tiff -r300 -painters [[newdirname,'\'],newfilename]

close

end