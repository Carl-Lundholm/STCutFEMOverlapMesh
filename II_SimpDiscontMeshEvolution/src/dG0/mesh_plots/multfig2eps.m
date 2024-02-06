function multfig2eps(pofn)

%newdirname = 'testepsfiles';
newdirname = 'epsfiles_nummatpub';

%set(gcf,'renderer','Painters')

mkdir(newdirname);

files = dir(['*',pofn,'*']);

for file = files'

    mlfig = file.name;
    
    open(mlfig);

    lines = findobj(gcf,'Type','Line');
    for i = 1:numel(lines)
        lines(i).LineWidth = 2.0;
    end
    
    set(gcf,'renderer','Painters')

    xlabel('$x$', 'interpreter', 'latex', 'Fontsize', 24);
    ylabel('$t$', 'interpreter', 'latex', 'Fontsize', 24);
    zlabel('$u_h$', 'interpreter', 'latex', 'Fontsize', 24);
    
    newfilename = strrep(mlfig,'fig','eps');
    %exportgraphics(gcf,[[newdirname,'\'],newfilename],'ContentType','vector')
    %saveas(gcf,[[newdirname,'\'],newfilename],'epsc');
    print('-depsc','-painters',[[newdirname,'\'],newfilename]);   
    %print -depsc2 -tiff -r300 -painters [[newdirname,'\'],newfilename]

close

end