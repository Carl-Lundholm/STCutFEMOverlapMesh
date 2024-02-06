function multfig2fususide(pofn)

files = dir(['*',pofn,'*']);

for file = files'

    mlfig = file.name;
    
    open(mlfig);
    
    % View settings for the fususide figures
    az = 173;
    el = 22;
    axis([0 1 0 3 0 1])
    view(az,el)
    
    newfilename = strrep(mlfig,'fusuab','fususide');
    saveas(gcf, newfilename);
   
    close

end



