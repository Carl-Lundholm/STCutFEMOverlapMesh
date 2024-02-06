function PlotMesh(G, x1, x2, x1r, x1ga, x1gb, x2_nm1, ...
                  tn, tnm1, tna, tnm1a, tnb, tnm1b)

global I1 I2

xGl = G(1);
xGr = G(end);

plot(x1,tn*ones(I1,1),'bo-');
hold on
if ~isempty(x2)
    plot(x2,tn*ones(I2,1),'rx-');
end
if ~isempty(x1r)
    plot([x1r; x1r],[tn, tnm1],'b-')
end
if ~isempty(x1ga)
    plot([x1ga; x1ga], [tna(2:end-1); tnm1a(2:end-1)],'b-')
end
if ~isempty(x1gb)
    plot([x1gb; x1gb], [tnb(2:end-1); tnm1b(2:end-1)],'b-')
end

plot([G; G],[tn, tnm1],'r-')
plot([xGl, xGl],[tn, tnm1],'k-')
plot([xGr, xGr],[tn, tnm1],'k-')




    