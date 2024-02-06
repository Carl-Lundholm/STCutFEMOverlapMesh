function PlotMesh(x1, x2, x1r, x1o, x1ga, x1gb, x2_nm1, ...
                  tn, tnm1, tna, tnm1a, tnb, tnm1b, tnao, tnm1ao, tnbo, tnm1bo)

global a b anm1 bnm1 I1 I2 x0_init x0_fin t0 G Gnm1

plot(x1,tn*ones(I1,1),'bo-');
hold on
if ~isempty(x2)
    plot(x2,tn*ones(I2,1),'rx-');
end
if ~isempty(x1r)
    plot([x1r; x1r],[tn, tnm1],'b-')
end
if ~isempty(x1o)
    plot([x1o; x1o],[tn, tnm1],'b:')
end
if ~isempty(x1ga)
    plot([x1ga; x1ga], [tna(2:end-1); tnm1a(2:end-1)],'b-')
end
if ~isempty(x1gb)
    plot([x1gb; x1gb], [tnb(2:end-1); tnm1b(2:end-1)],'b-')
end
if ~isempty(x1ga)
    plot([x1ga; x1ga], [tnao(2:end-1); tnm1ao(2:end-1)],'b:')
end
if ~isempty(x1gb)
    plot([x1gb; x1gb], [tnbo(2:end-1); tnm1bo(2:end-1)],'b:')
end


plot([G; Gnm1],[tn, tnm1],'r-')
plot([a, anm1],[tn, tnm1],'k-')
plot([b, bnm1],[tn, tnm1],'k-')




    