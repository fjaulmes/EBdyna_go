function values = lininterp3(map,indexes,slopes,expr)

if nargin==3
    expr=true(numel(slopes.x),1);
end
% Linear interpolation
sx=slopes.x(expr);
sy=slopes.y(expr);
sz=slopes.z(expr);
sxsy=sx.*sy;
sxsz=sx.*sz;
% sysz=sy.*sz;
sxsysz=sxsy.*sz;
ds1=sxsysz-sy.*sz;
ds2=sxsysz-sxsz;
ds3=sxsz-sz;
ds4=sxsy-sy;

values=      (1-sx+ds4+ds3-ds1).*map(indexes.ldb(expr))...   %left  down  bottom corner
            +(-ds3+ds1)        .*map(indexes.ldt(expr))...   %left  down  top    corner
            +(-ds4+ds1)        .*map(indexes.lub(expr))...   %left  upper bottom corner
            +(-ds1)            .*map(indexes.lut(expr))...   %left  upper top    corner
            +(sx-sxsy+ds2)     .*map(indexes.rdb(expr))...   %right down  bottom corner
            +(-ds2)            .*map(indexes.rdt(expr))...   %right down  top    corner
            +(sxsy-sxsysz)     .*map(indexes.rub(expr))...   %right upper bottom corner
            +(sxsysz)          .*map(indexes.rut(expr));     %right upper top    corner
end