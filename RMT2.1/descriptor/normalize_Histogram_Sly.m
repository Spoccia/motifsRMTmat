function sd = normalize_Histogram_Sly(descr)
    norm=0 ;
    for i = 1: size(descr,1)
        norm = norm + (descr(i,1)) * (descr(i,1)) ;
    end
    norm = sqrt(norm) ;
    sd = descr/norm;
%    v = descr/norm(descr);
end