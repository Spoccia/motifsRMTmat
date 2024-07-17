function overlap = overlapQueryFeature(depdin1,depdin2)
depdin1 = depdin1(depdin1~=0);
depdin2 = depdin2(depdin2~=0);
ov = intersect(depdin1(1,:),depdin2(1,:));
overlap = size(ov,2)/size(depdin1,2);
end