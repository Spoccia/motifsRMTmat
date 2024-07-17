function overlap = overlapJaccardSimilarity(depdin1,depdin2)
depdin1 = depdin1(depdin1~=0);
depdin2 = depdin2(depdin2~=0);
intersecFeatures = intersect(depdin1(:,1),depdin2(:,1));
unionFeatures = union(depdin1(:,1),depdin2(:,1));
overlap=0;
if(size(intersecFeatures, 2)==0)
    
else
    overlap = size(intersecFeatures, 1)/size(unionFeatures, 1);
end
end
