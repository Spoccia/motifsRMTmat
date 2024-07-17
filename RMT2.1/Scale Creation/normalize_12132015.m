function H = normalize_12132015(H)
Ht = H';
for i=1:size(H,1) % = N
    if(sum(H(i,:))>0)
        H(i,:) = H(i,:)/sum(H(i,:));
    end
end
end