function quantizedMatrix = singlevariateQuantization(Matrice)
quantizedMatrix = zeros(size(Matrice));

m=min(Matrice);
M=max(Matrice);
maximum = 15;%255;
for i= 1: size(Matrice,2)
    quantizedMatrix(:,i)=round((Matrice(:,i) - m(i)) * (maximum / (M(i)-m(i))))+1;
end
end