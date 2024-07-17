function quantizedMatrix = singlevariateQuantization(Matrice)
quantizedMatrix = zeros(size(Matrice));
maxvalue=15;%64;%128;%255;
m=min(Matrice(:));
M=max(Matrice(:));
quantizedMatrix =round((Matrice - m) * (maxvalue / (M-m)))+1;
% for i= 1: size(Matrice,2)
%     quantizedMatrix(:,i)=round((Matrice(:,i) - m(i)) * (255 / (M(i)-m(i))));
% end
end