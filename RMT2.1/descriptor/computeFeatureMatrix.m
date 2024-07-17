function fgss = computeFeatureMatrix(octave,smin, oframes,H)
[M N S] = size(octave);
xi = max(floor(oframes(1,1)+0.5),1);
xi = min(xi,N);
yi = max(floor(oframes(2,1)+0.5),1);
yi = min(yi,M);
si = min(S,floor(oframes(3,1)+0.5) - smin);
si = max(1,si);
%for j=si:si+1
% construct matrix on the middle
fgss = zeros(M,N,S);
minv = max(1,si-1);
maxv = min(si+1,S);
for j=minv:maxv
fgss(:,xi,j) = octave(:,xi,si);
% construct matrix on right side
for i=xi+1:N
    temp = (H^(i-xi)*octave(:,:,si)')';
    fgss(:,i,j) = temp(:,xi);
end
% construct matrix on left side
for i=1:xi-1
    temp = (pinv(H)^(xi-i)*octave(:,:,si)')';
    fgss(:,i,j) = temp(:,xi);
end
end
