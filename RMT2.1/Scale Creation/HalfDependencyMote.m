function [ HalfDependencyMatrix, LocM, IDM ] = HalfDependencyMote(LocM, I)
%% kmeans
% Smooth half on the dependency dimension. 
% I: Input matrix, LocM: distance information between angles. 
% Num = floor(size(LocM,1)/2);
Num = floor(size(LocM,1)/2);
% HalfDependencyMatrix = zeros(size(I,1),size(I,2)/2);
HalfDependencyMatrix = zeros(size(I,1),floor(size(I,2)/2));
% HalfDependencyMatrix = zeros(size(I,1),floor(size(I,2)/2));
M = findinit(LocM, Num);
% fprintf('\n size of clustered matrix: %d, %d \n', size(M,1), size(M, 2));
[IDM, LocM] = kmeans(LocM, Num,'start', M,'emptyaction','singleton');
for i=1:Num
    t = find(IDM==i);
    for j=1:length(t)
        HalfDependencyMatrix(:,i) = HalfDependencyMatrix(:,i)+I(:,t(j));
    end
    HalfDependencyMatrix(:,i) = HalfDependencyMatrix(:,i)/length(t);
end

function M = findinit(A, k)
% A : Location Matirx
% k : number of clusters

s = size(A,1);
step = floor(s/k);
for i=1:k
    M(i,:) = A((i-1)*step+1,:);
end


