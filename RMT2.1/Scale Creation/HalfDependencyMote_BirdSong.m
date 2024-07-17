function [ HalfDependencyMatrix, IDM ] = HalfDependencyMote_BirdSong(LocM, I, IDM)
%%% EMG Data No Kmeans Clustering Used
% Instead use predefined parameters
% Smooth half on the dependency dimension. 
% I: Input matrix, LocM: distance information between angles. 

% Num = floor(size(LocM,1)/2);
Num = max(IDM(:));

HalfDependencyMatrix = zeros(size(I,1), Num);

% fprintf('\n size of clustered matrix: %d, %d \n', size(M,1), size(M, 2));
for i=1:Num
    t = find(IDM==i);
    for j=1:length(t)
        HalfDependencyMatrix(:,i) = HalfDependencyMatrix(:,i)+I(:,t(j));
    end
    HalfDependencyMatrix(:,i) = HalfDependencyMatrix(:,i)/length(t);
end