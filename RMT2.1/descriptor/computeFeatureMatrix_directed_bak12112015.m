function [fgss] = computeFeatureMatrix_directed_bak12112015(octave, oframes,H)
% each cell is a matrix
[M, N, Three, Four] = size(octave);
% octave = reshape(octave,M,N,Three*Four);
% S = Three*Four;

% construct matrix on the middle
fgss = cell(1,size(oframes,2));
for w=1:size(fgss,2) % time index
    % fgss{w} = zeros(M,N,S);
    fgss{w} = zeros(M,N,Three);
end
% minv = 1;
% maxv = S;



minv = min(oframes(3,:))+1;
maxv = max(oframes(3,:))+3;
for j=minv:maxv
    clear temp;
    
    %% initialize cell w
    for w = 1:size(oframes,2) % feature index
        xi = oframes(1,w); % sensor index
        fgss{w}(:,xi,j) = octave(:,xi,j,j);
    end
    
    %% construct matrix on right side
    H1 = H;
    for i=1:N % sensor index
        temp = (H1*octave(:,:,j,j)')';
        for w = 1:size(oframes,2) % feature index
            xi = oframes(1,w);
            if ((xi+i)<=N)
                fgss{w}(:,xi+i,j) = temp(:,xi);
            end
        end
        H1 = H1*H;
    end
    
    %% construct matrix on left side
    H1 = H';
    for i=1:N
        temp = (H1*octave(:,:,j,j)')';
        for w = 1:size(oframes,2)
            xi = oframes(1,w);
            
            if ((xi-i)>0)
                fgss{w}(:,xi-i,j) = temp(:,xi);
            end
        end
        H1 = H1*H';
    end
end


