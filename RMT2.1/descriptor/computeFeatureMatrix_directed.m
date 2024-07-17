function fgss = computeFeatureMatrix_directed(octave,smin, oframes,H)
[M N S] = size(octave);

%for j=si:si+1
% construct matrix on the middle
for w=1:size(oframes,2)
    fgss{w} = zeros(M,N,S);
end

minv = min(oframes(3,:))+1;
maxv = max(oframes(3,:))+3;

for j=minv:maxv
    clear temp;
    for w = 1:size(oframes,2)
        xi = oframes(1,w);
        yi = oframes(2,w);
        si = oframes(3,w)-smin+1;
        fgss{w}(:,xi,j) = octave(:,xi,j);
    end
    % construct matrix on right side
    H1 = H;
    for i=1:N
        temp = (H1*octave(:,:,j)')';
        for w = 1:size(oframes,2)
            xi = oframes(1,w);
            yi = oframes(2,w);
            si = oframes(3,w)-smin+1;
            if ((xi+i)<=N)
                fgss{w}(:,xi+i,j) = temp(:,xi);
            end
        end
        H1 = H1*H;
    end
    % construct matrix on left side
    H1 = H';
    for i=1:N
        temp = (H1*octave(:,:,j)')';
        for w = 1:size(oframes,2)
            xi = oframes(1,w);
            yi = oframes(2,w);
            si = oframes(3,w)-smin+1;
            if ((xi-i)>0)
                fgss{w}(:,xi-i,j) = temp(:,xi);
            end
        end
        H1 = H1*H';
    end
end


