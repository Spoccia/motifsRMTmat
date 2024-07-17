function [ idx ] = SiftLocalMax_Directed( input, N, M, Stime, Sdepd)
%% compute local extrema
% Stime --- row
% Sdepd --- column
% N  -- outgoing neighbours
% M -- incoming neighbours

% nneighbours = 3^size(size(input),1)-1; % in this case 26
% Skip the first and the last

% dss.octave{OTime, ODepd}{1} = RowVector; % depd difference
% dss.octave{OTime, ODepd}{2} = ColumnVector; % time difference
% dss.octave{OTime, ODepd}{3} = DiaVectorDown; % both directions

% input{2} -- time diff
% input{1} -- depd diff
% input{3} -- both diff
ttScale = size(input{2}, 3);
tdScale = size(input{2}, 4);

dtScale = size(input{1}, 3);
ddScale = size(input{1}, 4);

btScale = size(input{3}, 3);
bdScale = size(input{3}, 4);

idx = [];
timeStep = size(input{2}, 1);
depdStep = size(input{2}, 2);
timeScale = Stime-1;
depdScale = Sdepd-1;
tmax = 1;
for i = 1 : timeScale
    for j = 1 : depdScale
        for x = 2 : timeStep - 1
            for y = 2 : depdStep - 1
                % key points index: (x, y, i, j)
                temp = [abs(input{1}(x, y, i, j)), abs(input{2}(x, y, i, j)), abs(input{3}(x, y, i, j))];
                % temp = [input{1}(x, y, i, j), input{2}(x, y, i, j), input{3}(x, y, i, j)];
                tmax = max(temp);
                % tmax = min(temp);
                
                % temp = [input{1}(x, y, i, j), input{2}(x, y, i, j), input{3}(x, y, i, j)];
                % tmax = max(temp);
                if(i == 1) % first row
                    if(j == 1) % upper-left corner
                        neighbour = zeros(1,51);
                        %% current scale - 2
                        neighbour(1) = input{2}((x-1), :, i, j)*M(:, y);
                        neighbour(2) = input{2}(x, :, i, j)*M(:, y);
                        neighbour(3) = input{2}((x+1), :, i, j)*M(:, y);
                        
                        neighbour(4) = input{2}((x-1), y, i, j);
                        neighbour(5) = input{2}((x+1), y, i, j);
                        
                        neighbour(6) = input{2}((x-1), :, i, j)*N(:,y);
                        neighbour(7) = input{2}(x, :, i, j)*N(:,y);
                        neighbour(8) = input{2}((x+1), :, i, j)*N(:,y);
                        %% next scale on time direction
                        if(i+1 > ttScale)
                        else
                            neighbour(9) = input{2}((x-1), :, i+1, j)*M(:, y);
                            neighbour(10) = input{2}(x, :, i+1, j)*M(:, y);
                            neighbour(11) = input{2}((x+1), :, i+1, j)*M(:, y);
                            
                            neighbour(12) = input{2}((x-1), y, i+1, j);
                            neighbour(13) = input{2}(x, y, i+1, j);
                            neighbour(14) = input{2}((x+1), y, i+1, j);
                            
                            neighbour(15) = input{2}((x-1), :, i+1, j)*N(:,y);
                            neighbour(16) = input{2}(x, :, i+1, j)*N(:,y);
                            neighbour(17) = input{2}((x+1), :, i+1, j)*N(:,y);
                        end
                        %% current scale - 1 on depd
                        neighbour(18) = input{1}((x-1), :, i, j)*M(:, y);
                        neighbour(19) = input{1}(x, :, i, j)*M(:, y);
                        neighbour(20) = input{1}((x+1), :, i, j)*M(:, y);
                        
                        neighbour(21) = input{1}((x-1), y, i, j);
                        neighbour(22) = input{1}((x+1), y, i, j);
                        
                        neighbour(23) = input{1}((x-1), :, i, j)*N(:,y);
                        neighbour(24) = input{1}(x, :, i, j)*N(:,y);
                        neighbour(25) = input{1}((x+1), :, i, j)*N(:,y);
                        %% next scale on depd direction
                        if(j+1 > ddScale)
                        else
                            neighbour(26) = input{1}((x-1), :, i, j+1)*M(:, y+1);
                            neighbour(27) = input{1}(x, :, i, j+1)*M(:, y+1);
                            neighbour(28) = input{1}((x+1), :, i, j+1)*M(:, y+1);
                            
                            neighbour(29) = input{1}((x-1), y, i, j+1);
                            neighbour(30) = input{1}(x, y, i, j+1);
                            neighbour(31) = input{1}((x+1), y, i, j+1);
                            
                            neighbour(32) = input{1}((x-1), :, i, j+1)*N(:,y+1);
                            neighbour(33) = input{1}(x, :, i, j+1)*N(:,y+1);
                            neighbour(34) = input{1}((x+1), :, i, j+1)*N(:,y+1);
                        end
                        %% current scale on both dimension
                        neighbour(35) = input{3}((x-1), :, i, j)*M(:, y);
                        neighbour(36) = input{3}(x, :, i, j)*M(:, y);
                        neighbour(37) = input{3}((x+1), :, i, j)*M(:, y);
                        
                        neighbour(38) = input{3}((x-1), y, i, j);
                        neighbour(39) = input{3}((x+1), y, i, j);
                        
                        neighbour(40) = input{3}((x-1), :, i, j)*N(:,y);
                        neighbour(41) = input{3}(x, :, i, j)*N(:,y);
                        neighbour(42) = input{3}((x+1), :, i, j)*N(:,y);
                        %% next scale on both direction
                        % use or
                        if(i+1 > btScale || j + 1 > bdScale)
                        else
                            neighbour(43) = input{3}((x-1), :, i+1, j+1)*M(:, y+1);
                            neighbour(44) = input{3}(x, :, i+1, j+1)*M(:, y+1);
                            neighbour(45) = input{3}((x+1), :, i+1, j+1)*M(:, y+1);
                            
                            neighbour(46) = input{3}((x-1), y, i+1, j+1);
                            neighbour(47) = input{3}(x, y, i+1, j+1);
                            neighbour(48) = input{3}((x+1), y, i+1, j+1);
                            
                            neighbour(49) = input{3}((x-1), :, i+1, j+1)*N(:,y+1);
                            neighbour(50) = input{3}(x, :, i+1, j+1)*N(:,y+1);
                            neighbour(51) = input{3}((x+1), :, i+1, j+1)*N(:,y+1);
                        end
                    elseif( j == depdScale) % upper-right corner
                        neighbour = zeros(1,34);
                        %% current scale on time direction
                        neighbour(1) = input{2}((x-1), :, i, j)*M(:, y);
                        neighbour(2) = input{2}(x, :, i, j)*M(:, y);
                        neighbour(3) = input{2}((x+1), :, i, j)*M(:, y);
                        
                        neighbour(4) = input{2}((x-1), y, i, j);
                        neighbour(5) = input{2}((x+1), y, i, j);
                        
                        neighbour(6) = input{2}((x-1), :, i, j)*N(:,y);
                        neighbour(7) = input{2}(x, :, i, j)*N(:,y);
                        neighbour(8) = input{2}((x+1), :, i, j)*N(:,y);
                        %% next scale on time direction
                        if(i+1 > ttScale)
                        else
                            neighbour(9) = input{2}((x-1), :, i+1, j)*M(:, y);
                            neighbour(10) = input{2}(x, :, i+1, j)*M(:, y);
                            neighbour(11) = input{2}((x+1), :, i+1, j)*M(:, y);
                            
                            neighbour(12) = input{2}((x-1), y, i+1, j);
                            neighbour(13) = input{2}(x, y, i+1, j);
                            neighbour(14) = input{2}((x+1), y, i+1, j);
                            
                            neighbour(15) = input{2}((x-1), :, i+1, j)*N(:,y);
                            neighbour(16) = input{2}(x, :, i+1, j)*N(:,y);
                            neighbour(17) = input{2}((x+1), :, i+1, j)*N(:,y);
                        end
                        %% current scale on depd direction
                        neighbour(18) = input{1}((x-1), :, i, j)*M(:, y);
                        neighbour(19) = input{1}(x, :, i, j)*M(:, y);
                        neighbour(20) = input{1}((x+1), :, i, j)*M(:, y);
                        
                        neighbour(21) = input{1}((x-1), y, i, j);
                        neighbour(22) = input{1}((x+1), y, i, j);
                        
                        neighbour(23) = input{1}((x-1), :, i, j)*N(:,y);
                        neighbour(24) = input{1}(x, :, i, j)*N(:,y);
                        neighbour(25) = input{1}((x+1), :, i, j)*N(:,y);
                        %% previous scale on depd direction
                        if(j-1<=0)
                        else
                            neighbour(26) = input{1}((x-1), :, i, j-1)*M(:, j-1);
                            neighbour(27) = input{1}(x, :, i, j-1)*M(:, j-1);
                            neighbour(28) = input{1}((x+1), :, i, j-1)*M(:, j-1);
                            
                            neighbour(29) = input{1}((x-1), y, i, j-1);
                            neighbour(30) = input{1}(x, y, i, j-1);
                            neighbour(31) = input{1}((x+1), y, i, j-1);
                            
                            neighbour(32) = input{1}((x-1), :, i, j-1)*N(:,j-1);
                            neighbour(33) = input{1}(x, :, i, j-1)*N(:,j-1);
                            neighbour(34) = input{1}((x+1), :, i, j-1)*N(:,j-1);
                        end
                    else
                        neighbour = zeros(1,60);
                        %% current scale on time direction
                        neighbour(1) = input{2}((x-1), :, i, j)*M(:, y);
                        neighbour(2) = input{2}(x, :, i, j)*M(:, y);
                        neighbour(3) = input{2}((x+1), :, i, j)*M(:, y);
                        
                        neighbour(4) = input{2}((x-1), y, i, j);
                        neighbour(5) = input{2}((x+1), y, i, j);
                        
                        neighbour(6) = input{2}((x-1), :, i, j)*N(:,y);
                        neighbour(7) = input{2}(x, :, i, j)*N(:,y);
                        neighbour(8) = input{2}((x+1), :, i, j)*N(:,y);
                        %% next scale on time direction
                        if(i+1 > ttScale)
                        else
                            neighbour(9) = input{2}((x-1), :, i+1, j)*M(:, y);
                            neighbour(10) = input{2}(x, :, i+1, j)*M(:, y);
                            neighbour(11) = input{2}((x+1), :, i+1, j)*M(:, y);
                            
                            neighbour(12) = input{2}((x-1), y, i+1, j);
                            neighbour(13) = input{2}(x, y, i+1, j);
                            neighbour(14) = input{2}((x+1), y, i+1, j);
                            
                            neighbour(15) = input{2}((x-1), :, i+1, j)*N(:,y);
                            neighbour(16) = input{2}(x, :, i+1, j)*N(:,y);
                            neighbour(17) = input{2}((x+1), :, i+1, j)*N(:,y);
                        end
                        %% current scale on depd direction
                        neighbour(18) = input{1}((x-1), :, i, j)*M(:, y);
                        neighbour(19) = input{1}(x, :, i, j)*M(:, y);
                        neighbour(20) = input{1}((x+1), :, i, j)*M(:, y);
                        
                        neighbour(21) = input{1}((x-1), y, i, j);
                        neighbour(22) = input{1}((x+1), y, i, j);
                        
                        neighbour(23) = input{1}((x-1), :, i, j)*N(:,y);
                        neighbour(24) = input{1}(x, :, i, j)*N(:,y);
                        neighbour(25) = input{1}((x+1), :, i, j)*N(:,y);
                        %% next scale on depd direction
                        if(j+1 > ddScale)
                        else
                            neighbour(26) = input{1}((x-1), :, i, j+1)*M(:, y+1);
                            neighbour(27) = input{1}(x, :, i, j+1)*M(:, y+1);
                            neighbour(28) = input{1}((x+1), :, i, j+1)*M(:, y+1);
                            
                            neighbour(29) = input{1}((x-1), y, i, j+1);
                            neighbour(30) = input{1}(x, y, i, j+1);
                            neighbour(31) = input{1}((x+1), y, i, j+1);
                            
                            neighbour(32) = input{1}((x-1), :, i, j+1)*N(:,y+1);
                            neighbour(33) = input{1}(x, :, i, j+1)*N(:,y+1);
                            neighbour(34) = input{1}((x+1), :, i, j+1)*N(:,y+1);
                        end
                        %% current scale on both direction
                        neighbour(35) = input{3}((x-1), :, i, j)*M(:, y);
                        neighbour(36) = input{3}(x, :, i, j)*M(:, y);
                        neighbour(37) = input{3}((x+1), :, i, j)*M(:, y);
                        
                        neighbour(38) = input{3}((x-1), y, i, j);
                        neighbour(39) = input{3}((x+1), y, i, j);
                        
                        neighbour(40) = input{3}((x-1), :, i, j)*N(:,y);
                        neighbour(41) = input{3}(x, :, i, j)*N(:,y);
                        neighbour(42) = input{3}((x+1), :, i, j)*N(:,y);
                        %% next scale on both direction
                        if(i+1 > btScale || j+1> bdScale)
                        else
                            neighbour(43) = input{3}((x-1), :, i+1, j+1)*M(:, y+1);
                            neighbour(44) = input{3}(x, :, i+1, j+1)*M(:, y+1);
                            neighbour(45) = input{3}((x+1), :, i+1, j+1)*M(:, y+1);
                            
                            neighbour(46) = input{3}((x-1), y, i+1, j+1);
                            neighbour(47) = input{3}(x, y, i+1, j+1);
                            neighbour(48) = input{3}((x+1), y, i+1, j+1);
                            
                            neighbour(49) = input{3}((x-1), :, i+1, j+1)*N(:,y+1);
                            neighbour(50) = input{3}(x, :, i+1, j+1)*N(:,y+1);
                            neighbour(51) = input{3}((x+1), :, i+1, j+1)*N(:,y+1);
                        end
                        %% previous scale on depd direction
                        if(j-1 <=0)
                        else
                            neighbour(52) = input{1}((x-1), :, i, j-1)*M(:, j-1);
                            neighbour(53) = input{1}(x, :, i, j-1)*M(:, j-1);
                            neighbour(54) = input{1}((x+1), :, i, j-1)*M(:, j-1);
                            
                            neighbour(55) = input{1}((x-1), y, i, j-1);
                            neighbour(56) = input{1}(x, y, i, j-1);
                            neighbour(57) = input{1}((x+1), y, i, j-1);
                            
                            neighbour(58) = input{1}((x-1), :, i, j-1)*N(:,j-1);
                            neighbour(59) = input{1}(x, :, i, j-1)*N(:,j-1);
                            neighbour(60) = input{1}((x+1), :, i, j-1)*N(:,j-1);
                        end
                    end
                elseif(j == 1) % first column
                    if(i == 1) % upper-left corner
                    elseif( i == timeScale) % lower-left corner
                        neighbour = zeros(1, 34);
                        %% current scale on time direction
                        neighbour(1) = input{2}((x-1), :, i, j)*M(:, y);
                        neighbour(2) = input{2}(x, :, i, j)*M(:, y);
                        neighbour(3) = input{2}((x+1), :, i, j)*M(:, y);
                        
                        neighbour(4) = input{2}((x-1), y, i, j);
                        neighbour(5) = input{2}((x+1), y, i, j);
                        
                        neighbour(6) = input{2}((x-1), :, i, j)*N(:,y);
                        neighbour(7) = input{2}(x, :, i, j)*N(:,y);
                        neighbour(8) = input{2}((x+1), :, i, j)*N(:,y);
                        
                        %% previous scale on time direction
                        if(i-1 <=0)
                        else
                            neighbour(9) = input{2}((x-1), :, i-1, j)*M(:, y);
                            neighbour(10) = input{2}(x, :, i-1, j)*M(:, y);
                            neighbour(11) = input{2}((x+1), :, i-1, j)*M(:, y);
                            
                            neighbour(12) = input{2}((x-1), y, i-1, j);
                            neighbour(13) = input{2}(x, y, i-1, j);
                            neighbour(14) = input{2}((x+1), y, i-1, j);
                            
                            neighbour(15) = input{2}((x-1), :, i-1, j)*N(:,y);
                            neighbour(16) = input{2}(x, :, i-1, j)*N(:,y);
                            neighbour(17) = input{2}((x+1), :, i-1, j)*N(:,y);
                        end
                        %% current scale on depd direction
                        neighbour(18) = input{1}((x-1), :, i, j)*M(:, y);
                        neighbour(19) = input{1}(x, :, i, j)*M(:, y);
                        neighbour(20) = input{1}((x+1), :, i, j)*M(:, y);
                        
                        neighbour(21) = input{1}((x-1), y, i, j);
                        neighbour(22) = input{1}((x+1), y, i, j);
                        
                        neighbour(23) = input{1}((x-1), :, i, j)*N(:,y);
                        neighbour(24) = input{1}(x, :, i, j)*N(:,y);
                        neighbour(25) = input{1}((x+1), :, i, j)*N(:,y);
                        %% next scale on depd direction
                        if(j+1 > ddScale)
                        else
                            neighbour(26) = input{1}((x-1), :, i, j+1)*M(:, y+1);
                            neighbour(27) = input{1}(x, :, i, j+1)*M(:, y+1);
                            neighbour(28) = input{1}((x+1), :, i, j+1)*M(:, y+1);
                            
                            neighbour(29) = input{1}((x-1), y, i, j+1);
                            neighbour(30) = input{1}(x, y, i, j+1);
                            neighbour(31) = input{1}((x+1), y, i, j+1);
                            
                            neighbour(32) = input{1}((x-1), :, i, j+1)*N(:,y+1);
                            neighbour(33) = input{1}(x, :, i, j+1)*N(:,y+1);
                            neighbour(34) = input{1}((x+1), :, i, j+1)*N(:,y+1);
                        end
                    else
                        neighbour = zeros(1,60);
                        %% current scale on time direction
                        neighbour(1) = input{2}((x-1), :, i, j)*M(:, y);
                        neighbour(2) = input{2}(x, :, i, j)*M(:, y);
                        neighbour(3) = input{2}((x+1), :, i, j)*M(:, y);
                        
                        neighbour(4) = input{2}((x-1), y, i, j);
                        neighbour(5) = input{2}((x+1), y, i, j);
                        
                        neighbour(6) = input{2}((x-1), :, i, j)*N(:,y);
                        neighbour(7) = input{2}(x, :, i, j)*N(:,y);
                        neighbour(8) = input{2}((x+1), :, i, j)*N(:,y);
                        %% next scale on time direction
                        if(i+1 > ttScale)
                        else
                            neighbour(9) = input{2}((x-1), :, i+1, j)*M(:, y);
                            neighbour(10) = input{2}(x, :, i+1, j)*M(:, y);
                            neighbour(11) = input{2}((x+1), :, i+1, j)*M(:, y);
                            
                            neighbour(12) = input{2}((x-1), y, i+1, j);
                            neighbour(13) = input{2}(x, y, i+1, j);
                            neighbour(14) = input{2}((x+1), y, i+1, j);
                            
                            neighbour(15) = input{2}((x-1), :, i+1, j)*N(:,y);
                            neighbour(16) = input{2}(x, :, i+1, j)*N(:,y);
                            neighbour(17) = input{2}((x+1), :, i+1, j)*N(:,y);
                        end
                        %% previous scale on time direction
                        if(i-1 <=0)
                        else
                            neighbour(18) = input{2}((x-1), :, i-1, j)*M(:, y);
                            neighbour(19) = input{2}(x, :, i-1, j)*M(:, y);
                            neighbour(20) = input{2}((x+1), :, i-1, j)*M(:, y);
                            
                            neighbour(21) = input{2}((x-1), y, i-1, j);
                            neighbour(22) = input{2}(x, y, i-1, j);
                            neighbour(23) = input{2}((x+1), y, i-1, j);
                            
                            neighbour(24) = input{2}((x-1), :, i-1, j)*N(:,y);
                            neighbour(25) = input{2}(x, :, i-1, j)*N(:,y);
                            neighbour(26) = input{2}((x+1), :, i-1, j)*N(:,y);
                        end
                        
                        %% current scale on depd direction
                        neighbour(27) = input{1}((x-1), :, i, j)*M(:, y);
                        neighbour(28) = input{1}(x, :, i, j)*M(:, y);
                        neighbour(29) = input{1}((x+1), :, i, j)*M(:, y);
                        
                        neighbour(30) = input{1}((x-1), y, i, j);
                        neighbour(31) = input{1}((x+1), y, i, j);
                        
                        neighbour(32) = input{1}((x-1), :, i, j)*N(:,y);
                        neighbour(33) = input{1}(x, :, i, j)*N(:,y);
                        neighbour(34) = input{1}((x+1), :, i, j)*N(:,y);
                        %% next scale on depd direction
                        if(j+1 > ddScale)
                        else
                            neighbour(35) = input{1}((x-1), :, i, j+1)*M(:, y+1);
                            neighbour(36) = input{1}(x, :, i, j+1)*M(:, y+1);
                            neighbour(37) = input{1}((x+1), :, i, j+1)*M(:, y+1);
                            
                            neighbour(38) = input{1}((x-1), y, i, j+1);
                            neighbour(39) = input{1}(x, y, i, j+1);
                            neighbour(40) = input{1}((x+1), y, i, j+1);
                            
                            neighbour(41) = input{1}((x-1), :, i, j+1)*N(:,y+1);
                            neighbour(42) = input{1}(x, :, i, j+1)*N(:,y+1);
                            neighbour(43) = input{1}((x+1), :, i, j+1)*N(:,y+1);
                        end
                        %% current scale on both direction
                        neighbour(44) = input{3}((x-1), :, i, j)*M(:, y);
                        neighbour(45) = input{3}(x, :, i, j)*M(:, y);
                        neighbour(46) = input{3}((x+1), :, i, j)*M(:, y);
                        
                        neighbour(47) = input{3}((x-1), y, i, j);
                        neighbour(48) = input{3}((x+1), y, i, j);
                        
                        neighbour(49) = input{3}((x-1), :, i, j)*N(:,y);
                        neighbour(50) = input{3}(x, :, i, j)*N(:,y);
                        neighbour(51) = input{3}((x+1), :, i, j)*N(:,y);
                        %% next scale on both direction
                        if(i+1 >btScale || j+1>bdScale)
                        else
                            neighbour(52) = input{3}((x-1), :, i+1, j+1)*M(:, y+1);
                            neighbour(53) = input{3}(x, :, i+1, j+1)*M(:, y+1);
                            neighbour(54) = input{3}((x+1), :, i+1, j+1)*M(:, y+1);
                            
                            neighbour(55) = input{3}((x-1), y, i+1, j+1);
                            neighbour(56) = input{3}(x, y, i+1, j+1);
                            neighbour(57) = input{3}((x+1), y, i+1, j+1);
                            
                            neighbour(58) = input{3}((x-1), :, i+1, j+1)*N(:,y+1);
                            neighbour(59) = input{3}(x, :, i+1, j+1)*N(:,y+1);
                            neighbour(60) = input{3}((x+1), :, i+1, j+1)*N(:,y+1);
                        end
                    end
                elseif(i == timeScale) % last row
                    if( j == 1) % lower-left corner
                    elseif( j == depdScale) % lower-right corner
                        neighbour = zeros(1,51);
                        %% current scale on time direction
                        neighbour(1) = input{2}((x-1), :, i, j)*M(:, y);
                        neighbour(2) = input{2}(x, :, i, j)*M(:, y);
                        neighbour(3) = input{2}((x+1), :, i, j)*M(:, y);
                        
                        neighbour(4) = input{2}((x-1), y, i, j);
                        neighbour(5) = input{2}((x+1), y, i, j);
                        
                        neighbour(6) = input{2}((x-1), :, i, j)*N(:,y);
                        neighbour(7) = input{2}(x, :, i, j)*N(:,y);
                        neighbour(8) = input{2}((x+1), :, i, j)*N(:,y);
                        
                        %% previous scale on time direction
                        if(i-1<=0)
                        else
                            neighbour(9) = input{2}((x-1), :, i-1, j)*M(:, y);
                            neighbour(10) = input{2}(x, :, i-1, j)*M(:, y);
                            neighbour(11) = input{2}((x+1), :, i-1, j)*M(:, y);
                            
                            neighbour(12) = input{2}((x-1), y, i-1, j);
                            neighbour(13) = input{2}(x, y, i-1, j);
                            neighbour(14) = input{2}((x+1), y, i-1, j);
                            
                            neighbour(15) = input{2}((x-1), :, i-1, j)*N(:,y);
                            neighbour(16) = input{2}(x, :, i-1, j)*N(:,y);
                            neighbour(17) = input{2}((x+1), :, i-1, j)*N(:,y);
                        end
                        
                        %% current scale on depd direction
                        neighbour(18) = input{1}((x-1), :, i, j)*M(:, y);
                        neighbour(19) = input{1}(x, :, i, j)*M(:, y);
                        neighbour(20) = input{1}((x+1), :, i, j)*M(:, y);
                        
                        neighbour(21) = input{1}((x-1), y, i, j);
                        neighbour(22) = input{1}((x+1), y, i, j);
                        
                        neighbour(23) = input{1}((x-1), :, i, j)*N(:,y);
                        neighbour(24) = input{1}(x, :, i, j)*N(:,y);
                        neighbour(25) = input{1}((x+1), :, i, j)*N(:,y);
                        %% previous scale on depd direction
                        if(j-1<=0)
                        else
                            neighbour(26) = input{1}((x-1), :, i, j-1)*M(:, j-1);
                            neighbour(27) = input{1}(x, :, i, j-1)*M(:, j-1);
                            neighbour(28) = input{1}((x+1), :, i, j-1)*M(:, j-1);
                            
                            neighbour(29) = input{1}((x-1), y, i, j-1);
                            neighbour(30) = input{1}(x, y, i, j-1);
                            neighbour(31) = input{1}((x+1), y, i, j-1);
                            
                            neighbour(32) = input{1}((x-1), :, i, j-1)*N(:,j-1);
                            neighbour(33) = input{1}(x, :, i, j-1)*N(:,j-1);
                            neighbour(34) = input{1}((x+1), :, i, j-1)*N(:,j-1);
                        end
                        
                        %% current scale on both direction
                        neighbour(35) = input{3}((x-1), :, i, j)*M(:, y);
                        neighbour(36) = input{3}(x, :, i, j)*M(:, y);
                        neighbour(37) = input{3}((x+1), :, i, j)*M(:, y);
                        
                        neighbour(38) = input{3}((x-1), y, i, j);
                        neighbour(39) = input{3}((x+1), y, i, j);
                        
                        neighbour(40) = input{3}((x-1), :, i, j)*N(:,y);
                        neighbour(41) = input{3}(x, :, i, j)*N(:,y);
                        neighbour(42) = input{3}((x+1), :, i, j)*N(:,y);
                        
                        %% previous scale on both direction
                        if(i-1<=0||j-1<=0)
                        else
                            neighbour(43) = input{3}((x-1), :, i-1, j-1)*M(:, j-1);
                            neighbour(44) = input{3}(x, :, i-1, j-1)*M(:, j-1);
                            neighbour(45) = input{3}((x+1), :, i-1, j-1)*M(:, j-1);
                            
                            neighbour(46) = input{3}((x-1), y, i-1, j-1);
                            neighbour(47) = input{3}(x, y, i-1, j-1);
                            neighbour(48) = input{3}((x+1), y, i-1, j-1);
                            
                            neighbour(49) = input{3}((x-1), :, i-1, j-1)*N(:,j-1);
                            neighbour(50) = input{3}(x, :, i-1, j-1)*N(:,j-1);
                            neighbour(51) = input{3}((x+1), :, i-1, j-1)*N(:,j-1);
                        end
                    else
                        neighbour = zeros(1,60);
                        %% current scale on time direction
                        neighbour(1) = input{2}((x-1), :, i, j)*M(:, y);
                        neighbour(2) = input{2}(x, :, i, j)*M(:, y);
                        neighbour(3) = input{2}((x+1), :, i, j)*M(:, y);
                        
                        neighbour(4) = input{2}((x-1), y, i, j);
                        neighbour(5) = input{2}((x+1), y, i, j);
                        
                        neighbour(6) = input{2}((x-1), :, i, j)*N(:,y);
                        neighbour(7) = input{2}(x, :, i, j)*N(:,y);
                        neighbour(8) = input{2}((x+1), :, i, j)*N(:,y);
                        
                        %% previous scale on time direction
                        if(i-1<=0)
                        else
                            neighbour(9) = input{2}((x-1), :, i-1, j)*M(:, y);
                            neighbour(10) = input{2}(x, :, i-1, j)*M(:, y);
                            neighbour(11) = input{2}((x+1), :, i-1, j)*M(:, y);
                            
                            neighbour(12) = input{2}((x-1), y, i-1, j);
                            neighbour(13) = input{2}(x, y, i-1, j);
                            neighbour(14) = input{2}((x+1), y, i-1, j);
                            
                            neighbour(15) = input{2}((x-1), :, i-1, j)*N(:,y);
                            neighbour(16) = input{2}(x, :, i-1, j)*N(:,y);
                            neighbour(17) = input{2}((x+1), :, i-1, j)*N(:,y);
                        end
                        
                        %% current scale on depd direction
                        neighbour(18) = input{1}((x-1), :, i, j)*M(:, y);
                        neighbour(19) = input{1}(x, :, i, j)*M(:, y);
                        neighbour(20) = input{1}((x+1), :, i, j)*M(:, y);
                        
                        neighbour(21) = input{1}((x-1), y, i, j);
                        neighbour(22) = input{1}((x+1), y, i, j);
                        
                        neighbour(23) = input{1}((x-1), :, i, j)*N(:,y);
                        neighbour(24) = input{1}(x, :, i, j)*N(:,y);
                        neighbour(25) = input{1}((x+1), :, i, j)*N(:,y);
                        %% next scale on depd direction
                        if(j+1 > ddScale)
                        else
                            neighbour(26) = input{1}((x-1), :, i, j+1)*M(:, y+1);
                            neighbour(27) = input{1}(x, :, i, j+1)*M(:, y+1);
                            neighbour(28) = input{1}((x+1), :, i, j+1)*M(:, y+1);
                            
                            neighbour(29) = input{1}((x-1), y, i, j+1);
                            neighbour(30) = input{1}(x, y, i, j+1);
                            neighbour(31) = input{1}((x+1), y, i, j+1);
                            
                            neighbour(32) = input{1}((x-1), :, i, j+1)*N(:,y+1);
                            neighbour(33) = input{1}(x, :, i, j+1)*N(:,y+1);
                            neighbour(34) = input{1}((x+1), :, i, j+1)*N(:,y+1);
                        end
                        %% previous scale on depd direction
                        if(j-1 <=0)
                        else
                            neighbour(35) = input{1}((x-1), :, i, j-1)*M(:, j-1);
                            neighbour(36) = input{1}(x, :, i, j-1)*M(:, j-1);
                            neighbour(37) = input{1}((x+1), :, i, j-1)*M(:, j-1);
                            
                            neighbour(38) = input{1}((x-1), y, i, j-1);
                            neighbour(39) = input{1}(x, y, i, j-1);
                            neighbour(40) = input{1}((x+1), y, i, j-1);
                            
                            neighbour(41) = input{1}((x-1), :, i, j-1)*N(:,j-1);
                            neighbour(42) = input{1}(x, :, i, j-1)*N(:,j-1);
                            neighbour(43) = input{1}((x+1), :, i, j-1)*N(:,j-1);
                        end
                        
                        %% current scale on both direction
                        neighbour(44) = input{3}((x-1), :, i, j)*M(:, y);
                        neighbour(45) = input{3}(x, :, i, j)*M(:, y);
                        neighbour(46) = input{3}((x+1), :, i, j)*M(:, y);
                        
                        neighbour(47) = input{3}((x-1), y, i, j);
                        neighbour(48) = input{3}((x+1), y, i, j);
                        
                        neighbour(49) = input{3}((x-1), :, i, j)*N(:,y);
                        neighbour(50) = input{3}(x, :, i, j)*N(:,y);
                        neighbour(51) = input{3}((x+1), :, i, j)*N(:,y);
                        
                        %% previous scale on both direction
                        if(i-1 <=0 || j-1<=0)
                        else
                            neighbour(52) = input{3}((x-1), :, i-1, j-1)*M(:, j-1);
                            neighbour(53) = input{3}(x, :, i-1, j-1)*M(:, j-1);
                            neighbour(54) = input{3}((x+1), :, i-1, j-1)*M(:, j-1);
                            
                            neighbour(55) = input{3}((x-1), y, i-1, j-1);
                            neighbour(56) = input{3}(x, y, i-1, j-1);
                            neighbour(57) = input{3}((x+1), y, i-1, j-1);
                            
                            neighbour(58) = input{3}((x-1), :, i-1, j-1)*N(:,j-1);
                            neighbour(59) = input{3}(x, :, i-1, j-1)*N(:,j-1);
                            neighbour(60) = input{3}((x+1), :, i-1, j-1)*N(:,j-1);
                        end
                    end
                elseif(j == depdScale) % last column
                    if(i == 1) % upper-right corner
                    elseif(i == timeScale) % lower-right corner
                    else
                        neighbour = zeros(1,60);
                        %% current scale on time direction
                        neighbour(1) = input{2}((x-1), :, i, j)*M(:, y);
                        neighbour(2) = input{2}(x, :, i, j)*M(:, y);
                        neighbour(3) = input{2}((x+1), :, i, j)*M(:, y);
                        
                        neighbour(4) = input{2}((x-1), y, i, j);
                        neighbour(5) = input{2}((x+1), y, i, j);
                        
                        neighbour(6) = input{2}((x-1), :, i, j)*N(:,y);
                        neighbour(7) = input{2}(x, :, i, j)*N(:,y);
                        neighbour(8) = input{2}((x+1), :, i, j)*N(:,y);
                        %% previous scale on time direction
                        if(i-1<=0)
                        else
                            neighbour(9) = input{2}((x-1), :, i-1, j)*M(:, y);
                            neighbour(10) = input{2}(x, :, i-1, j)*M(:, y);
                            neighbour(11) = input{2}((x+1), :, i-1, j)*M(:, y);
                            
                            neighbour(12) = input{2}((x-1), y, i-1, j);
                            neighbour(13) = input{2}(x, y, i-1, j);
                            neighbour(14) = input{2}((x+1), y, i-1, j);
                            
                            neighbour(15) = input{2}((x-1), :, i-1, j)*N(:,y);
                            neighbour(16) = input{2}(x, :, i-1, j)*N(:,y);
                            neighbour(17) = input{2}((x+1), :, i-1, j)*N(:,y);
                        end
                        %% next scale on time direction
                        if(i + 1>ttScale)
                        else
                            neighbour(18) = input{2}((x-1), :, i+1, j)*M(:, y);
                            neighbour(19) = input{2}(x, :, i+1, j)*M(:, y);
                            neighbour(20) = input{2}((x+1), :, i+1, j)*M(:, y);
                            
                            neighbour(21) = input{2}((x-1), y, i+1, j);
                            neighbour(22) = input{2}(x, y, i+1, j);
                            neighbour(23) = input{2}((x+1), y, i+1, j);
                            
                            neighbour(24) = input{2}((x-1), :, i+1, j)*N(:,y);
                            neighbour(25) = input{2}(x, :, i+1, j)*N(:,y);
                            neighbour(26) = input{2}((x+1), :, i+1, j)*N(:,y);
                        end
                        
                        %% current scale on depd direction
                        neighbour(27) = input{1}((x-1), :, i, j)*M(:, y);
                        neighbour(28) = input{1}(x, :, i, j)*M(:, y);
                        neighbour(29) = input{1}((x+1), :, i, j)*M(:, y);
                        
                        neighbour(30) = input{1}((x-1), y, i, j);
                        neighbour(31) = input{1}((x+1), y, i, j);
                        
                        neighbour(32) = input{1}((x-1), :, i, j)*N(:,y);
                        neighbour(33) = input{1}(x, :, i, j)*N(:,y);
                        neighbour(34) = input{1}((x+1), :, i, j)*N(:,y);
                        %% previous scale on depd direction
                        if(j-1<=0)
                        else
                            neighbour(35) = input{1}((x-1), :, i, j-1)*M(:, j-1);
                            neighbour(36) = input{1}(x, :, i, j-1)*M(:, j-1);
                            neighbour(37) = input{1}((x+1), :, i, j-1)*M(:, j-1);
                            
                            neighbour(38) = input{1}((x-1), y, i, j-1);
                            neighbour(39) = input{1}(x, y, i, j-1);
                            neighbour(40) = input{1}((x+1), y, i, j-1);
                            
                            neighbour(41) = input{1}((x-1), :, i, j-1)*N(:,j-1);
                            neighbour(42) = input{1}(x, :, i, j-1)*N(:,j-1);
                            neighbour(43) = input{1}((x+1), :, i, j-1)*N(:,j-1);
                        end
                        %% current scale on both direction
                        neighbour(44) = input{3}((x-1), :, i, j)*M(:, y);
                        neighbour(45) = input{3}(x, :, i, j)*M(:, y);
                        neighbour(46) = input{3}((x+1), :, i, j)*M(:, y);
                        
                        neighbour(47) = input{3}((x-1), y, i, j);
                        neighbour(48) = input{3}((x+1), y, i, j);
                        
                        neighbour(49) = input{3}((x-1), :, i, j)*N(:,y);
                        neighbour(50) = input{3}(x, :, i, j)*N(:,y);
                        neighbour(51) = input{3}((x+1), :, i, j)*N(:,y);
                        
                        %% previous scale on both direction
                        if(i-1<=0||j-1<=0)
                        else
                            neighbour(52) = input{3}((x-1), :, i-1, j-1)*M(:, j-1);
                            neighbour(53) = input{3}(x, :, i-1, j-1)*M(:, j-1);
                            neighbour(54) = input{3}((x+1), :, i-1, j-1)*M(:, j-1);
                            
                            neighbour(55) = input{3}((x-1), y, i-1, j-1);
                            neighbour(56) = input{3}(x, y, i-1, j-1);
                            neighbour(57) = input{3}((x+1), y, i-1, j-1);
                            
                            neighbour(58) = input{3}((x-1), :, i-1, j-1)*N(:,j-1);
                            neighbour(59) = input{3}(x, :, i-1, j-1)*N(:,j-1);
                            neighbour(60) = input{3}((x+1), :, i-1, j-1)*N(:,j-1);
                        end
                    end
                else % internal condition
                    neighbour = zeros(1,78);
                    %% current scale on time direction
                    neighbour(1) = input{2}((x-1), :, i, j)*M(:, y);
                    neighbour(2) = input{2}(x, :, i, j)*M(:, y);
                    neighbour(3) = input{2}((x+1), :, i, j)*M(:, y);
                    
                    neighbour(4) = input{2}((x-1), y, i, j);
                    neighbour(5) = input{2}((x+1), y, i, j);
                    
                    neighbour(6) = input{2}((x-1), :, i, j)*N(:,y);
                    neighbour(7) = input{2}(x, :, i, j)*N(:,y);
                    neighbour(8) = input{2}((x+1), :, i, j)*N(:,y);
                    %% previous scale on time direction
                    if(i-1<=0)
                    else
                        neighbour(9) = input{2}((x-1), :, i-1, j)*M(:, y);
                        neighbour(10) = input{2}(x, :, i-1, j)*M(:, y);
                        neighbour(11) = input{2}((x+1), :, i-1, j)*M(:, y);
                        
                        neighbour(12) = input{2}((x-1), y, i-1, j);
                        neighbour(13) = input{2}(x, y, i-1, j);
                        neighbour(14) = input{2}((x+1), y, i-1, j);
                        
                        neighbour(15) = input{2}((x-1), :, i-1, j)*N(:,y);
                        neighbour(16) = input{2}(x, :, i-1, j)*N(:,y);
                        neighbour(17) = input{2}((x+1), :, i-1, j)*N(:,y);
                    end
                    %% next scale on time direction
                    if(i+1 > ttScale)
                    else
                        neighbour(18) = input{2}((x-1), :, i+1, j)*M(:, y);
                        neighbour(19) = input{2}(x, :, i+1, j)*M(:, y);
                        neighbour(20) = input{2}((x+1), :, i+1, j)*M(:, y);
                        
                        neighbour(21) = input{2}((x-1), y, i+1, j);
                        neighbour(22) = input{2}(x, y, i+1, j);
                        neighbour(23) = input{2}((x+1), y, i+1, j);
                        
                        neighbour(24) = input{2}((x-1), :, i+1, j)*N(:,y);
                        neighbour(25) = input{2}(x, :, i+1, j)*N(:,y);
                        neighbour(26) = input{2}((x+1), :, i+1, j)*N(:,y);
                    end
                    
                    %% current scale on depd direction
                    neighbour(27) = input{1}((x-1), :, i, j)*M(:, y);
                    neighbour(28) = input{1}(x, :, i, j)*M(:, y);
                    neighbour(29) = input{1}((x+1), :, i, j)*M(:, y);
                    
                    neighbour(30) = input{1}((x-1), y, i, j);
                    neighbour(31) = input{1}((x+1), y, i, j);
                    
                    neighbour(32) = input{1}((x-1), :, i, j)*N(:,y);
                    neighbour(33) = input{1}(x, :, i, j)*N(:,y);
                    neighbour(34) = input{1}((x+1), :, i, j)*N(:,y);
                    %% previous scale on depd direction
                    if(j-1 <=0)
                    else
                        neighbour(35) = input{1}((x-1), :, i, j-1)*M(:, j-1);
                        neighbour(36) = input{1}(x, :, i, j-1)*M(:, j-1);
                        neighbour(37) = input{1}((x+1), :, i, j-1)*M(:, j-1);
                        
                        neighbour(38) = input{1}((x-1), y, i, j-1);
                        neighbour(39) = input{1}(x, y, i, j-1);
                        neighbour(40) = input{1}((x+1), y, i, j-1);
                        
                        neighbour(41) = input{1}((x-1), :, i, j-1)*N(:,j-1);
                        neighbour(42) = input{1}(x, :, i, j-1)*N(:,j-1);
                        neighbour(43) = input{1}((x+1), :, i, j-1)*N(:,j-1);
                    end
                    %% next scale on depd direction
                    if(j+1 > ddScale)
                    else
                        neighbour(44) = input{1}((x-1), :, i, j+1)*M(:, y+1);
                        neighbour(45) = input{1}(x, :, i, j+1)*M(:, y+1);
                        neighbour(46) = input{1}((x+1), :, i, j+1)*M(:, y+1);
                        
                        neighbour(47) = input{1}((x-1), y, i, j+1);
                        neighbour(48) = input{1}(x, y, i, j+1);
                        neighbour(49) = input{1}((x+1), y, i, j+1);
                        
                        neighbour(50) = input{1}((x-1), :, i, j+1)*N(:,y+1);
                        neighbour(51) = input{1}(x, :, i, j+1)*N(:,y+1);
                        neighbour(52) = input{1}((x+1), :, i, j+1)*N(:,y+1);
                    end
                    %% current scale on both direction
                    neighbour(53) = input{3}((x-1), :, i, j)*M(:, y);
                    neighbour(54) = input{3}(x, :, i, j)*M(:, y);
                    neighbour(55) = input{3}((x+1), :, i, j)*M(:, y);
                    
                    neighbour(56) = input{3}((x-1), y, i, j);
                    neighbour(57) = input{3}((x+1), y, i, j);
                    
                    neighbour(58) = input{3}((x-1), :, i, j)*N(:,y);
                    neighbour(59) = input{3}(x, :, i, j)*N(:,y);
                    neighbour(60) = input{3}((x+1), :, i, j)*N(:,y);
                    
                    %% previous scale on both direction
                    if(i-1<=0||j-1<=0)
                    else
                        neighbour(61) = input{3}((x-1), :, i-1, j-1)*M(:, j-1);
                        neighbour(62) = input{3}(x, :, i-1, j-1)*M(:, j-1);
                        neighbour(63) = input{3}((x+1), :, i-1, j-1)*M(:, j-1);
                        
                        neighbour(64) = input{3}((x-1), y, i-1, j-1);
                        neighbour(65) = input{3}(x, y, i-1, j-1);
                        neighbour(66) = input{3}((x+1), y, i-1, j-1);
                        
                        neighbour(67) = input{3}((x-1), :, i-1, j-1)*N(:,j-1);
                        neighbour(68) = input{3}(x, :, i-1, j-1)*N(:,j-1);
                        neighbour(69) = input{3}((x+1), :, i-1, j-1)*N(:,j-1);
                    end
                    
                    %% next scale on both direction
                    if(i+1 > btScale || j+1 > bdScale)
                    else
                        neighbour(70) = input{3}((x-1), :, i+1, j+1)*M(:, y+1);
                        neighbour(71) = input{3}(x, :, i+1, j+1)*M(:, y+1);
                        neighbour(72) = input{3}((x+1), :, i+1, j+1)*M(:, y+1);
                        
                        neighbour(73) = input{3}((x-1), y, i+1, j+1);
                        neighbour(74) = input{3}(x, y, i+1, j+1);
                        neighbour(75) = input{3}((x+1), y, i+1, j+1);
                        
                        neighbour(76) = input{3}((x-1), :, i+1, j+1)*N(:,y+1);
                        neighbour(77) = input{3}(x, :, i+1, j+1)*N(:,y+1);
                        neighbour(78) = input{3}((x+1), :, i+1, j+1)*N(:,y+1);
                    end
                end
                if(tmax>0.9*abs(neighbour(:)))% this is a key point
                % if(tmax>abs(neighbour(:)))% this is a key point
                % if(tmax>neighbour(:))% this is a key point
                    idx = cat(1, idx, [x, y, i, j, tmax]);
                end
                clear neighbour
            end
        end
    end
end
idx = idx';
end
