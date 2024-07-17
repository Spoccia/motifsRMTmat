function [ H ] = NormalizeByRow( H )
%NORMALIZEH Summary of this function goes here
%   Detailed explanation goes here

for i = 1:size(H,1)
    if( sum(H(i,:))>0 )
        H(i, :) = H(i, :) / sum(H(i, :));
    end
end

end

