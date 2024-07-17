function [ HH ] = NormalizeB( H )

HH = H;
% HH = tril(H);
for i = 1:size(HH,1)
    if( sum(HH(i,:))>0 )
        HH(i, :) = HH(i, :) / sum(HH(i, :));
    end
end

end


