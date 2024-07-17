function E = entropy_Vectorial_Data(data)
% data= column = variate rows = timestamp
[Alphabet,IA,IC] = unique(data,'rows');
Frequency = zeros(size(Alphabet,1),1);
%for i =1:size(data,1)
    for symbol = 1:size(Alphabet,1)
        pippo= repmat(Alphabet(symbol,:),size(data,1),1);
        Dist= data-pippo;
        if(size(Alphabet,2)==1)
            Frequency(symbol)=sum(sum(Dist'==0));
        else
            Frequency(symbol)=sum(sum(Dist')==0);
        end
        
    end
%end
divisor=sum(Frequency);
P = Frequency / divisor;
E = -sum(P .* log2(P)); 
end