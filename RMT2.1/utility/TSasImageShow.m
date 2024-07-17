function [ output_args ] = TSasImageShow( TS )
%TSASIMAGESHOW Summary of this function goes here
%   Detailed explanation goes here

    manip_f = TS - min(TS(:));
    %% normalize 0-1
    manip_f = manip_f/(max(TS(:))-min(TS(:)));
    %normalize 0-255
    manip2_f = manip_f*255;
    %%Increase the number of variates to see what is going on

    finimage=zeros(800,size(manip2_f,2));
    i=1;
    for j=1:1600
        finimage(j,:)=manip2_f(i,:);
        %%each pizel become  100 pixel
        resto= mod(j,200);
        if (resto==0)
            i=i+1;
        end
        
    end
    imshow(uint8(finimage));
end

