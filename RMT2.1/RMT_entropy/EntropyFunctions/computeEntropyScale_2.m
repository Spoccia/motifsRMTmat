function [Y]=computeEntropyScale_2(I,sigma,Smatrix,threshold)
% convolution of time
%threshold=0;
[time,variate] = size(I);
Y=zeros(time,variate);
for i=1:variate
     M_mask=Smatrix(i, :)>threshold;%0;
    for j=1:time
        timeslice = round(j-3*sigma) : round(j+3*sigma);
        %Y(j,i) = EntropySingVariate_mex(I(timeslice(timeslice>0 & timeslice<time),M_mask),-inf);
        Y(j,i) = entropy(uint8(I(timeslice(timeslice>0 & timeslice<time),M_mask))');
    end
end
% J = imsmooth(I,sigma) ;
% % convolution of dependency
% Y = (Smatrix*J')';