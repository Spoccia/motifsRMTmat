function [Y]=computeEntropyScale(I,sigma)
% convolution of time
[time,variate] = size(I);
Y=zeros(time,variate);
for i=1:variate
    for j=1:time
        timeslice = round(j-3*sigma) : round(j+3*sigma);
        Y(j,i) = EntropySingVariate_mex(I(timeslice(timeslice>0 & timeslice<time),i),-inf);
    end
end
% J = imsmooth(I,sigma) ;
% % convolution of dependency
% Y = (Smatrix*J')';