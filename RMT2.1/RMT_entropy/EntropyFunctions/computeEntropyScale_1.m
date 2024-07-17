function [Y]=computeEntropyScale_1(I,sigma,Smatrix,threshold)
% convolution of time
%threshold=0.015;
[time,variate] = size(I);
Y=zeros(time,variate);
for i=1:variate
     M_mask=Smatrix(i, :)>threshold;%0;
    for j=1:time
        timeslice = round(j-3*sigma) : round(j+3*sigma);
        subTS = I(timeslice(timeslice>0 & timeslice<time),M_mask);
        quantizesubTS= singlevariateQuantization(subTS );
        Y(j,i) = entropy(quantizesubTS/max(quantizesubTS(:)));%EntropySingVariate_mex(quantizesubTS,-inf);%EntropySingVariate_mex(subTS,-inf);%entropy(subTS/max(subTS(:)));%
    end
end
% J = imsmooth(I,sigma) ;
% % convolution of dependency
% Y = (Smatrix*J')';