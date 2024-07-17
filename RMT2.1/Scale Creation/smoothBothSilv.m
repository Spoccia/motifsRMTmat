function [Y]=smoothBothSilv(I,sigma,Smatrix)
% convolution of time
J = imsmooth(I,sigma) ;
% convolution of dependency
Y = (Smatrix*J')';
end