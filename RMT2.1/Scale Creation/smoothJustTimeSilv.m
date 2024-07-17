function [Y]=smoothJustTimeSilv(I,sigma,Smatrix)
% convolution of time
Y = imsmooth(I,sigma) ;
end