function [Y,S]=smooth(I,S,sigma)
% convolution of time
J = imsmooth(I,sigma) ;
% convolution of dependency
Y = (S*J')';
 

        