function [Y]=smoothJustDependencySilv(I,sigma,Smatrix)
% convolution of time
Y = (Smatrix*I')';
end