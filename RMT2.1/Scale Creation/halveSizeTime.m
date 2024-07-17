function [ J ] = halveSizeTime( I )
%HALVESIZE Summary of this function goes here
%   Detailed explanation goes here

%% Shirnk the matrix by half on rows
J = I(1:2:end,:);

end