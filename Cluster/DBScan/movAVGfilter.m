function [ filtered ] = movAVGfilter(data,windowSize)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
b = (1/windowSize)*ones(1,windowSize);
a = 1;
paddata = padarray(data,windowSize-1,'pre','replicate');
filtered = filter(b,a,paddata);
filtered = filtered(windowSize+floor(windowSize/2):end);
 end
% value = floor(windowSize/2);
% paddata = padarray(data,value,'both','replicate');
% filtered = filter(b,a,paddata);
% filtered = filtered(value:end-value);
