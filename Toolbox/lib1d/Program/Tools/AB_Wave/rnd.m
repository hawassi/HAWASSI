function [y] = rnd(x,N)
% 
% function rounds with N-digits
% is used for bandpass filter (bandpass.m)

n = 10^N;
y = round(n*x)/n;