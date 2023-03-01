function [stderrmean] = SEM(signal,dimension)

% Created by: Anya Krok
% Created on: 26 June 2018
% Description: calculate standard error of mean
%
% [sem] = SEM(signal,dimension)
%
% INPUT
%   signal: matrix to analyze
%   dimension: 1 or 2
%       if dim = 1, then SEM returns row vector containing SEM of the elements in each column
%       if dim = 2, then SEM returns column vector containing SEM of elements in each row
%
% OUTPUT
%   stderrmean: scaler or vector of standard error of mean along specified dimension

stderrmean = nanstd(signal,0,dimension)/sqrt(size(signal,dimension));

end