function [ output_args ] = statisticalAnalysis(otimes,feat,beatq)
%STATISTICALANALYSIS Summary of this function goes here
%   Detailed explanation goes here

estimate5 = estimateCO_v3(otimes,feat,beatq,5,100);
estimate14 = estimateCO_v3(otimes,feat,beatq,5,100);

error = estimate14-estimate5;

plot(error);

end

