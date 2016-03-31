clear all
clc

results=dlmread('input_fint.res');

results(100:length(results(:,1)),:)=[];
results(:,2)=(results(:,2)/200)*7;
results(:,1)=results(:,1);

plot(results(:,2),results(:,1));