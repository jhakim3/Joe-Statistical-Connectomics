function [co, to, told, fea] = code2()

ABPgross=load('D:\Downloads D\s00020-2567-03-30-17-47_ABP.txt');
abp=ABPgross(:,2);
wabpres=wabp(abp);
featres=abpfeature(abp,wabpres);
jsqires=jSQI(featres,wabpres,abp);
[co, to, told, fea]=estimateCO_v3(wabpres,featres,jsqires,5,0);