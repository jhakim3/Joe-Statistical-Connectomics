% This is the driver for the program
% The variables abpfiles and nfiles contain the file name list
% Usage: rename abpfiles and nfiles to local file location and run

% define location of abp and n files
%patientnums=[20,708];
patientnums=[20];

abpfiles = {  
    './s00020-2567-03-30-17-47_ABP.txt';
    './s00708-2567-03-30-17-47_ABP.txt';
    
};
nfiles = {
    './s00020-2567-03-30-17-47n.txt';...
    './s00020-2567-03-30-17-47n.txt';
    };


for x = 1:length(patientnums) %loop through patients 
    % load abp into matlab
    %%abpfile = load(abpfiles{x});
    ntable=readtable(nfiles{x}, 'Delimiter','tab');
    
    % process data
    %%abp = abpfile(:,2);
    %%otimes = wabp(abp);
    %%feat = abpfeature(abp, otimes);
    %%beatq = jSQI(feat, otimes, abp);
    
    %COvals=part1(feat, otimes, ntable);
    COvals=zeros(size(feat));
    part2(abp,feat, otimes, ntable,COvals);
 
end