% This is the driver for the program
% The variables abpfiles and nfiles contain the file name list
% Usage: rename abpfiles and nfiles to local file location and run

% define location of abp and n files
patientnums=[20,708];

abpfiles = {  
    './s00020-2567-03-30-17-47_ABP.txt';
    './s00708-2868-08-02-12-43_ABP.txt';
    
};
nfiles = {
    './s00020-2567-03-30-17-47n.txt';...
    './s00708-2868-08-02-12-43n.txt';
    };


for index = 1:length(patientnums) %loop through patients 
    % load abp into matlab
    abpfile = load(abpfiles{index});
    ntable=readtable(nfiles{index}, 'Delimiter','tab');
    
    % process data
    abp = abpfile(:,2);
    otimes = wabp(abp);
    feat = abpfeature(abp, otimes);
    beatq = jSQI(feat, otimes, abp);
    
    
    [COVals,CVals]=part1(patientnums,index,abp,feat,otimes,ntable);
    part2(51,abp,feat,otimes,COVals,CVals)
end