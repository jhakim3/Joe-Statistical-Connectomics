% This is the driver for the program
% The variables abpfiles and nfiles contain the file name list
% Usage: rename abpfiles and nfiles to local file location and run

% define location of abp and n files
patientnums=[20,708,214,1502];

abpfiles = {  
    'D:\Downloads D\s00020-2567-03-30-17-47_ABP.txt';...
    'D:\Downloads D\s00708-2868-08-02-12-43_ABP.txt';...
    'D:\Downloads D\s00214-3084-11-28-16-23_ABP.txt';...
    'D:\Downloads D\s01502-2719-07-27-20-26_ABP.txt';
};
nfiles = {
    'D:\Downloads D\s00020-2567-03-30-17-47n.txt';...
    'D:\Downloads D\s00214-3084-11-28-16-23n.txt';
    };


for x = 1:length(abpfiles) %loop through patients 
    % load abp into matlab
    patientabp = abpfiles{x}; 
    abpfile = load(patientabp);
    
    % process data
    abp = abpfile(:,2);
    otimes = wabp(abp);
    feat = abpfeature(abp, otimes);
    beatq = jSQI(feat, otimes, abp);
    
    % runs part 3
    part3(patientnums,x,abp,otimes,feat);
    
    % runs only first 2 patients for part 4
    if x == 1 || x == 2
        patientn = nfiles{x};
        part4(patientnums,x,patientn,otimes,feat,beatq);
    end
 
end