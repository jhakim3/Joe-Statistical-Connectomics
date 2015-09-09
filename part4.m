sub20 = 'C:\Users\Steven\Documents\School Stuff\College\JUNIOR FALL 15\Computational Biology\Patient Data\s00020\s00020-2567-03-30-17-47_ABP.txt';
nfile = 'C:\Users\Steven\Documents\School Stuff\College\JUNIOR FALL 15\Computational Biology\Patient Data\s00020\s00020-2567-03-30-17-47n.txt';

% sub20 = 'D:\Users\Richard\Documents\s00020-2567-03-30-17-47_ABP.txt';
% nfile = 'D:\Users\Richard\Documents\s00020-2567-03-30-17-47n.txt';

addpath('C:\Users\Steven\Documents\MATLAB\School Stuff\Computational Medicine\PhysioToolkitCardiacOutput_MatlabCode\2analyze');
addpath('C:\Users\Steven\Documents\MATLAB\School Stuff\Computational Medicine\PhysioToolkitCardiacOutput_MatlabCode\3estimate');
addpath('C:\Users\Steven\Documents\MATLAB\School Stuff\Computational Medicine');

% Produce estimates without calibration
abpfile = importdata(sub20);
abp = abpfile(:,2);

otimes = wabp(abp);
feat = abpfeature(abp, otimes);
beatq = jSQI(feat, otimes, abp);
[estimates, timemin] = estimateCO_v3(otimes, feat, beatq, 5, 100);

% Calibrate estimates against CO from thermodilution
ntable = readtable(nfile, 'Delimiter', 'tab');
COdouble = str2double(table2cell(ntable(2:end,16)));

% Finding the first nonzero COtd value and its time
firstind = 1;
while COdouble(firstind) == 0
    firstind = firstind+1;
end
COtdtime = (firstind-1); %first nonzero value time (min)

% Finding the closest estimates time point to the first COtd measurement
n = 1;
min = intmax;
while abs(timemin(n) - COtdtime) < min
    min = abs(timemin(n) - COtdtime);
    n = n+1;
end
C2ratio = COdouble(firstind)/estimates(n-1);

% Performing calibration calculations and plotting
estimates = estimates*C2ratio;
