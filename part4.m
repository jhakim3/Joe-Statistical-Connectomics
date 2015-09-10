sub20 = '/Users/hatle/Documents/Intro Comp Med/Patient Data/OneDrive/s00020-2567-03-30-17-47_ABP.txt';
nfile = '/Users/hatle/Documents/Intro Comp Med/Patient Data/OneDrive/s00020-2567-03-30-17-47n.txt';

% sub20 = 'D:\Users\Richard\Documents\s00020-2567-03-30-17-47_ABP.txt';
% nfile = 'D:\Users\Richard\Documents\s00020-2567-03-30-17-47n.txt';

addpath('/Users/hatle/Documents/Intro Comp Med/Project_Assignments/PhysioToolkitCardiacOutput_MatlabCode/2analyze');
addpath('/Users/hatle/Documents/Intro Comp Med/Project_Assignments/PhysioToolkitCardiacOutput_MatlabCode/3estimate');

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

%Get PP, MAP and HR from features file
PP = feat(:,5);
MAP = feat(:,6);
HR = feat(:,7);

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
n = 1;
min = intmax;
while abs(timemin(n) - 720) < min
    min = abs(timemin(n) - 720);
    n = n+1;
end

figure
ax1 = subplot(4,1,1);
ax2 = subplot(4,1,2);
ax3 = subplot(4,1,3);
ax4 = subplot(4,1,4);
plot(ax1,timemin(1:n), estimates(1:n));

xlabel(ax1,'time (min)')
ylabel(ax1,'Cardiac Output (L/min)')

plot(ax2,timemin(1:n), PP(1:n));

xlabel(ax2,'time (min)')
ylabel(ax2,'PP (mmHg)')

plot(ax3,timemin(1:n), MAP(1:n));

xlabel(ax3,'time (min)')
ylabel(ax3,'Mean Arterial Pressure (mmHg)')

plot(ax4,timemin(1:n), HR(1:n));
axis(ax1,[0 720 2 7])
axis(ax2,[0 720 20 120])
axis(ax3,[0 720 40 100])
axis(ax4,[0 720 50 100])
xlabel(ax4,'time (min)')
ylabel(ax4,'Heart Rate (bpm)')


