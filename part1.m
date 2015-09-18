function [COVals, CVals] = part1(patientnums,index,abp,feat,otimes,ntable)

Tn = feat(:,7);
MAP = feat(:,6);
% PP = 2*(MAP-feat(:,5));
PP = feat(:,5);

delP = zeros(length(feat),1);
for i = 1:length(feat)
    delP(i) = abp(otimes(i+1))-abp(otimes(i));
end

tau_naive = Tn.*MAP./(PP-delP);
win = 1001; % tau estimate window size
tau_window = zeros(length(tau_naive)-(win-1),1);
for i = 1:length(tau_window)
    tau_window(i) = mean(tau_naive(i:i+win-1));
end

%finding C(n)
% Calibrate estimate against CO from thermodilution
COdouble = str2double(table2cell(ntable(2:end,16)));

% Finding the first nonzero COtd value and its time
COtdtimescell = {};
for i = 1:12*60
    if COdouble(i) ~= 0
        COtdtimescell{length(COtdtimescell)+1} = i-1;
    end
end

% Moving COtd data into numeric array so that it can be plotted
COtdtimes = zeros(length(COtdtimescell));
COtdvals = zeros(length(COtdtimescell));
for i = 1:length(COtdtimescell)
    COtdtimes(i) = COtdtimescell{i};
    COtdvals(i) = COdouble(COtdtimes(i)+1);
end

COtdMAPvector= [];
COtdMAPvector(:,1) = COtdtimes(:,1);
COtdMAPvector(:,2) = COtdvals(:,1);
n=1;
min = .1; %minimum difference between the converted beat time to whole minute times
for i = 11:length(otimes) 
    if n <= length(COtdMAPvector)
     if abs(((otimes(i,1))/(125*60)) - COtdMAPvector(n,1)) < min 
       
        COtdMAPvector(n,3)= otimes(i,1); %place the onset time (in minutes) in the 3rd column
        COtdMAPvector(n,4)= feat(i,7); %place the MAP from feat in the 4th column
        n = n+1;
        %min = big number again
     end
    end
end
min = 20;
for n = 1:length(otimes)
    otimes(n,2) = (n);
end
n=1;
tau_window(:,2) = (otimes(((win+1)/2):(end-1-((win-1)/2)),1));
for i = 1:length(tau_window) 
    if n <= size(COtdMAPvector,1)
        if abs(tau_window(i,2) - COtdMAPvector(n,3)) < min 
            COtdMAPvector(n,5)= tau_window(i,2); %place the time of window (in otimes) in the 5th column
            COtdMAPvector(n,6)= tau_window(i,1); %place the tau from tau_window in the 6th column
            COtdMAPvector(n,7) = otimes(((win-1)/2)+i,2); %beat number
            n = n+1;
        end

    end
end

%find Ctds
for n = 1:size(COtdMAPvector,1)
    beatnum = COtdMAPvector(n,7);
COtdMAPvector(n,8) = COtdMAPvector(n,2)/((delP(beatnum)/Tn(beatnum)) + (MAP(beatnum)/COtdMAPvector(n,6))); 
end
X = [ones(size(COtdMAPvector,1),1) COtdMAPvector(:,8)];
gamma = X\COtdMAPvector(:,8); %Uses the  \ operator to perform the least squares mean regression to find Y intercept and the slope of the C(n) linear regression based off of the COtd datapoints

CVals=zeros(length(tau_window),1);
for n = 1:length(tau_window)
CVals(n) = gamma(1) + gamma(2)*MAP(n+(win-1)/2,1);
end

%Find COestimate
COVals = zeros(length(tau_window),1);
for n = 1:length(tau_window)
COVals(n) = CVals(n)*((delP(n+((win-1)/2))/Tn(n+((win-1)/2))) + (MAP(n+((win-1)/2))/tau_window(n)));
end

% Plotting Parlikar
timemin = otimes(1+((win-1)/2):end-((win-1)/2)-1,1)/125/60;
n = 1;
min = intmax;
while abs(timemin(n) - 720) < min
    min = abs(timemin(n) - 720);
    n = n+1;
end
figure;
suptitle(sprintf('Patient %d Parlikar Estimator',patientnums(index)));
ax1 = subplot(3,1,1);
ax2 = subplot(3,1,2);
ax3 = subplot(3,1,3);

plot(ax1,timemin(1:n), COVals(1:n));
hold(ax1, 'on');
plot(ax1, COtdtimes, COtdvals, 'o');
xlabel(ax1,'time (min)');
ylabel(ax1,'CO (L/min)');
axis(ax1,[0 720 -100 500]);
hold off;

plot(ax2,timemin(1:n), PP(1:n));
hold(ax2, 'on');
xlabel(ax2,'time (min)');
ylabel(ax2,'PP (mmHg)');
axis(ax2,[0 720 10 120]);
hold off;

plot(ax3,timemin(1:n), MAP(1:n));
hold(ax3, 'on');
xlabel(ax3,'time (min)');
ylabel(ax3,'MAP (mmHg)');
axis(ax3,[0 720 40 120]);
hold off;

% Plotting Liljestrand
[estimate5,min] = estimateCO_v3(otimes, feat, beatq, 5, 100);
estimate5 = estimate5*COtdMAPvector(1,2)./estimate5(find(min>COtdMAPvector(1,1),1));
figure;
suptitle(sprintf('Patient %d Liljestrand Estimator',patientnums(index)));
ax1 = subplot(3,1,1);
ax2 = subplot(3,1,2);
ax3 = subplot(3,1,3);

plot(ax1,timemin(1:n), estimate5(1:n));
hold(ax1, 'on');
plot(ax1, COtdtimes, COtdvals, 'o');
xlabel(ax1,'time (min)');
ylabel(ax1,'CO (L/min)');
xlim(ax1,[0 720]);
hold off;

plot(ax2,timemin(1:n), PP(1:n));
hold(ax2, 'on');
xlabel(ax2,'time (min)');
ylabel(ax2,'PP (mmHg)');
axis(ax2,[0 720 10 120]);
hold off;

plot(ax3,timemin(1:n), MAP(1:n));
hold(ax3, 'on');
xlabel(ax3,'time (min)');
ylabel(ax3,'MAP (mmHg)');
axis(ax3,[0 720 40 120]);
hold off;

% Plotting Liljestrand/Parlikar Estimator
figure;
% suptitle(sprintf('Patient %d',patientnums(index)));
ax1 = subplot(2,1,1);
ax2 = subplot(2,1,2);


plot(ax1,timemin(1:n), COVals(1:n));
hold(ax1, 'on');
plot(ax1, COtdtimes, COtdvals, 'o');
xlabel(ax1,'time (min)');
ylabel(ax1,'CO (L/min)');
axis(ax1,[0 720 -100 500]);
title(ax1,'Parlikar');
hold off;

plot(ax2,timemin(1:n), estimate5(1:n));
hold(ax2, 'on');
plot(ax2, COtdtimes, COtdvals, 'o');
xlabel(ax2,'time (min)');
ylabel(ax2,'CO (L/min)');
xlim(ax2,[0 720]);
title(ax2,'Liljestrand');
hold off;

end