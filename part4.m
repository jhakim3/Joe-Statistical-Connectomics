function part4(patientnumbers,index,patientn,otimes,feat,beatq)

%This function creates the CO estimation using the estimateCO function
%provided, and uses the patient data ('n' files, in patientn input) to
%retrieve the thermodilution data to calibrate the estimation. We use
%method C2 for calibration.
%Inputs: Patient ID numbers for the graphs, the index of that list,
%onsettimes from wabp, features from abpfeature, and beatq is the
%beatquality assesment from jSQI.

[estimate5, timemin] = estimateCO_v3(otimes, feat, beatq, 5, 100); %using estimators id 5, 6 and 7
estimate6 = estimateCO_v3(otimes, feat, beatq, 6, 100);
estimate7 = estimateCO_v3(otimes, feat, beatq, 7, 100);
estimates = [estimate5, estimate6, estimate7]; %estimates(data index)(estimator index)

% Calibrate estimate against CO from thermodilution
ntable = readtable(patientn, 'Delimiter', 'tab');
COdouble = str2double(table2cell(ntable(2:end,16)));

% Get PP, MAP and HR from features file
PP = feat(:,5);
MAP = feat(:,6);
HR = feat(:,7);

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

% Finding the closest estimate time point to the first COtd measurement
n = 1;
min = intmax;
while abs(timemin(n) - COtdtimes(i)) < min
    min = abs(timemin(n) - COtdtimes(i));
    n = n+1;
end

%Scaling the data
for y = 1:3
    estimates(:,y) = estimates(:,y)*COtdvals(1)/estimates(1, y);
end

% Plotting
n = 1;
min = intmax;
while abs(timemin(n) - 720) < min
    min = abs(timemin(n) - 720);
    n = n+1;
end

%3x2 grid of plots

figure
% suptitle(sprintf('Patient %d part 4 graphs',patientnumbers(index)));
ax1 = subplot(3,2,1);
ax2 = subplot(3,2,2);
ax3 = subplot(3,2,3);
ax4 = subplot(3,2,4);
ax5 = subplot(3,2,5);
ax6 = subplot(3,2,6);

plot(ax1,timemin(1:n), estimates(1:n, 1));
hold(ax1, 'on');
plot(ax1, COtdtimes, COtdvals, 'o');
xlabel(ax1,'time (min)')
ylabel(ax1,'CO (L/min)')

plot(ax3, timemin(1:n), estimates(1:n, 2));
hold(ax3, 'on');
plot(ax3, COtdtimes, COtdvals, 'o');
xlabel(ax3,'time (min)')
ylabel(ax3,'CO (L/min)')

plot(ax5, timemin(1:n), estimates(1:n, 3));
hold(ax5, 'on');
plot(ax5, COtdtimes, COtdvals, 'o');
xlabel(ax5,'time (min)')
ylabel(ax5,'CO (L/min)')

plot(ax2,timemin(1:n), PP(1:n));
xlabel(ax2,'time (min)')
ylabel(ax2,'PP (mmHg)')

plot(ax4,timemin(1:n), MAP(1:n));
xlabel(ax4,'time (min)')
ylabel(ax4,'MAP (mmHg)')

plot(ax6,timemin(1:n), HR(1:n));
xlim(ax1,[0 720]);
axis(ax2,[0 720 10 120]);
xlim(ax3,[0 720]);
axis(ax4,[0 720 40 120]);
xlim(ax5,[0 720]);
axis(ax6,[0 720 40 120]);

ylimmax=max([ylim(ax1),ylim(ax3),ylim(ax5)]);

ylim(ax1, [0 ylimmax]);
ylim(ax3, [0 ylimmax]);
ylim(ax5, [0 ylimmax]);

xlabel(ax6,'time (min)')
ylabel(ax6,'Heart Rate (bpm)')