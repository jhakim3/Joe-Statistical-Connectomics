addpath('C:\Users\Steven\Documents\School Stuff\College\JUNIOR FALL 15\Computational Biology\Project_Assignments\PhysioToolkitCardiacOutput_MatlabCode\2analyze');
addpath('C:\Users\Steven\Documents\School Stuff\College\JUNIOR FALL 15\Computational Biology\Project_Assignments\PhysioToolkitCardiacOutput_MatlabCode\3estimate');
addpath('C:\Users\Steven\Documents\School Stuff\College\JUNIOR FALL 15\Computational Biology\Project_Assignments');

load('patient20.mat');

Tn = feat(11:end,7);
PP = feat(11:end,5);
MAP = feat(11:end,6);
b2bP = zeros(length(feat)-10,1);
for i = 1:length(feat)-10
    b2bP(i) = abp(otimes(i+11))-abp(otimes(i+10));
end

tau_naive = Tn.*MAP./(PP-b2bP);
m = 51; %window size
tau_window = zeros(length(tau_naive)-(m-1),1); 
for j = 1:length(tau_window)
    tau_window(j) = mean(tau_naive(j:j+m-1)); %mean or least squares?
end
plot(tau_naive);
figure;
plot(tau_window);

%finding C(n)
% Calibrate estimate against CO from thermodilution
ntable = readtable(patientn, 'Delimiter', 'tab');
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
min = 0.1; %minimum difference between the converted hertz time to whole minute times
for i = 1:length(MAP) 
    if n <= length(COtdMAPvector)
     if abs(otimes(i,1)/(125*60)) - COtdMAPvector(n,1)) < min 
        min = abs((feat(i,1)/(125*60)) - COtdMAPvector(n,1));
        COtdMAPvector(n,3)= feat(i,1); %place the time of feat (in hertz) in the 3rd column
        COtdMAPvector(n,4)= feat(i,7); %place the MAP from feat in the 4th column
        n = n+1;
     end
    end
end

tau_window(:,2) = (otimes(((m+1)/2):end-((m+1)/2),1));
for i = 1:length(tau_window) 
    if n <= length(COtdMAPvector)
     if abs(tau_window(i - COtdMAPvector(n,1) < min 
        min = abs(otimes(i,1)/(125*60)) - COtdMAPvector(n,1);
        COtdMAPvector(n,5)= feat(i,1); %place the time of feat (in hertz) in the 3rd column
        COtdMAPvector(n,6)= feat(i,7); %place the MAP from feat in the 4th column
        n = n+1;
     end
    end
end
%place tau in relation to COtds to find Ctds
for i = 1:length(tau_window)
    
end

X = [ones(length(COtdMAPvector),1) COtdMAPvector(:,4)];
gamma = X\COtdMAPvector(:,2); %Uses the  \ operator to perform the least squares mean regression to find Y intercept and the slope of the C(n) linear regression based off of the COtd datapoints

C=[];
for n = 1:length(MAP)
C(n,1) = gamma(1) + gamma(2)*MAP(n,1);
end



% adjustimg time scale for CO
% Tsample = zeros(length(Tn),1);
% for k = 1:length(Tsample)
%     Tsample(k) = sum(Tn(1:k));
% end
% Tsamplemin = Tsample/125/60;
% % Calibrate estimate5 against CO from thermodilution
% ntable = readtable(patientn, 'Delimiter', 'tab');
% COdouble = str2double(table2cell(ntable(2:end,16)));
% 
% % Finding the first nonzero COtd value and its time
% COtdtimescell = {};
% for i = 1:12*60
%     if COdouble(i) ~= 0
%         COtdtimescell{length(COtdtimescell)+1} = i-1;
%     end
% end
% COtdtimescell = COtdtimescell';
% 
% % Moving COtd data into numeric array so that it can be plotted
% COtdtimes = zeros(length(COtdtimescell),1);
% COtdvals = zeros(length(COtdtimescell),1);
% for i = 1:length(COtdtimescell)
%     COtdtimes(i) = COtdtimescell{i};
%     COtdvals(i) = COdouble(COtdtimes(i)+1);
% end

% CO = (b2bP(1+(m-1)/2:end-((m-1)/2))./Tn(1+(m-1)/2:end-((m-1)/2))) + (PP(1+(m-1)/2:end-((m-1)/2))./tau_window);

% C = mean(COtdvals)/mean(CO(COtdtimes));
% figure;
% plot(Tsamplemin(1+(m-1)/2:end-((m-1)/2)),CO);
% hold on;
% plot(COtdtimes,COtdvals,'ro');