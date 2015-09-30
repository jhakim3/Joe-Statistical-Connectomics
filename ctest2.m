%%HERE IS WHERE YOU NEED TO ENTER THE COLUMNS OF NUM YOUR FINAL MODEL
%%INCLUDES!
inc_clin_cols = [4];
inc_wave_cols = [2,3,5,6];
%%%%%%%%%%%%%%%

%import clinical data - this is the basis of the static model
clinical_data = load('clinical_data_training.mat');
waveform_data = load('waveform_data_training.mat');

%extract variables for data structures
Fn = waveform_data.Fn;%non-septic waveform data
IDn = waveform_data.IDn;%non-septic pateint ID for each waveform time point
IDn_unique = unique(IDn); %list of non-septic patient ID numbers

Fs = waveform_data.Fs;%septic waveform data
IDs = waveform_data.IDs;%septic pateint ID for each waveform time point
IDs_unique = unique(IDs); %list of septic patient ID numbers

%generate covariates for complex model
Ybasic = clinical_data.num(:,2);%patient's septic/non-septic status
Xbasic = clinical_data.num(:,inc_clin_cols);%grabs some demographic data for model
IDbasic = clinical_data.num(:,1);%patient's ID numbers

%preparing a septic/non-septic vector for glmfit
%note all septic data will be first, followed by all non-septic data
Y_s = ones(length(IDs),1);
Y_n = zeros(length(IDn),1);

%initialize variables that will change in size
%memory inefficent but convenient
X_s = [];
X_n = [];

%create covariate matrix including both demographic info and waveform data
for k = 1:length(IDs_unique)%for septic data
    %tile demographic data to length of patients waveform data
    X_s_demo_add = repmat(Xbasic(IDbasic==IDs_unique(k) & Ybasic,:),sum(IDs == IDs_unique(k)),1);
    %grab some waveform data
    %make sure the rows grabbed here match equivalent line in next for loop
   
    
    X_s_wave_add = (Fs(2,IDs==IDs_unique(k)))'.^3;
    X_s_wave_add = [X_s_wave_add , Fs(3,IDs==IDs_unique(k)).^.15'];
    X_s_wave_add = [X_s_wave_add , Fs(5,IDs==IDs_unique(k)).^2'];
    X_s_wave_add = [X_s_wave_add , Fs(6,IDs==IDs_unique(k)).^7'];
       
    
    
    
    X_s_add = [X_s_demo_add X_s_wave_add];%combine demographic and wave data
    X_s = [X_s;X_s_add];%add current patient to growing covariate matrix 
end

for k = 1:length(IDn_unique)%for non-septic data
    X_n_demo_add = repmat(Xbasic(IDbasic==IDn_unique(k)& ~Ybasic,:) ,sum(IDn == IDn_unique(k)),1);
   
    
    X_n_wave_add = (Fn(2,IDn==IDn_unique(k)))'.^3;
    X_n_wave_add = [X_n_wave_add , Fn(3,IDn==IDn_unique(k)).^.15'];
    X_n_wave_add = [X_n_wave_add , Fn(5,IDn==IDn_unique(k)).^2'];
    X_n_wave_add = [X_n_wave_add , Fn(6,IDn==IDn_unique(k)).^7'];
    
    
    
    X_n_add = [X_n_demo_add X_n_wave_add];
    X_n = [X_n;X_n_add]; 
end

%training model
X = [X_n;X_s];
Y = [Y_n;Y_s];


X = [X, X(:,1).*X(:,3),X(:,1).*X(:,5)];

modelspec = 'y ~ x1 + x2 + x3 + x4 + x5 + x1:x3 + x1:x5';
md1 = fitglm(X,Y,modelspec,'Distribution','binomial')
B=md1.Coefficients.Estimate;

Phat = 1./(1+exp(-[ones(size(X,1),1) X]*B));
[thresh] = test_performance(Phat, Y);

% Found threshold for determining single data points in waveform, now
% scaling up to the patient level to find a threshold value for the patient
septic_wave = (Phat > thresh);
P_Y = [zeros(length(IDn_unique),1); ones(length(IDs_unique),1)];
Phat_P = [];

C = 1;
for i = 1:length(IDn_unique)
    a = sum(IDn == IDn_unique(i));
    atemp = septic_wave(C:C+a-1);
    C = C+a;
    percent_septic = length(find(atemp == 1))/length(atemp);
    Phat_P = [Phat_P; percent_septic];
end
for i = 1:length(IDs_unique)
    a = sum(IDs == IDs_unique(i));
    atemp = septic_wave(C:C+a-1);
    C = C+a;
    percent_septic = length(find(atemp == 1))/length(atemp);
    Phat_P = [Phat_P; percent_septic];
end

thresh_P = test_performance(Phat_P, P_Y);

Y_test_bestguess = Phat_P > thresh_P;
PercentCorrect = (1 - sum(abs(P_Y-Y_test_bestguess))/length(P_Y))*100;
disp(PercentCorrect);


%%%%%%%%%%%%%%%%%% Testing block %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%REQUIRE TESTING DATA





%bring in testing data
clinical_data_t = load('clinical_data_testing.mat');
waveform_data_t = load('waveform_data_testing.mat');

%extract variables for data structures
Fn_t = waveform_data_t.Fn;%non-septic waveform data
IDn_t = waveform_data_t.IDn;%non-septic pateint ID for each waveform time point
IDn_unique_t = unique(IDn_t); %list of non-septic patient ID numbers

Fs_t = waveform_data_t.Fs;%septic waveform data
IDs_t = waveform_data_t.IDs;%septic pateint ID for each waveform time point
IDs_unique_t = unique(IDs_t); %list of septic patient ID numbers

%generate covariates for complex model
Ybasic_t = clinical_data_t.num(:,2);%patient's septic/non-septic status
Xbasic_t = clinical_data_t.num(:,inc_clin_cols);%grabs some demographic data for model
IDbasic_t = clinical_data_t.num(:,1);%patient's ID numbers

%preparing a septic/non-septic vector for glmfit
%note all septic data will be first, followed by all non-septic data
Y_s_t = ones(length(IDs_t),1);
Y_n_t = zeros(length(IDn_t),1);

%initialize variables that will change in size
%memory inefficent but convenient
X_s_t = [];
X_n_t = [];

%create covariate matrix including both demographic info and waveform data
for k = 1:length(IDs_unique_t)
    X_s_demo_add_t = repmat(Xbasic_t(IDbasic_t==IDs_unique_t(k) & Ybasic_t,:),sum(IDs_t == IDs_unique_t(k)),1);
   
    
    
    X_s_wave_add = (Fs(2,IDs==IDs_unique(k)))'.^3;
    X_s_wave_add = [X_s_wave_add , Fs(3,IDs==IDs_unique(k)).^.15'];
    X_s_wave_add = [X_s_wave_add , Fs(5,IDs==IDs_unique(k)).^2'];
    X_s_wave_add = [X_s_wave_add , Fs(6,IDs==IDs_unique(k)).^7'];
       
    
    
    X_s_add_t = [X_s_demo_add_t X_s_wave_add_t];%combine demographic and wave data
    X_s_t = [X_s_t;X_s_add_t];%add current patient to growing covariate matrix 
end

for k = 1:length(IDn_unique_t)
    X_n_demo_add_t = repmat(Xbasic_t(IDbasic_t==IDn_unique_t(k)& ~Ybasic_t,:) ,sum(IDn_t == IDn_unique_t(k)),1);
    
    
    X_n_wave_add = (Fn(2,IDn==IDn_unique(k)))'.^3;
    X_n_wave_add = [X_n_wave_add , Fn(3,IDn==IDn_unique(k)).^.15'];
    X_n_wave_add = [X_n_wave_add , Fn(5,IDn==IDn_unique(k)).^2'];
    X_n_wave_add = [X_n_wave_add , Fn(6,IDn==IDn_unique(k)).^7'];
        
    
    X_n_add_t = [X_n_demo_add_t X_n_wave_add_t];
    X_n_t = [X_n_t;X_n_add_t]; 
end

X_test = [X_s_t; X_n_t];

X = [X, X(:,1).*X(:,3),X(:,1).*X(:,5)];

Phat_test = 1./(1+exp(-[ones(size(X_test,1),1) X_test]*B));

septic_wave = (Phat_test > thresh);
P_Y = [zeros(length(IDn_unique_t),1); ones(length(IDs_unique_t),1)];
Phat_P = [];

C = 1;
for i = 1:length(IDn_unique_t)
    a = sum(IDn_t == IDn_unique_t(i));
    atemp = septic_wave(C:C+a-1);
    C = C+a;
    percent_septic = length(find(atemp == 1))/length(atemp);
    Phat_P = [Phat_P; percent_septic];
end
for i = 1:length(IDs_unique_t)
    a = sum(IDs_t == IDs_unique_t(i));
    atemp = septic_wave(C:C+a-1);
    C = C+a;
    percent_septic = length(find(atemp == 1))/length(atemp);
    Phat_P = [Phat_P; percent_septic];
end

Phat_test = 1./(1+exp(-[ones(size(X_test,1),1) X_test]*B));
Y_test_bestguess = Phat_test>thresh_P;
PercentCorrect = (1 - sum(abs(P_Y-Y_test_bestguess))/length(Y_test))*100;

disp(PercentCorrect)