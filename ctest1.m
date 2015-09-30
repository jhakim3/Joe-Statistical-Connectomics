%%HERE IS WHERE YOU NEED TO ENTER THE COLUMNS OF NUM YOUR FINAL MODEL
%%INCLUDES!

%import clinical data - this is the basis of the static model
clinical_data = load('clinical_data_training.mat');

%generate your final static glm here (replace the code below with your
%model)



q=0.001;

Y = clinical_data.num(:,2);
X = [clinical_data.num(:,3) , clinical_data.num(:,4).^q];

%compute glm USING NEW FUNCTION 
modelspec = 'linear';
md1 = fitglm(X,Y,modelspec,'Distribution','binomial');

B=md1.Coefficients.Estimate;
% %construct phat from parameters and X 

%X = [X, X(:,1).*X(:,2)]; %since were using interactions, need to concatenate X

Phat = 1./(1+exp(-[ones(size(X,1),1) X]*B));

[thresh] = test_performance(Phat, Y);

%bring in testing data
clinical_data_test = load('clinical_data_testing.mat');

%modify line below to incoporate your final covariates
X_test = [clinical_data_test.num(:,3) , clinical_data_test.num(:,4).^q];

X_test = [X_test, X_test(:,1).*X_test(:,2)];
Phat_test = 1./(1+exp(-[ones(size(X_test,1),1) X_test]*B));
Y_test = clinical_data_test.num(:,2);

Y_test_bestguess = Phat_test>thresh;


PercentCorrect = (1 - sum(abs(Y_test-Y_test_bestguess))/length(Y_test))*100

%You can add code to output sensitivity and specificity below

