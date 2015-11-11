theta_in=30;

close all

%location of each pixel
nX = 200; %number of columns
nY = 200; %number of rows
xj = 1:nX; %x location of each column
yi = 1:nY; %y location of each row
[xij,yij] = meshgrid (xj,yi); % x,y location of each pixel

%generate landmarks
X = [125 100; 150 50];
%generate rotation variable
R = theta_in;%degrees;
%define target land marks by rotating X by rotation variable
for a = 1:length(X)
    magnitude = sqrt(X(a,1)^2 + X(a,2)^2);
    initialTheta = atand(X(a,2)/X(a,1));
    Y(a,1) = magnitude*cosd(initialTheta+R);
    Y(a,2) = magnitude*sind(initialTheta+R);
end

%Plot Landmarks and Grid
figure;
scatter(X(:,1),X(:,2),40,'MarkerEdgeColor','k',...
    'MarkerFaceColor','c',...
    'LineWidth',1.5);
hold on;
scatter(Y(:,1),Y(:,2),40,'MarkerEdgeColor','k',...
    'MarkerFaceColor','r',...
    'LineWidth',1.5);


%plot unstransformed grid
down = 10; %down sampling is important so you can see things clearly
xijdown = xij(1:down:end, 1:down:end);
yijdown = yij(1:down:end,1:down:end);
%this is a trick to plot a grid easily
%We actually plot a 3D surface.
%but view down directly from aboveso it lookes 2D
surf(xijdown,yijdown,ones(size(xijdown)),'facecolor','none','edgecolor','k');

%2.4 calculate optimal 2x2 matrix transformation
A = (((transpose(X)*X)^(-1))*(transpose(X)*Y)); %from last weeks homework; also makes sense because rotation should be of the form  [cos(theta) sin(theta);-sin(theat) cos(theta)];

%2.5
AX = X*A;
for i = 1:length(xij)
    for j = 1:length(yij)
        Aij = A*[xij(i,j); yij(i,j)];
        Axij(i,j) = Aij(2);
        Ayij(i,j) = Aij(1);
    end
end
figure;
scatter(X(:,1),X(:,2),40,'MarkerEdgeColor','k',...
    'MarkerFaceColor','c',...
    'LineWidth',1.5);
hold on;
scatter(Y(:,1),Y(:,2),40,'MarkerEdgeColor','k',...
    'MarkerFaceColor','r',...
    'LineWidth',1.5);
scatter(AX(:,1),AX(:,2),40,'MarkerEdgeColor','k',...
    'MarkerFaceColor','b',...
    'LineWidth',1.5);
Axijdown = Axij(1:down:end,1:down:end);
Ayijdown = Ayij(1:down:end,1:down:end);
surf(Axijdown,Ayijdown,ones(size(Axijdown)),'facecolor','none','edgecolor','k');

%2.6 Calculate the Jacobian
%This is 1 everywhere. Since the gradient of AX is just A, calculate as

det_Jacobian = ones(200,200).*det(A);

%2.7 Calculate an optimal Guassian kernel transformation



LHS=[(Y(1,:)-X(1,:))' ; (Y(2,:)-X(2,:))'];




function K_out = K(x,y)

id = [1 0;0 1];
K_out=id*exp((-1/5000)*((norm(x-y))^2));

end
K_mat=[K(X(1,:),X(1,:)) K(X(1,:),X(2,:));K(X(2,:),X(1,:)) K(X(2,:),X(2,:))];

P_Joe = K_mat\LHS;

%2.8

% initialize to identity, we will add the displacment v
phix = xij;
phiy = yij;

grad1=zeros(200,200);
grad2=zeros(200,200);
grad=zeros(200,200);
for i = 1 : nY
    for j = 1 : nX
        % add the displacement for each p(k) in the sum
        for k = 1 : size(X,2) % number of landmarks
            id = [1 0;0 1];
            Kij = K([j,i],X(k,:));
            % ... implement this, the kernel evaluated at (j,i) - X(k)
            temp=(Kij*[P_Joe(1);P_Joe(3)]) + (Kij*[P_Joe(2);P_Joe(4)]);
            phix(i,j) = phix(i,j) + temp(1); % ... add the x component for p(k)
            
            phiy(i,j) = phiy(i,j) + temp(2); % ... add the y component for p(k)
            
            grad1(i,j)=temp(1);
            grad2(i,j)=temp(2);
            
            %grad(i,j)=det([temp(1)+1,
            
        end
    end
end

figure;
scatter(X(:,1),X(:,2),40,'MarkerEdgeColor','k',...
    'MarkerFaceColor','c',...
    'LineWidth',1.5);
hold on;
scatter(Y(:,1),Y(:,2),40,'MarkerEdgeColor','k',...
    'MarkerFaceColor','r',...
    'LineWidth',1.5);

V_X = [([LHS(1),LHS(2)]+X(1,:));[LHS(3),LHS(4)]+X(2,:)];

scatter(V_X(:,1),V_X(:,2),40,'MarkerEdgeColor','k',...
    'MarkerFaceColor','b',...
    'LineWidth',1.5);

phixdown = phix(1:down:end,1:down:end);
phiydown = phiy(1:down:end,1:down:end);
surf(phixdown,phiydown,ones(size(phixdown)),'facecolor','none','edgecolor','k');


xij(2,2)

q=gradient(phix);
%2.9
for i = 1:nY-1
    for j = 1:nX-1
        %[a,b]=gradient(phix
        %grad(i,j)=det([phix(i,j+1)-phix(i,j),phix(i+1,j)-phix(i,j);phiy(i,j+1)-phiy(i,j),phiy(i+1,j)-phiy(i,j)]);
        %grad(i,j)=det([q(i+1,j+1),q(i+1,j);q(i,j+1),q(i,j)]);
    end
end

figure
contourf(grad,ones(200),200,'linestyle','none')
