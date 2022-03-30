%##########################################################################
% ========================= Hammerstein & ARMAX ===========================
% ====================Presented by: Reza Saadatyar ========================
% ============================== 2017 =====================================
% Please Run Code Main.m, the Software will run automatically and there is no need to run the Functions folder codes
clc
clear
close all;
%% System: y(t)+a1y(t-1)+a2y(t-2)=b1u(t-1)+b2u(t-2)+v(t)+d1v(t-1)
N = 3000 ; % Itertion
x=zeros(N,1);x(1:2,1)=randn(2,1);
w1=zeros(N,1);w1(1:2,1)=randn(2,1);
w2=zeros(N,1);w2(1:2,1)=randn(2,1);
for i=1:2
if i==1
    va=0.5;
    v = va*randn(N,1); % noise
    for j=3:N
        w1(j)=1.6*w1(j-1)-0.8*w1(j-2)+v(j)-0.64*v(j-1);
    end
else
    va=2;
    v = va*randn(N,1); % noise
    for j=3:N
        w2(j)=1.6*w2(j-1)-0.8*w2(j-2)+v(j)-0.64*v(j-1);
    end
end
u=randn(N,1); % Input 
y = zeros(N,1) ; y(1:2,:) = 1 ; %output
teta=[1.6;0.8;0.85;0.425;0.2125;0.65;0.325;0.1625;0.64];%teta=[a1;a2;b1;c2b1;c3b1;b2;c2b2;c3b3;d1]
Nc=9 ; % number Parmeters
tetah=zeros(Nc,N);tetah(:,1)=randn(Nc,1) ;
Error=zeros(N,1);
PE=[];
P = 100 * eye(Nc);
for i = 3:N
    x(i)=1.6*x(i-1)-0.8*x(i-2)+0.85*u(i-1)+0.85*0.5*u(i-1)^2+...
     0.85*0.25*u(i-1)^3+0.65*u(i-2)+0.65*0.5*u(i-2)^2+0.65*0.25*u(i-2)^3;
    % y(t)+a1y(t-1)+a2y(t-2)=b1u(t-1)+b2u(t-2)+v(t)+d1v(t-1)
    y(i)=teta(1)*y(i-1)-teta(2)*y(i-2)+teta(3)*u(i-1)+teta(4)*u(i-1)^2+teta(5)*u(i-1)^3+...
    teta(6)*u(i-2)+ teta(7)*u(i-2)^2+ teta(8)*u(i-2)^3+v(i)-teta(9)*v(i-1);
    phi=[y(i-1), -y(i-2), u(i-1), u(i-1)^2, u(i-1)^3, u(i-2), u(i-2)^2, u(i-2)^3,-v(i-1)]';
    [tetah(: ,i) , P] = RLS(phi , y(i)  , tetah(: ,i-1) , P , Nc) ;
    % tetah
    if (i==100)||(i==200)||(i==300)||(i==500)||(i==1000)||(i==1500)||(i==2000)||(i==2500)||(i==3000)
    a1=tetah(1 ,i);a2=tetah(2 ,i);b1=tetah(3 ,i);b2=tetah(6 ,i);
    c2=tetah(4 ,i)/b1;c3=tetah(5 ,i)/b1;d1=tetah(9 ,i)/b1;
    num=norm(teta-tetah(:,i),2);
    dem=norm(teta,2);
    Error(i)=num/dem;
    E= Error(i);
    PE=[PE;i -a1 a2 b1 b2 c2 c3 -d1 E*100];%#ok
    end
    % Error
    num=norm(teta-tetah(:,i),2);
    dem=norm(teta,2);
    Error(i)=num/dem;
    Error(1:2,1)=Error(3);
end

T = array2table(PE,'VariableNames',{'t','a1','a2','b1','b2','c2','c3','d1','zeta'})%#ok

figure(1)
plot(Error , 'LineWidth' , 2) ;grid on ;hold on
ax=gca;ax.FontSize=10;ax.FontWeight='bold';
xlabel('t','FontWeight','bold','FontSize',14,'FontName','Times New Roman') ; 
ylabel('\delta','FontWeight','bold','FontSize',14,'FontName','Times New Roman') ; 
end
delta1=sqrt(var(w1)/var(x))*100;delta2=sqrt(var(w2)/var(x))*100;
legend({['\delta_n_s:' num2str(delta1),';  \sigma^2: 0.5'];[ '\delta_n_s:' num2str(delta2),';  \sigma^2: 2']},'FontWeight','bold','FontSize',12,'FontName','Times New Roman') ;
title('The parameter estimation errors \delta vs t')
xlim([0 800])