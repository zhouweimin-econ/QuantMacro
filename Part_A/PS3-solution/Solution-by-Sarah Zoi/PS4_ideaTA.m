close all
clear all
clc
% Heterogeneity in the initial endowment and in the second period
% labor productivity
% Initialize parameters
beta=0.99;
r=0.02;
sigma=3;
k=4;
nu=4;
eps=0.05;
% initial vector of labor productivities
eta=[1 1.5 2.5 3];
% initial distribution of wealth
n=100;
% choose 0 for asymmetric distribution, 1 for uniform and 2 for pareto
wd=0;
% choose initial wealth distribution
if wd==0
y0=linspace(0.001,0.009,n);
y0(y0>=0.0055 & y0<=0.0087)=0.001;

elseif wd==1
    y0=linspace(0.001/2,0.009/2,n);
elseif wd==2
    f=1/(length(eta)*n);
  y0=f*ones(1,n);  
end
y0=sort(y0);
% figure
% hist(y0,25)
% title('Distribution of initial wealth')
% xlabel('y_{0}')
%% Solving the economy (General Equilibrium) without taxes and with flat tax rate
T1=0;
T2=0;
tau=0;
%if want the case with taxes put tau1 different from zero
tau1=0.115;
options=optimoptions('fsolve','MaxFunEvals',200000,'MaxIter',100000);
global optch rev
optch=nan(n,3,length(eta));
rev=nan(n,2,length(eta));


                r0=0.005;
                rstar=fsolve(@(r) totalassets(y0,r,eta,sigma,beta,tau,T1,T2,eps,k,nu),r0,options);
                r=rstar;
                C0=nan(n,1,length(eta));
                C1=nan(n,2,length(eta));
                %consumption growth
                cg=nan(n,3,length(eta));
                SR=nan(n,1,length(eta));
                %store labor income before and after taxes
                wh0=nan(n,3,length(eta));
                wh1=nan(n,6,length(eta));
                % before tax labor income growth
                whg=nan(n,3,length(eta));
                U=nan(n,2,length(eta));
                
                atot=sum(sum(optch(:,1,:)));
                lc=[0 0 1.0;0 0 0;0.8 0 0;0 0.8 0.2];
                
                for i=1:length(eta)
                    % compute consumption and store it
                C0(:,:,i)=y0'-optch(:,1,i)+(1-tau)*eta(i)*optch(:,2,i)+T1;
                % c1 high
                C1(:,1,i)=optch(:,1,i)*(1+r)+(1-tau)*(eta(i)+eps)*optch(:,3,i)+T2;
                % c1 low
                C1(:,2,i)=optch(:,1,i)*(1+r)+(1-tau)*(eta(i)-eps)*optch(:,3,i)+T2;
                % consumption growth
                cg(:,1,i)=(C1(:,1,i)-C0(:,:,i))./C0(:,:,i);
                cg(:,2,i)=(C1(:,2,i)-C0(:,:,i))./C0(:,:,i);
                cg(:,3,i)=(cg(:,1,i)+cg(:,2,i))/2 ;
                % utility today
                U(:,1,i)=( (C0(:,:,i)).^(1-sigma) )./(1-sigma) - k* ((optch(:,2,i)).^(1+1/nu) )./(1+1/nu);
                % discounted expected utility tomorrow
                U(:,2,i)=beta* 0.5* [ ( (C1(:,1,i)).^(1-sigma) )./(1-sigma) - k* ((optch(:,3,i)).^(1+1/nu) )./(1+1/nu)]+...
                         beta* 0.5* [ ( (C1(:,2,i)).^(1-sigma) )./(1-sigma) - k* ((optch(:,3,i)).^(1+1/nu) )./(1+1/nu)]
                % Saving rate
                SR(:,1,i)=optch(:,1,i)./y0';
                 % pretax labor income t=0
                wh0(:,1,i)=optch(:,2,i)*eta(i);
                % after-tax labor income t=0 
                wh0(:,2,i)=optch(:,2,i)*eta(i)*(1-tau);
                % total after-tax income t=0
                wh0(:,3,i)=wh0(:,2,i)+y0'+T1;
                
                % pretax labor income t=1 eps=+
                wh1(:,1,i)=optch(:,3,i)*(eta(i)+eps);
                % after-tax labor income t=1 eps=+
                wh1(:,2,i)=optch(:,3,i)*(eta(i)+eps)*(1-tau);
                % total after-tax income t=1 eps=+
                wh1(:,3,i)=wh1(:,2,i)+optch(:,1,i)*(1+r)+T2;
                % pretax labor income t=1 eps=-
                wh1(:,4,i)=optch(:,3,i)*(eta(i)-eps);
                % after-tax labor income t=1 eps=-
                wh1(:,5,i)=optch(:,3,i)*(eta(i)-eps)*(1-tau);
                % total after-tax income t=1 eps=-
                wh1(:,6,i)=wh1(:,5,i)+optch(:,1,i)*(1+r)+T2;
                
                % before tax labor income growth high shock
                whg(:,1,i)=(wh1(:,1,i)-wh0(:,1,i))./wh0(:,1,i);
                % before tax labor income growth low shock
                whg(:,2,i)=(wh1(:,4,i)-wh0(:,1,i))./wh0(:,1,i);
                % expected before tax labor income 
                whg(:,3,i)=(whg(:,1,i)+whg(:,2,i))/2 ;
                figure(i)
                
                
                lci=lc(i,:);
                subplot(3,3,1)
                
                plot(y0,optch(:,1,i),'Color',[lci])
                grid on
                axis tight
                xlabel('y_{0}')
                ylabel('a^{*}')

                subplot(3,3,2)
                title(['\eta = ' num2str(eta(i)) ],'FontSize', 12)
                plot(y0,optch(:,2,i),'Color',[lci])
                grid on
                axis tight
                xlabel('y_{0}')
                ylabel('h_{0}^{*}')
                title(['\eta = ' num2str(eta(i)) ' \epsilon = ' num2str(eps) ...
                    ' r_{\tau=0}^{*} = ' num2str(round(r,3))],'FontSize', 14)


                subplot(3,3,3)
                plot(y0,optch(:,3,i),'Color',[lci])
                grid on
                axis tight
                xlabel('y_{0}')
                ylabel('h_{1}^{*}')


                subplot(3,3,4)
                plot(y0,C0(:,:,i),'Color',[lci])
                grid on
                axis tight
                xlabel('y_{0}')
                ylabel('c_{0}^{*}')

                subplot(3,3,5)
                plot(y0,C1(:,1,i),'Color',[lci])
                hold on
                plot(y0,C1(:,2,i),'Color',[0.8 0.6 0])
                grid on
                axis tight
                xlabel('y_{0}')
                ylabel('c_1^{*}')
                legend ('\epsilon_{H}','\epsilon_{L}')
                
                subplot(3,3,6)
                plot(y0,cg(:,1,i),'Color',[lci])
                hold on
                plot(y0,cg(:,2,i),'Color',[0.8 0.6 0])
                grid on
                axis tight
                xlabel('y_{0}')
                ylabel('g_{C}')
                
                subplot(3,3,7)
                plot(y0,SR(:,1,i),'Color',[lci])
                grid on
                axis tight
                xlabel('y_{0}')
                ylabel('a^{*}/y_{0}')
                
                subplot(3,3,8)
                plot(y0,whg(:,1,i),'Color',[lci])
                hold on
                plot(y0,whg(:,2,i),'Color',[0.8 0.6 0])
                grid on
                axis tight
                xlabel('y_{0}')
                ylabel('g_{wh}')
                
                subplot(3,3,9)
                plot(y0,cg(:,3,i)./whg(:,3,i),'Color',[lci])               
                grid on
                axis tight
                xlabel('y_{0}')
                ylabel('E[g_{c}]/E[g_{wh}]')
                
                
                end
                
 % welfare function, as function of initial wealth (across productivities)
 W0=sum(U(:,1,:),3);
 W1=sum(U(:,2,:),3);
 figure (length(eta)+1)
set(0,'DefaultAxesColorOrder',[ 0 0 1.0;0 0 0;0.8 0 0;0 0.8 0.2],...
                    'defaultLineLineWidth',1.3,'DefaultAxesLineStyleOrder','-|-.|:') 
 subplot(2,3,1)
 plot(y0,permute(U(:,1,:),[1,3,2]))
 grid on
 axis tight
 xlabel('y_{0}')
 title('U (c_{0},h_{0})')
 legend(['\eta = ' num2str(eta(1))],['\eta = ' num2str(eta(2))],...
                    ['\eta = ' num2str(eta(3))],['\eta = ' num2str(eta(4))])
 subplot(2,3,2)
 plot(y0,permute(U(:,2,:),[1,3,2]))
 grid on
 axis tight
 title('\beta* E [U (c_{1},h_{1})]')
 subplot(2,3,3)
 plot(y0,permute(U(:,1,:),[1,3,2])+permute(U(:,2,:),[1,3,2]))
 grid on
 axis tight
 xlabel('y_{0}')
 title('V(y_{0},\eta)')
 subplot(2,3,4)
 plot(y0,W0)
 grid on
 axis tight
 xlabel('y_{0}')
 title('W_{0}(y_{0})')
 subplot(2,3,5)
 plot(y0,W1)
 grid on
 axis tight
 xlabel('y_{0}')
 title('E [W_{1}(y_{0})]')
 subplot(2,3,6)
 plot(y0,W0+W1)
 grid on
 axis tight
 xlabel('y_{0}')
 title('W (y_{0})')
 
                [out]=graph1(y0,optch, C0, C1, SR, tau, eps, eta, r)
                [out]=graph2(y0,wh0,wh1,tau, eps, eta, r, T1, T2)
                [out]=graph3(y0,cg,whg, tau, eps, eta, r, T1, T2)
                display('--------------------------------------------------------')
                display(['The equilibrium interest rate is ' num2str(r)])
                display(['At this interest rate net assets are ' num2str(atot)])
                display('--------------------------------------------------------')
rtau0=r;
results=[rtau0 T1 T2];
savtau0=sum(sum(optch(optch(:,1,:)>=0)));
savings=[savtau0];
%Lshare has pretax LS in t=0, after tax LS in t=0, pretax LS in t=1, after tax LS in t=1
[SRaggtau0, Lsharetau0, Htau0]=aggstat(y0,optch,eta,tau,r,T1,T2);
%% General Equilibrium with flat tax rate
tau=tau1;
                r0=0.1;
                T10=0.2;
                T20=0.2;
                rT0=[r0 T10 T20];
                %realization of the shock
                shockhit=rand(n,length(eta));
                % index of people hit by the shock
                shockhit=(shockhit<=0.5);
                etarealized=repmat(eta,n,1).*ones(n,length(eta))+eps*shockhit-(ones(n,length(eta))-shockhit)*eps;
                
                
                optch=nan(n,3,length(eta));
                rTstar=fsolve(@(rT) GE(y0,rT(1),eta,sigma,beta,tau,rT(2),rT(3),eps,k,nu,etarealized),rT0,options);
                r=rTstar(1);
                T1=rTstar(2);
                T2=rTstar(3);
                C0=nan(n,1,length(eta));
                C1=nan(n,2,length(eta));
                SR=nan(n,1,length(eta));
                %store labor income before and after taxes
                wh0=nan(n,3,length(eta));
                wh1=nan(n,6,length(eta));
                % before tax labor income growth
                whg=nan(n,3,length(eta));
                U=nan(n,2,length(eta));
                %consumption growth
                cg=nan(n,3,length(eta));
                lc=[0 1 1; 0.6 0.8 0.8; 1 0.2 1; 0 0.4 0.4];
                atot=sum(sum(optch(:,1,:)));
                govbal1=sum(sum(rev(:,1,:)))-T1*n*length(eta);
                govbal2=sum(sum(rev(:,2,:)))-T2*n*length(eta);
                % realized labor income in the second period
                whreal=nan(n,3,length(eta));
                for i=1:length(eta)
                    %before tax
                    whreal(:,1,i)=optch(:,3,i).*etarealized(:,i);
                    %after tax
                    whreal(:,2,i)=whreal(:,1,i)*(1-tau);
                    % total income after tax
                    whreal(:,3,i)=whreal(:,2,i)+T2+optch(:,1,i)*(1+r);
                end
                
                
                %govbal=sum(sum(rev(:,1,:)))-T*n*length(eta)+(sum(sum(rev(:,2,:)))-T*n*length(eta))/(1+r);
                for i=1:length(eta)
                    % compute consumption and store it
                C0(:,:,i)=y0'-optch(:,1,i)+(1-tau)*eta(i)*optch(:,2,i)+T1;
                % c1 high
                C1(:,1,i)=optch(:,1,i)*(1+r)+(1-tau)*(eta(i)+eps)*optch(:,3,i)+T2;
                % c1 low
                C1(:,2,i)=optch(:,1,i)*(1+r)+(1-tau)*(eta(i)-eps)*optch(:,3,i)+T2;
                % consumption growth
                cg(:,1,i)=(C1(:,1,i)-C0(:,:,i))./C0(:,:,i);
                cg(:,2,i)=(C1(:,2,i)-C0(:,:,i))./C0(:,:,i);
                cg(:,3,i)= (cg(:,1,i)+cg(:,2,i))/2 ;
                % utility today  
                U(:,1,i)=( (C0(:,:,i)).^(1-sigma) )./(1-sigma) - k* ((optch(:,2,i)).^(1+1/nu) )./(1+1/nu);
                % discounted expected utility tomorrow
                U(:,2,i)=beta* 0.5* [ ( (C1(:,1,i)).^(1-sigma) )./(1-sigma) - k* ((optch(:,3,i)).^(1+1/nu) )./(1+1/nu)]+...
                         beta* 0.5* [ ( (C1(:,2,i)).^(1-sigma) )./(1-sigma) - k* ((optch(:,3,i)).^(1+1/nu) )./(1+1/nu)];
                % Saving rate
                SR(:,1,i)=optch(:,1,i)./y0';
                % pretax labor income t=0
                wh0(:,1,i)=optch(:,2,i)*eta(i);
                % after-tax labor income t=0 
                wh0(:,2,i)=optch(:,2,i)*eta(i)*(1-tau);
                % total after-tax income t=0
                wh0(:,3,i)=wh0(:,2,i)+y0'+T1;
                
                
                % pretax labor income t=1 eps=+
                wh1(:,1,i)=optch(:,3,i)*(eta(i)+eps);
                % after-tax labor income t=1 eps=+
                wh1(:,2,i)=optch(:,3,i)*(eta(i)+eps)*(1-tau);
                % total after-tax income t=1 eps=+
                wh1(:,3,i)=wh1(:,2,i)+optch(:,1,i)*(1+r)+T2;
                % pretax labor income t=1 eps=-
                wh1(:,4,i)=optch(:,3,i)*(eta(i)-eps);
                % after-tax labor income t=1 eps=-
                wh1(:,5,i)=optch(:,3,i)*(eta(i)-eps)*(1-tau);
                % total after-tax income t=1 eps=-
                wh1(:,6,i)=wh1(:,5,i)+optch(:,1,i)*(1+r)+T2;
                
                % before tax labor income growth high shock
                whg(:,1,i)=(wh1(:,1,i)-wh0(:,1,i))./wh0(:,1,i);
                % before tax labor income growth low shock
                whg(:,2,i)=(wh1(:,4,i)-wh0(:,1,i))./wh0(:,1,i);
                % expected before tax labor income 
                whg(:,3,i)=(whg(:,1,i)+whg(:,2,i))/2 ;
                
                figure(i)
                lci=lc(i,:);
                subplot(3,3,1)
                hold on
                plot(y0,optch(:,1,i),'Color',[lci])
                grid on
                axis tight
                xlabel('y_{0}')
                ylabel('a^{*}')
                legend('\tau = 0',['\tau = ' num2str(tau)])
                
                subplot(3,3,2)
                hold on
                title(['\eta = ' num2str(eta(i)) ],'FontSize', 14)
                plot(y0,optch(:,2,i),'Color',[lci])
                grid on
                axis tight
                xlabel('y_{0}')
                ylabel('h_{0}^{*}')
                title({['\eta = ' num2str(eta(i)) ' \epsilon= ' num2str(eps) ' r_{\tau=0}^{*} = ' num2str(round(rtau0,3))]; ...
                    [' r^{*}_{\tau=0.115}= ' num2str(round(r,3)) ' T1_{\tau=0.115}^{*}= ' num2str(round(T1,3)) ' T2_{\tau=0.115}^{*}= ' num2str(round(T2,3))]},'FontSize', 12)


                subplot(3,3,3)
                hold on
                plot(y0,optch(:,3,i),'Color',[lci])
                grid on
                axis tight
                xlabel('y_{0}')
                ylabel('h_{1}^{*}')


                subplot(3,3,4)
                hold on
                plot(y0,C0(:,:,i),'Color',[lci])
                grid on
                axis tight
                xlabel('y_{0}')
                ylabel('c_{0}^{*}')

                subplot(3,3,5)
                hold on
                plot(y0,C1(:,1,i),'Color',[lci])
                hold on
                plot(y0,C1(:,2,i),'Color',[0.6 1 0.2])
                grid on
                axis tight
                xlabel('y_{0}')
                ylabel('c_1^{*}')
                legend ('\tau = 0, \epsilon_{H}','\tau = 0, \epsilon_{L}',...
                    ['\tau = ' num2str(tau) ' \epsilon_{H}'],['\tau = ' num2str(tau) ' \epsilon_{L}'])
                
                subplot(3,3,6)
                plot(y0,cg(:,1,i),'Color',[lci])
                hold on
                plot(y0,cg(:,2,i),'Color',[0.6 1 0.2])
                grid on
                axis tight
                xlabel('y_{0}')
                ylabel('g_{c}')
                
                subplot(3,3,7)
                hold on
                plot(y0,SR(:,1,i),'Color',[lci])
                grid on
                axis tight
                xlabel('y_{0}')
                ylabel('a^{*}/y_{0}')
                
                subplot(3,3,8)
                hold on
                plot(y0,whg(:,1,i),'Color',[lci])
                hold on
                plot(y0,whg(:,2,i),'Color',[0.6 1 0.2])
                grid on
                axis tight
                xlabel('y_{0}')
                ylabel('g_{wh}')
                
                subplot(3,3,9)
                hold on
                plot(y0,cg(:,3,i)./whg(:,3,i),'Color',[lci])               
                grid on
                axis tight
                xlabel('y_{0}')
               ylabel('E[g_{c}]/E[g_{wh}]')
                end
                % welfare function, as function of initial wealth (across productivities)
 
 W0=sum(U(:,1,:),3);
 W1=sum(U(:,2,:),3);
 figure (length(eta)+1)
 %set(0,'DefaultAxesColorOrder',[ 0 1 1; 0.6 0.8 0.8; 1 0.2 1; 0 0.4 0.4;0 1 1; 0.6 0.8 0.8; 1 0.2 1; 0 0.4 0.4],...
 %                   'defaultLineLineWidth',1.3,'DefaultAxesLineStyleOrder','-')
 subplot(2,3,1)
 hold on
 plot(y0,permute(U(:,1,:),[1,3,2]))
 legend(['\eta = ' num2str(eta(1)) ', \tau =0'],['\eta = ' num2str(eta(2)) ', \tau =0'],...
        ['\eta = ' num2str(eta(3)) ', \tau =0'],['\eta = ' num2str(eta(4)) ', \tau =0'],...
        ['\eta = ' num2str(eta(1)) ', \tau =' num2str(tau)],['\eta = ' num2str(eta(2)) ', \tau =' num2str(tau)],...
        ['\eta = ' num2str(eta(3)) ', \tau =' num2str(tau)],['\eta = ' num2str(eta(4)) ', \tau =' num2str(tau)])
 subplot(2,3,2)
 hold on
 plot(y0,permute(U(:,2,:),[1,3,2]))
 subplot(2,3,3)
 hold on
 plot(y0,permute(U(:,1,:),[1,3,2])+permute(U(:,2,:),[1,3,2]))
 subplot(2,3,4)
 hold on
 plot(y0,W0,'Color',[0 1 1])
 legend('\tau = 0',['\tau =' num2str(tau)])
 subplot(2,3,5)
 hold on
 plot(y0,W1,'Color',[0 1 1])
 subplot(2,3,6)
 hold on
 plot(y0,W0+W1,'Color',[0 1 1])
 
[out]=graph1(y0,optch, C0, C1, SR, tau, eps, eta, r, T1, T2)
[out]=graph2(y0,wh0,wh1,tau, eps, eta, r, T1,T2)
[out]=graph3(y0,cg,whg, tau, eps, eta, r, T1,T2)
                display('--------------------------------------------------------')
                display(['The equilibrium interest rate for tau= ' num2str(tau) ' is ' num2str(r)])
                display(['The equilibrium Transfer in t=1 for tau= ' num2str(tau) ' is ' num2str(T1)])
                display(['The equilibrium Transfer in t=2 for tau= ' num2str(tau) ' is ' num2str(T2)])
                display(['At this interest rate and Transfers net assets are ' num2str(atot)])
                display(['At this interest rate and Transfers governement budget in period 1 is ' num2str(govbal1)])
                display(['At this interest rate and Transfers governement budget in period 2 is ' num2str(govbal2)])
                %display(['The intertemporal governement budget is ' num2str(govbal)])
                display('--------------------------------------------------------')
rtauflat=r;
Ttauflat1=T1;
Ttauflat2=T2;
    results=[results; rtauflat Ttauflat1 Ttauflat2];  
    savtauflat=sum(sum(optch(optch(:,1,:)>=0)));
savings=[savings; savtauflat];
[SRaggtauflat, Lsharetauflat, Htauflat]=aggstat(y0,optch,eta,tau,r,T1,T2);
%% Solving the economy with pregressive tax rate
theta=0.18;% check HSV
lambda=0.15;
delta=1-lambda;
optch=nan(n,3,length(eta));
r0=0.05;
T10=0.2;
T20=0.2;
rT0=[r0 T10 T20];
rprotax=fsolve(@(rT) GEprotax(y0,rT(1),eta,sigma,beta,theta,delta,rT(2),rT(3),eps,k,nu,etarealized),rT0,options);
r=rprotax(1);
T1=rprotax(2);
T2=rprotax(3);
C0=nan(n,1,length(eta));
C1=nan(n,2,length(eta));
%consumption growth
cg=nan(n,3,length(eta));
% before tax labor income growth
whg=nan(n,3,length(eta));
U=nan(n,2,length(eta));
SR=nan(n,1,length(eta));

                %store labor income before and after taxes
                wh0=nan(n,3,length(eta));
                wh1=nan(n,6,length(eta));
                lc=[0 1 1; 0.6 0.8 0.8; 1 0.2 1; 0 0.4 0.4];
atot=sum(sum(optch(:,1,:)));
govbal1=sum(sum(rev(:,1,:)))-T1*n*length(eta);
govbal2=sum(sum(rev(:,2,:)))-T2*n*length(eta);
for i=1:length(eta)
    tau1=ones(n,1)-delta*( eta(i)*optch(:,2,i) ).^(-theta);
    tau2h=ones(n,1)-delta*( (eta(i)+eps)*optch(:,3,i)).^(-theta);
    tau2l=ones(n,1)-delta*( (eta(i)-eps)*optch(:,3,i)).^(-theta);
                    % compute consumption and store it
                C0(:,:,i)=y0'-optch(:,1,i)+(ones(n,1)-tau1)*eta(i).*optch(:,2,i)+T1;
                % c1 high
                C1(:,1,i)=optch(:,1,i)*(1+r)+(ones(n,1)-tau2h)*(eta(i)+eps).*optch(:,3,i)+T2;
                % c1 low
                C1(:,2,i)=optch(:,1,i)*(1+r)+(ones(n,1)-tau2l)*(eta(i)-eps).*optch(:,3,i)+T2;
                % consumption growth
                cg(:,1,i)=(C1(:,1,i)-C0(:,:,i))./C0(:,:,i);
                cg(:,2,i)=(C1(:,2,i)-C0(:,:,i))./C0(:,:,i);
                cg(:,3,i)= (cg(:,1,i)+cg(:,2,i))/2 ;
                % utility today  
                U(:,1,i)=( (C0(:,:,i)).^(1-sigma) )./(1-sigma) - k* ((optch(:,2,i)).^(1+1/nu) )./(1+1/nu);
                % discounted expected utility tomorrow
                U(:,2,i)=beta* 0.5* [ ( (C1(:,1,i)).^(1-sigma) )./(1-sigma) - k* ((optch(:,3,i)).^(1+1/nu) )./(1+1/nu)]+...
                         beta* 0.5* [ ( (C1(:,2,i)).^(1-sigma) )./(1-sigma) - k* ((optch(:,3,i)).^(1+1/nu) )./(1+1/nu)];
                % Saving rate
                SR(:,1,i)=optch(:,1,i)./y0';
                
                % pretax labor income t=0
                wh0(:,1,i)=optch(:,2,i)*eta(i);
                % after-tax labor income t=0 
                wh0(:,2,i)=optch(:,2,i)*eta(i).*(ones(n,1)-tau1);
                % total after-tax income t=0
                wh0(:,3,i)=wh0(:,2,i)+y0'+T1;
                
                % pretax labor income t=1 eps=+
                wh1(:,1,i)=optch(:,3,i)*(eta(i)+eps);
                % after-tax labor income t=1 eps=+
                wh1(:,2,i)=optch(:,3,i)*(eta(i)+eps).*(ones(n,1)-tau2h);
                % total after-tax income t=1 eps=+
                wh1(:,3,i)=wh1(:,2,i)+optch(:,1,i)*(1+r)+T2;
                % pretax labor income t=1 eps=-
                wh1(:,4,i)=optch(:,3,i)*(eta(i)-eps);
                % after-tax labor income t=1 eps=-
                wh1(:,5,i)=optch(:,3,i)*(eta(i)-eps).*(ones(n,1)-tau2l);
                % total after-tax income t=1 eps=-
                wh1(:,6,i)=wh1(:,5,i)+optch(:,1,i)*(1+r)+T2;
                
                % realized labor income in the second period
                whreal=nan(n,3,length(eta));
                for i=1:length(eta)
                    %before tax
                    whreal(:,1,i)=optch(:,3,i).*etarealized(:,i);
                    %after tax
                    whreal(:,2,i)=whreal(:,1,i)*(1-tau);
                    % total income after tax
                    whreal(:,3,i)=whreal(:,2,i)+T2+optch(:,1,i)*(1+r);
                end
                % before tax labor income growth high shock
                whg(:,1,i)=(wh1(:,1,i)-wh0(:,1,i))./wh0(:,1,i);
                % before tax labor income growth low shock
                whg(:,2,i)=(wh1(:,4,i)-wh0(:,1,i))./wh0(:,1,i);
                % expected before tax labor income 
                whg(:,3,i)=(whg(:,1,i)+whg(:,2,i))/2 ;
                
end               
 for i=1:length(eta)            
                figure(i)
                lci=lc(i,:);
                subplot(3,3,1)
                hold on
                plot(y0,optch(:,1,i),'Color',[lci],'LineStyle','--')
                grid on
                axis tight
                xlabel('y_{0}')
                ylabel('a^{*}')
                legend('\tau = 0',['\tau = ' num2str(tau)],['\lambda = ' num2str(round(lambda,3)),'\theta = ' num2str(theta) ])

                subplot(3,3,2)
                title(['\eta = ' num2str(eta(i)) ],'FontSize', 12)
                hold on
                plot(y0,optch(:,2,i),'Color',[lci],'LineStyle','--')
                grid on
                axis tight
                xlabel('y_{0}')
                ylabel('h_{0}^{*}')
                title({['\eta = ' num2str(eta(i)) ' \epsilon= ' num2str(eps) ' r_{\tau=0}^{*} = ' num2str(round(rtau0,3))];...
                    [' r^{*}_{\tau=0.115}= ' num2str(round(rtauflat,3)) ' T1_{\tau=0.115}^{*}= ' num2str(round(Ttauflat1,3)) ' T2_{\tau=0.15}^{*}= ' num2str(round(Ttauflat2,3))];...
                    [' r^{*}_{\lambda, \theta}^{*}= ' num2str(round(r,3)) ' T_{1}^{*}_{\lambda, \theta}^{*}= ' num2str(round(T1,3)) ' T_{2}^{*}_{\lambda, \theta}^{*}= ' num2str(round(T2,3))]},'FontSize', 12)


                subplot(3,3,3)
                hold on
                plot(y0,optch(:,3,i),'Color',[lci],'LineStyle','--')
                grid on
                axis tight
                xlabel('y_{0}')
                ylabel('h_{1}^{*}')


                subplot(3,3,4)
                hold on
                plot(y0,C0(:,:,i),'Color',[lci],'LineStyle','--')
                grid on
                axis tight
                xlabel('y_{0}')
                ylabel('c_{0}^{*}')

                subplot(3,3,5)
                hold on
                plot(y0,C1(:,1,i),'Color',[lci],'LineStyle','--')
                hold on
                plot(y0,C1(:,2,i),'Color',[0.6 1 0.2],'LineStyle','--')
                grid on
                axis tight
                xlabel('y_{0}')
                ylabel('c_1^{*}')
                 legend ('\tau = 0, \epsilon_{H}','\tau = 0, \epsilon_{L}',...
                    ['\tau = ' num2str(tau) ' \epsilon_{H}'],['\tau = ' num2str(tau) ' \epsilon_{L}'],...
                     ['\lambda = ' num2str(round(lambda,3))  ', \theta = ' num2str(theta) ', \epsilon_{H}'],...
                     ['\lambda = ' num2str(round(lambda,3)) ', \theta = ' num2str(theta) ', \epsilon_{L}'])
                 subplot(3,3,6)
                plot(y0,cg(:,1,i),'Color',[lci],'LineStyle','--')
                hold on
                plot(y0,cg(:,2,i),'Color',[0.6 1 0.2],'LineStyle','--')
                grid on
                axis tight
                xlabel('y_{0}')
                ylabel('g_{c}')
                
                subplot(3,3,7)
                hold on
                plot(y0,SR(:,1,i),'Color',[lci],'LineStyle','--')
                grid on
                axis tight
                xlabel('y_{0}')
                ylabel('a^{*}/y_{0}')
                
                subplot(3,3,8)
                hold on
                plot(y0,whg(:,1,i),'Color',[lci],'LineStyle','--')
                hold on
                plot(y0,whg(:,2,i),'Color',[0.6 1 0.2],'LineStyle','--')
                grid on
                axis tight
                xlabel('y_{0}')
                ylabel('g_{wh}')
                
                subplot(3,3,9)
                hold on
                plot(y0,cg(:,3,i)./whg(:,3,i),'Color',[lci],'LineStyle','--')               
                grid on
                axis tight
                xlabel('y_{0}')
                ylabel('E[g_{c}]/E[g_{wh}]')
                
end


 % welfare function, as function of initial wealth (across productivities)
 W0=sum(U(:,1,:),3);
 W1=sum(U(:,2,:),3);
 figure (length(eta)+1)
 subplot(2,3,1)
 hold on
 plot(y0,permute(U(:,1,:),[1,3,2]))
 legend(['\eta = ' num2str(eta(1)) ', \tau=0'],['\eta = ' num2str(eta(2)) ', \tau=0'],...
        ['\eta = ' num2str(eta(3)) ', \tau=0'],['\eta = ' num2str(eta(4)) ', \tau=0'],...
        ['\eta = ' num2str(eta(1)) ', \tau=' num2str(tau)],['\eta = ' num2str(eta(2)) ', \tau' num2str(tau)],...
        ['\eta = ' num2str(eta(3)) ', \tau' num2str(tau)],['\eta = ' num2str(eta(4)) ', \tau' num2str(tau)],...
    ['\eta = ' num2str(eta(1)) ', \lambda=' num2str(round(lambda,3)) ', \theta=' num2str(theta)],['\eta = ' num2str(eta(2)) ', \lambda=' num2str(delta) ', \theta=' num2str(theta)],...
    ['\eta = ' num2str(eta(3)) ', \lambda=' num2str(round(lambda,3)) ', \theta=' num2str(theta)],['\eta = ' num2str(eta(4)) ', \lambda=' num2str(delta) ', \theta=' num2str(theta)])
 subplot(2,3,2)
 hold on
 plot(y0,permute(U(:,2,:),[1,3,2]))
 subplot(2,3,3)
 hold on
 plot(y0,permute(U(:,1,:),[1,3,2])+permute(U(:,2,:),[1,3,2]))
 subplot(2,3,4)
 hold on
 plot(y0,W0,'Color',[0 0 1.0],'Linestyle',':')
 legend('\tau = 0',['\tau =' num2str(tau)],['\lambda=' num2str(lambda) ', \theta=' num2str(theta)])
 subplot(2,3,5)
 hold on
 plot(y0,W1,'Color',[0 0 1.0],'Linestyle',':')
 subplot(2,3,6)
 hold on
 plot(y0,W0+W1,'Color',[0 0 1.0],'Linestyle',':')

 
 
taupro=[delta;theta]; 
[out]=graph1(y0,optch, C0, C1, SR, taupro, eps, eta, r, T1,T2)
[out]=graph2(y0,wh0,wh1,taupro, eps, eta, r, T1,T2)
[out]=graph3(y0,cg,whg, taupro, eps, eta, r, T1,T2)
                display('--------------------------------------------------------')
                display(['The equilibrium interest rate for lambda= ' num2str(lambda) 'and theta= ' num2str(theta) ' is ' num2str(r)])
                display(['The equilibrium Transfer in t=0 for lambda= ' num2str(lambda) 'and theta= ' num2str(theta) ' is ' num2str(T1)])
                display(['The equilibrium Transfer in t=1 for lambda= ' num2str(lambda) 'and theta= ' num2str(theta) ' is ' num2str(T2)])
                display(['At this interest rate and Transfers net assets are ' num2str(atot)])
                display(['At this interest rate and Transfers governement budget in t=0 is ' num2str(govbal1)])
                display(['At this interest rate and Transfers governement budget in t=1 is ' num2str(govbal2)])
                display('--------------------------------------------------------')
rtaupro=r;
Ttaupro1=T1;
Ttaupro2=T2;
results=[results; rtaupro Ttaupro1 Ttaupro2];
savtaupro=sum(sum(optch(optch(:,1,:)>=0)));
savings=[savings; savtaupro];
[SRaggtaupro, Lsharetaupro, Htaupro]=aggstat(y0,optch,eta,taupro,r,T1,T2);
% progressivity
etah=[];
tauetah=[];
MTR0=[];
figure
for j=1:length(eta)
    etajh=optch(:,2,i).*eta(j);
    tauetajh=ones(n,1)-delta*etajh.^(-theta);
    MTRj=ones(n,1)-(1-theta)*delta*etajh.^(-theta);
    subplot(length(eta)/2,2,j)
    plot(etajh,tauetajh)
    hold on
    plot(etajh,MTRj)
    title(['\eta = ' num2str(eta(j)) '\lambda= ' num2str(lambda) '\theta= ' num2str(theta)])
    ylabel('\tau_{\eta},MTR_{\eta}')
    xlabel('\eta h^{*}')
    axis tight
    grid on
    legend('\tau','MTR')
   %etajh2=optch(:,3,i).*etarealized(j) 
    etah=[etah;etajh];
    tauetah=[tauetah;tauetajh];
    MTR0=[MTR0;MTRj];
end

figure
plot(etah,tauetah)
hold on
plot(etah,MTR0)
title(['\tau, t=0 ' ' \lambda= ' num2str(lambda) ' \theta= ' num2str(theta)])
grid on
axis tight
ylabel('\tau,MTR')
xlabel('\eta h^{*}')
legend('\tau','MTR')

aggstatistics=[SRaggtau0 Lsharetau0 Htau0;
    SRaggtauflat Lsharetauflat Htauflat;
    SRaggtaupro Lsharetaupro Htaupro];

