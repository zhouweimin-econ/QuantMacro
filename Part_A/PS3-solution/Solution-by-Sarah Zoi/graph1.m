function [out]=graph1(y0,optch, C0, C1, SR, tau, eps, eta, rstar, Tstar1, Tstar2)

if nargin>=9
    %figure 1 
out=figure
                set(0,'DefaultAxesColorOrder',[ 0 0 1.0;0 0 0;0.8 0 0;0 0.8 0.2],...
                    'defaultLineLineWidth',1.3,'DefaultAxesLineStyleOrder','-|--') %; 0 1 0.8;0.6 0.6 0.6; 1 0.4 0.4 ; 0.6 1 0 
                subplot(2,3,1)
                plot(y0,permute(optch(:,1,:),[1,3,2]))
                grid on
                axis tight
                xlabel('y_{0}')
                ylabel('a^{*}')
                legend(['\eta = ' num2str(eta(1))],['\eta = ' num2str(eta(2))],...
                    ['\eta = ' num2str(eta(3))],['\eta = ' num2str(eta(4))])
                subplot(2,3,2)
                plot(y0,permute(optch(:,2,:),[1,3,2]))
                grid on
                axis tight
                xlabel('y_{0}')
                ylabel('h_{0}^{*}')
               %legend(['\eta = ' num2str(eta(1))],['\eta = ' num2str(eta(2))],...
                   % ['\eta = ' num2str(eta(3))],['\eta = ' num2str(eta(4))])
                if nargin==9 && length(tau)==1
                title([' \tau= ' num2str(tau) ' \epsilon= ' num2str(eps)...
                    ' r^{*}= ' num2str(round(rstar,3))],'FontSize', 12)
                elseif nargin==11 && length(tau)==1
                    title([' \tau= ' num2str(tau) ' \epsilon= ' num2str(eps)...
                    ' r^{*}= ' num2str(round(rstar,3)) ' T_{1}^{*}= ' num2str(round(Tstar1,3)) ' T_{2}^{*}= ' num2str(round(Tstar2,3))],'FontSize', 12)
                elseif nargin==11 && length(tau)~=1
                   title([' \lambda= ' num2str(tau(1)) ' \theta= ' num2str(tau(2)) ' \epsilon= ' num2str(eps)...
                    ' r^{*}= ' num2str(round(rstar,3)) ' T_{1}^{*}= ' num2str(round(Tstar1,3)) ' T_{2}^{*}= ' num2str(round(Tstar2,3))],'FontSize', 12) 
                end
                subplot(2,3,3)
                plot(y0,permute(optch(:,3,:),[1,3,2]))
                grid on
                axis tight
                xlabel('y_{0}')
                ylabel('h_{1}^{*}')
%                 legend(['\eta = ' num2str(eta(1))],['\eta = ' num2str(eta(2))],...
%                     ['\eta = ' num2str(eta(3))],['\eta = ' num2str(eta(4))])
                subplot(2,3,4)
                plot(y0,permute(C0(:,1,:),[1,3,2]))
                grid on
                axis tight
                xlabel('y_{0}')
                ylabel('c_{0}^{*}')
%                 legend(['\eta = ' num2str(eta(1))],['\eta = ' num2str(eta(2))],...
%                     ['\eta = ' num2str(eta(3))],['\eta = ' num2str(eta(4))])
                
                subplot(2,3,5)
                plot(y0,permute(C1(:,1,:),[1,3,2]))
                hold on 
                plot(y0,permute(C1(:,2,:),[1,3,2]))
                grid on
                axis tight
                xlabel('y_{0}')
                ylabel('c_{0}^{*}')
                legend(['\eta = ' num2str(eta(1)), ',  \epsilon_{H}'],['\eta = ' num2str(eta(2)), ', \epsilon_{H}'],...
                    ['\eta = ' num2str(eta(3)), ', \epsilon_{H}'],['\eta = ' num2str(eta(4)), ', \epsilon_{H}' ],...
                    ['\eta = ' num2str(eta(1)), ', \epsilon_{L}'],['\eta = ' num2str(eta(2)), ', \epsilon_{L}'],...
                    ['\eta = ' num2str(eta(3)), ', \epsilon_{L}'],['\eta = ' num2str(eta(4)), ', \epsilon_{L}'] )
                subplot(2,3,6)
                hold on
                plot(y0,permute(SR(:,1,:),[1,3,2]))
                grid on
                axis tight
                xlabel('y_{0}')
                ylabel('a^{*}/y_{0}')
%                 legend(['\eta = ' num2str(eta(1))],['\eta = ' num2str(eta(2))],...
%                     ['\eta = ' num2str(eta(3))],['\eta = ' num2str(eta(4))])


end              
                
end