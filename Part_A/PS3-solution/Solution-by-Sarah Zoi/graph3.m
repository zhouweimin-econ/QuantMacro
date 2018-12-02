function [out]=graph3(y0,cg,whg, tau, eps, eta, rstar, Tstar1, Tstar2)

if nargin>=7
    %figure 2
    out=figure
    set(0,'DefaultAxesColorOrder',[ 0 0 1.0;0 0 0;0.8 0 0;0 0.8 0.2],...
                    'defaultLineLineWidth',1.3,'DefaultAxesLineStyleOrder','-|--|:') %; 0 1 0.8;0.6 0.6 0.6; 1 0.4 0.4 ; 0.6 1 0 
     
                subplot(2,3,1)
                plot(y0,permute(cg(:,3,:),[1,3,2]))
                hold on
                plot(y0,permute(cg(:,1,:),[1,3,2]))
                hold on
                plot(y0,permute(cg(:,2,:),[1,3,2]))
                grid on
                axis tight
                xlabel('y_{0}')
                ylabel('g_{c}')
                 legend(['\eta = ' num2str(eta(1)), 'E[g]'],['\eta = ' num2str(eta(2)), 'E[g]'],...
                    ['\eta = ' num2str(eta(3)), 'E[g]'],['\eta = ' num2str(eta(4)), 'E[g]' ],...
                    ['\eta = ' num2str(eta(1)), ',  \epsilon_{H}'],['\eta = ' num2str(eta(2)), ', \epsilon_{H}'],...
                    ['\eta = ' num2str(eta(3)), ', \epsilon_{H}'],['\eta = ' num2str(eta(4)), ', \epsilon_{H}' ],...
                    ['\eta = ' num2str(eta(1)), ', \epsilon_{L}'],['\eta = ' num2str(eta(2)), ', \epsilon_{L}'],...
                    ['\eta = ' num2str(eta(3)), ', \epsilon_{L}'],['\eta = ' num2str(eta(4)), ', \epsilon_{L}'] )
                
                subplot(2,3,2)
                plot(y0,permute(whg(:,3,:),[1,3,2]))
                hold on
                plot(y0,permute(whg(:,1,:),[1,3,2]))
                hold on
                plot(y0,permute(whg(:,2,:),[1,3,2]))
                grid on
                axis tight
                xlabel('y_{0}')
                ylabel('g_{wh}')
%                  legend(['\eta = ' num2str(eta(1)), 'E[g_{wh}]'],['\eta = ' num2str(eta(2)), 'E[g_{wh}]'],...
%                     ['\eta = ' num2str(eta(3)), 'E[g_{wh}]'],['\eta = ' num2str(eta(4)), 'E[g_{wh}]' ],...
%                     ['\eta = ' num2str(eta(1)), ',  \epsilon_{H}'],['\eta = ' num2str(eta(2)), ', \epsilon_{H}'],...
%                     ['\eta = ' num2str(eta(3)), ', \epsilon_{H}'],['\eta = ' num2str(eta(4)), ', \epsilon_{H}' ],...
%                     ['\eta = ' num2str(eta(1)), ', \epsilon_{L}'],['\eta = ' num2str(eta(2)), ', \epsilon_{L}'],...
%                     ['\eta = ' num2str(eta(3)), ', \epsilon_{L}'],['\eta = ' num2str(eta(4)), ', \epsilon_{L}'] )
   if nargin==9 && length(tau)==1
                title([' \tau= ' num2str(tau) ' \epsilon= ' num2str(eps)...
                    ' r^{*}= ' num2str(round(rstar,3)) ' T_{1}^{*}= ' num2str(round(Tstar1,3)) ' T_{2}^{*}= ' num2str(round(Tstar2,3))],'FontSize', 12)
   elseif nargin==7 && length(tau)==1
       title([' \tau= ' num2str(tau) ' \epsilon= ' num2str(eps)...
                    ' r^{*}= ' num2str(round(rstar,3))],'FontSize', 12)
   elseif nargin==9 && length(tau)==2
       title([' \lambda= ' num2str(tau(1)) ' \theta= ' num2str(tau(2)) ' \epsilon= ' num2str(eps)...
                    ' r^{*}= ' num2str(round(rstar,3)) ' T_{1}^{*}= ' num2str(round(Tstar1,3)) ' T_{2}^{*}= ' num2str(round(Tstar2,3))],'FontSize', 12)
   end
                subplot(2,3,3)
                plot(y0,permute(cg(:,3,:),[1,3,2])./permute(whg(:,3,:),[1,3,2]))
                grid on
                axis tight
                xlabel('y_{0}')
                ylabel('E[g_{c}]/E[g_{wh}]')
                
                
                subplot(2,3,4)
                plot(y0, ( permute(cg(:,1,:),[1,3,2])./permute(whg(:,1,:),[1,3,2]))./( permute(cg(:,3,:),[1,3,2])./permute(whg(:,3,:),[1,3,2])),'Linestyle','--' )
                hold on
                plot(y0, ( permute(cg(:,2,:),[1,3,2])./permute(whg(:,2,:),[1,3,2]))./( permute(cg(:,3,:),[1,3,2])./permute(whg(:,3,:),[1,3,2])),'Linestyle',':' )
                grid on
                axis tight
                xlabel('y_{0}')
                ylabel('{ g_{c}/g_{wh} } / { E [g_{c}] /E[ g_{wh}]}')
                 legend(['\eta = ' num2str(eta(1)), ',  \epsilon_{H}'],['\eta = ' num2str(eta(2)), ', \epsilon_{H}'],...
                    ['\eta = ' num2str(eta(3)), ', \epsilon_{H}'],['\eta = ' num2str(eta(4)), ', \epsilon_{H}' ],...
                    ['\eta = ' num2str(eta(1)), ', \epsilon_{L}'],['\eta = ' num2str(eta(2)), ', \epsilon_{L}'],...
                    ['\eta = ' num2str(eta(3)), ', \epsilon_{L}'],['\eta = ' num2str(eta(4)), ', \epsilon_{L}'] )
                



end