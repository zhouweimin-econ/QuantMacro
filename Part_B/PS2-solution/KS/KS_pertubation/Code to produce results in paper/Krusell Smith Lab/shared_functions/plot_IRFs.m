close all

%% Plotting

figurename=['IRF_Y_' casename  '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 1000 1000])
plot(1:mpar.maxlag-1,IRF_Y,'LineWidth',3.5)
ylabel('Percent','Interpreter','none','FontName','arial','FontSize',40)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',40)
hold on
plot(0:mpar.maxlag-1,zeros(1,mpar.maxlag),'--black')
set(gca, 'FontName','arial','FontSize',40);
set(gca,'Position',[0.2 0.15 0.75 0.8])
printpdf(gcf,['../latex/' figurename])

figurename=['IRF_C_' casename  '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 1000 1000])
plot(1:mpar.maxlag-1,IRF_C,'LineWidth',3.5)
ylabel('Percent','Interpreter','none','FontName','arial','FontSize',40)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',40)
hold on
plot(0:mpar.maxlag-1,zeros(1,mpar.maxlag),'--black')
set(gca, 'FontName','arial','FontSize',40); set(gca,'Position',[0.2 0.15 0.75 0.8]) 
printpdf(gcf,['../latex/' figurename])

figurename=['IRF_I_' casename  '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 1000 1000])
plot(1:mpar.maxlag-1,IRF_I,'LineWidth',3.5)
ylabel('Percent','Interpreter','none','FontName','arial','FontSize',40)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',40)
hold on
plot(0:mpar.maxlag-1,zeros(1,mpar.maxlag),'--black')
set(gca, 'FontName','arial','FontSize',40); set(gca,'Position',[0.2 0.15 0.75 0.8]) 
printpdf(gcf,['../latex/' figurename])

figurename=['IRF_K_' casename  '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 1000 1000])
plot(1:mpar.maxlag,IRF_K,'LineWidth',3.5)
ylabel('Percent','Interpreter','none','FontName','arial','FontSize',40)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',40)
hold on
plot(0:mpar.maxlag-1,zeros(1,mpar.maxlag),'--black')
set(gca, 'FontName','arial','FontSize',40); set(gca,'Position',[0.2 0.15 0.75 0.8]) 
printpdf(gcf,['../latex/' figurename])


figurename=['IRF_H_' casename  '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 1000 1000])
plot(1:mpar.maxlag-1,IRF_H,'LineWidth',3.5)
ylabel('Percent','Interpreter','none','FontName','arial','FontSize',40)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',40)
hold on
plot(0:mpar.maxlag-1,zeros(1,mpar.maxlag),'--black')
set(gca, 'FontName','arial','FontSize',40); set(gca,'Position',[0.2 0.15 0.75 0.8]) 
printpdf(gcf,['../latex/' figurename])

figurename=['IRF_S_' casename  '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 1000 1000])
plot(1:mpar.maxlag-1,IRF_S,'LineWidth',3.5)
ylabel('Percent','Interpreter','none','FontName','arial','FontSize',40)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',40)
hold on
plot(0:mpar.maxlag-1,zeros(1,mpar.maxlag),'--black')
set(gca, 'FontName','arial','FontSize',40); set(gca,'Position',[0.2 0.15 0.75 0.8]) 
printpdf(gcf,['../latex/' figurename])


figurename=['IRF_Q_' casename  '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 1000 1000])
plot(1:mpar.maxlag-1,IRF_Q,'LineWidth',3.5)
ylabel('Basis Points','Interpreter','none','FontName','arial','FontSize',40)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',40)
hold on
plot(0:mpar.maxlag-1,zeros(1,mpar.maxlag),'--black')
set(gca, 'FontName','arial','FontSize',40); set(gca,'Position',[0.2 0.15 0.75 0.8]) 
printpdf(gcf,['../latex/' figurename])

figurename=['IRF_R_' casename  '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 1000 1000])
plot(1:mpar.maxlag-1,IRF_R,'LineWidth',3.5)
ylabel('Basis Points','Interpreter','none','FontName','arial','FontSize',40)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',40)
hold on
plot(0:mpar.maxlag-1,zeros(1,mpar.maxlag),'--black')
set(gca, 'FontName','arial','FontSize',40); set(gca,'Position',[0.2 0.15 0.75 0.8]) 
printpdf(gcf,['../latex/' figurename])


figurename=['IRF_N_' casename  '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 1000 1000])
plot(1:mpar.maxlag-1,IRF_N,'LineWidth',3.5)
ylabel('Percent','Interpreter','none','FontName','arial','FontSize',40)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',40)
hold on
plot(0:mpar.maxlag-1,zeros(1,mpar.maxlag),'--black')
set(gca, 'FontName','arial','FontSize',40); set(gca,'Position',[0.2 0.15 0.75 0.8]) 
printpdf(gcf,['../latex/' figurename])

close all