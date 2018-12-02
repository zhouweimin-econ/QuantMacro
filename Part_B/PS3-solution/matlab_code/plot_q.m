function plot_q(X1, YMatrix1)

%  X1:  vector with the grid of debt issuance 
%  YMATRIX1:  matrix of debt prices for each level of y



% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create multiple lines using matrix input to plot
plot1 = plot(X1,YMatrix1,'LineWidth',6,'Parent',axes1);
set(plot1(1),'DisplayName','Low y');
set(plot1(2),'DisplayName','Medium y');
set(plot1(3),'DisplayName','High y');

% Create xlabel
xlabel('Debt issued b''','FontSize',20);

% Create ylabel
ylabel('Bond price q(b'',y)','FontSize',20);

box(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontSize',14,'XGrid','on','YGrid','on');
% Create legend
legend(axes1,'show');

