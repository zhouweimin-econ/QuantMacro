function sim_and_model_graph(YMatrix1)
%CREATEFIGURE(YMATRIX1)
%  YMATRIX1:  matrix with default prob data. 
% First the estimated probability



% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1,...
    'Position',[0.086996336996337 0.127897091722595 0.883699633699634 0.797102908277405]);
hold(axes1,'on');

% Create multiple lines using matrix input to plot
plot1 = plot(YMatrix1,'LineWidth',4);
set(plot1(1),'DisplayName','Estimated');
set(plot1(2),'DisplayName','Model-based','LineStyle','-.');

% Create xlabel
xlabel('Time');

% Create ylabel
ylabel('Default probability');

box(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontSize',20);
% Create legend
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.798878125902623 0.811764242320324 0.178010794604005 0.11744966442953],...
    'FontSize',20);

