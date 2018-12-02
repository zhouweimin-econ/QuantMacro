function Fit_model_graph(X1, Y1, X2, Y2)
%CREATEFIGURE(X1, Y1, S1, C1, X2, Y2)
%  X1:  Debt to DP ratio
%  Y1:  Default decision
%  X2:  Debt to GDP ratio grid
%  Y2:  estimated default probability

S1=20;

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create scatter
scatter(X1,Y1,S1,'red','DisplayName','Default decision',...
    'MarkerFaceColor',[1 0 0],...
    'MarkerEdgeColor',[1 0 0]);

% Create plot
plot(X2,Y2,'DisplayName','Estimated default probability','LineWidth',4,...
    'Color',[0 0 0]);

% Create xlabel
xlabel('Debt to GDP ratio');

% Create ylabel
ylabel('Default probability');

% Set the remaining axes properties
set(axes1,'FontSize',20);
% Create legend
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.157816438149007 0.668727906009565 0.530357142857142 0.125],...
    'FontSize',20);

