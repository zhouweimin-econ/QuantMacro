function Default_graph(yvector1, Y1)

%  YVECTOR1:  default events data
%  Y1:  Risk premia data

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1,...
    'Position',[0.13 0.11 0.761126025354213 0.815]);
hold(axes1,'on');

% Activate the left side of the axes
yyaxis(axes1,'left');
% Create area
area(yvector1,'DisplayName','Default','FaceColor',[1 0 0],...
    'EdgeColor',[1 0 0]);

% Create ylabel
ylabel('DEfault index','FontSize',20);

% Set the remaining axes properties
set(axes1,'YColor',[0 0 0]);
% Activate the right side of the axes
yyaxis(axes1,'right');
% Create plot
plot(Y1,'DisplayName','Risk premia','Color',[0 0 0]);

% Create ylabel
ylabel('Risk premia','FontSize',20);

% Set the remaining axes properties
set(axes1,'YColor',[0.85 0.325 0.098]);
% Create xlabel
xlabel('Time t','FontSize',20);

% Uncomment the following line to preserve the X-limits of the axes
% xlim(axes1,[0 500]);
box(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontSize',14);
% Create legend
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.655853840417598 0.175062972292191 0.109619686800895 0.132241813602015],...
    'FontSize',20);

