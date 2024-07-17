%% Hierarchical clustering
%  Version : v1.0
%  Author  : E. Ogier
%  Release : 30th april 2016
%
%  USAGE :
%
%  [FIGURE, TREE, CLUSTERS, ROOTS] = Hierarchical_clustering(MATRIX,METHOD,METRIC,PARAMETER,VALUE,COLORMAP);
%          [TREE, CLUSTERS, ROOTS] = Hierarchical_clustering(MATRIX,METHOD,METRIC,PARAMETER,VALUE,COLORMAP);
%
%  INPUT ARGUMENTS :
%  - MATRIX    : matrix of items values (matrix or vector, according to item type, scalar or vector)
%  - METHOD    : dissimilarity evaluation method : - 'WPGMA' (weighted pair group method with averaging)
%                                                  - 'UPGMA' (unweighted pair group method with averaging)
%  - METRIC    : metric anonymous function
%  - PARAMETER : clustering parameter :            - 'Number'
%                                                  - 'Limit'
%  - VALUE     : parameter value :                 - Number of clusters if PARAMETER is equal to 'Number'
%                                                  - Clustering dissimilarity limit if PARAMETER is equal to 'Limit'                                      
%  - COLORMAP  : colormap name among :              'parula', 'jet', 'hsv', 'hot', 'cool', 'spring', 'summer', 'autumn', 'winter', 'gray', 'bone', 'copper' and 'pink'
%
%  OUTPUT ARGUMENTS :
%
%  - FIGURE    : figure including dendrogram and clusters (optional, not defined if no graphic representation is required by user)
%  - TREE      : tree of hierarchical clusters
%  - CLUSTERS  : clusters structure containing the items of each cluster (dimension defined by user through PARAMETER argument)
%  - ROOTS     : list of dendrogram roots ranked according to the order of its branches
%
%  EXAMPLE 1 :
%
%  % Vector
%  % Method : UPGMA
%  % Metric : norm of difference
% 
%  % Setting
%  Matrix   = rand(50,1);
%  Method   = 'WPGMA';
%  Metric   = @(X,Y)abs(X-Y);
%  Limit    = 0.2;
%  Colormap = 'Cool';
%  
%  % Hierarchical clustering & dendrogram
%  [Figure, Tree, Clusters, Roots] = Hierarchical_clustering(Matrix,Method,Metric,'Limit',Limit,Colormap);
%
%  EXAMPLE 2 :
%
%  % Matrix of tridimensional coordinates
%  %  Method : UPGMA
%  %  Metric : norm of difference
%  
%  % Setting
%  Matrix   = rand(50,3);
%  Method   = 'UPGMA';
%  Metric   = @(X,Y)norm(X-Y);
%  Number   = 5;
%  Colormap = 'spring';
%  
%  % Abnormal point
%  Matrix(5,2) = 1.5;
%  
%  % Hierarchical clustering & dendrogram
%  [Figure, Tree, Clusters, Roots] = Hierarchical_clustering(Matrix,Method,Metric,'Number',Number,Colormap);
%
%  EXAMPLE 3 :
%
%  % Matrix of bidimensional coordinates
%  %  Method : UPGMA
%  %  Metric : norm of difference
%  
%  % Centroids coordinates
%  Centroids = [ 0  0;...
%                0 -1;...                    
%               -1 -1;...
%               -1  0;...
%               -1 +1;...
%                0  +1;...
%               +1 +1;...
%               +1  0;...
%               +1 -1];          
%  N = 50;                
%  M = size(Centroids,1);
%  Centroids = reshape(Centroids,1,2*M);                
%  Centroids = repmat(Centroids,N,1);
%  Centroids = reshape(Centroids,[],2);
%  
%  % Coordinates in presence of noises in a polar reference
%  r0 = 0.5;
%  Angle       = 2*pi*rand(N*M,1);
%  Radius      = 2*rand(N*M,1)-1;
%  Coordinates = r0*[Radius.*cos(Angle) Radius.*sin(Angle)];
%  
%  % Matrix
%  Matrix = Centroids+Coordinates;
%  
%  % Figure
%  Figure = figure('Color','w');
%  
%  % Full screen
%  drawnow;
%  warning('off','all');
%  jFrame = get(Figure,'JavaFrame');
%  jFrame.setMaximized(true);
%  warning('on','all');
%  pause(0.1);
%  
%  % Plot
%  title('Hierarchical clustering of bidimensional coordinates points',...
%        'Fontname',   'Times',...
%        'Fontsize',   12,...
%        'Fontweight', 'Demi',...
%        'FontAngle',  'Normal');
%  xlabel('Abscissa');
%  ylabel('Ordinate');
%  set(gca,'Xlim',1.5*[-1 +1],'Ylim',1.5*[-1 +1]);
%  box('on');
%  hold('on');
%  Points = [];
%  for n = 1:size(Matrix,1)
%      switch n
%          case 1
%              Points(1) = plot(Matrix(n,1),Matrix(n,2),'k.');          
%              legend(Points,'Points','Location','NorthEastOutside');
%          otherwise
%              plot(Matrix(n,1),Matrix(n,2),'k.');
%      end
%      drawnow();
%  end
%  
%  % Clustering setting
%  Method   = 'UPGMA';
%  Metric   = @(X,Y)norm(X-Y);
%  Number   = M;
%  Colormap = 'spring';
%  
%  % Hierarchical clustering
%  [Tree, Clusters, Roots] = Hierarchical_clustering(Matrix,Method,Metric,'Number',Number,Colormap);
%  
%  % Clusters display
%  Colors = colormap('hsv');
%  C = numel(Clusters);
%  for c = 1:C    
%      for i = 1:numel(Clusters{c})
%          Item = Clusters{c}(i);
%          Index = floor(c/C*(64-1)+1);
%          Color = Colors(Index,:);                
%          switch i
%              case 1                
%                  Points(c+1) = plot(Matrix(Item,1),Matrix(Item,2),'o','Color',Color);                 %#ok<SAGROW>
%                  Legend = arrayfun(@(c)sprintf('Cluster #%u',c),1:c,'UniformOutput',false);
%                  legend(Points,['Points', Legend],'Location','NorthEastOutside');
%              otherwise
%                  plot(Matrix(Item,1),Matrix(Item,2),'o','Color',Color);
%          end
%          drawnow();
%      end
%  end 
function varargout = hierarchicalcluster(Matrix,Method,Metric,Parameter,Value,Colormap)
% Data control
Matrix = Data_control(Matrix,Method,Metric,Parameter,Value,Colormap);
if nargout == 4    
    % Creation of the bars corresponding to values
    [Figure, Axes, Bars] = Creation_bars(Matrix,Colormap);
    
end
% Linkage
Tree = Linkage(Matrix,Method,Metric);
% Analysis of tree
[Roots, Clusters, Number] = Analysis(Tree,Parameter,Value);
if nargout == 4
    
    % Creation of the dendrogram
    Axes = Creation_dendrogram(Axes,Bars,Tree,Roots,Clusters,Matrix,Method,Metric,Number,Colormap);
    
    % X axes link
    linkaxes(Axes,'x');
    
end
% Output arguments
switch nargout
    case 3
        varargout{1} = Tree;
        varargout{2} = Clusters;
        varargout{3} = Roots;
    case 4
        varargout{1} = Figure;
        varargout{2} = Tree;
        varargout{3} = Clusters;
        varargout{4} = Roots;
end
end
% Data control
function Matrix = Data_control(Matrix,Method,Metric,Parameter,Value,Colormap)
% Matrix
[n,m] = size(Matrix);
% if n < m
%     Matrix = Matrix';
% end
% Method
switch upper(Method)
    case {'WPGMA','UPGMA'}
    otherwise, error('Clustering method must be ''WPGMA'' or ''UPGMA''.');
end
% Metric
switch class(Metric)
    case 'function_handle'
        switch nargin(Metric)
            case 2
            otherwise
                error('Metric must include two arguments (detected arguments: %u).',nargin(f));
        end
    otherwise, error('Metric must be an anonymous function (detected class: ''%s'').',class(Metric));
end
% Parameter
switch lower(Parameter)
    case {'number','limit'}
    otherwise
        error('Clustering parameter must be ''Number'' or ''Limit''.');
end
% Parameter value
if ~isnumeric(Value)
    error('Parameter ''%s'' is not numeric.',Parameter);
end
% Colormap
Colormaps = {'parula', 'jet', 'hsv', 'hot', 'cool', 'spring', 'summer', 'autumn', 'winter', 'gray', 'bone', 'copper','pink'};
if ~ismember(lower(Colormap),Colormaps)
    error('Colormap must be a string among %s ''%s'' and ''%s''.',cell2mat(arrayfun(@(s)sprintf('''%s'', ',Colormaps{s}),1:numel(Colormaps)-2,'UniformOutput',false)),Colormaps{end-1},Colormaps{end});
end
end
% Creation of the bars corresponding to values
function [Figure, Axes, Bars] = Creation_bars(Matrix,Colormap)
% Figure
Figure = figure('Color','w','Nextplot','Add');
% Full screen
drawnow;
warning('off','all');
jFrame = get(Figure,'JavaFrame');
jFrame.setMaximized(true);
warning('on','all');
pause(0.1);
% Axes
Axes(1) = subplot(4,1,4);
set(Axes(1),'TickLength',[0 0]);
hold(Axes(1),'on');
% Colormap
Colors = colormap(Colormap);
% Preallocztion
[N, M] = size(Matrix);
Bars = zeros(1,N);
% Bars    
switch M
    
    case 1
        
        % Extrema
        Min = min(Matrix);
        Max = max(Matrix);
        
        % Scalar colored bar
        for n = 1:N
            Index = floor((Matrix(n)-Min)/(Max-Min)*(64-1)+1);
            Color = Colors(Index,:);
            Bars(n) = bar(n,Matrix(n),'Facecolor',Color);
        end
        
        % Axes position
        P = get(Axes(1),'Position');
        
        % Colorbar        
        CB = colorbar('Position',[P(1)+P(3)+1e-2 P(2) 2e-2 P(4)]);
        set(CB,'Ytick',[0 1]+1,'YTickLabel',{sprintf('%.3f',Min) sprintf('%.3f',Max)} );
    
    otherwise 
        
        % Stacked monocolor multiple bars
        Bars = bar(Matrix,'stacked'); 
        
      	% Axes position
        P = get(Axes(1),'Position');
        
        % Legend
        Text = arrayfun(@(n)sprintf('Component #%u',n),1:M,'UniformOutput',false);        
        Legend = legend(Text,'Orientation','Horizontal');        
        Pl = get(Legend,'Position');
        set(Legend,'Position',[(1-Pl(3))/2 (P(2)-Pl(4))/2 Pl(3) Pl(4)]);
         
end
% Axes setting
set(Axes(1),...
    'Xlim',     [0 N+1],...
    'Fontsize', 9);
ylabel(Axes(1),'Value');
set(Axes(1),'Ygrid','on','YMinorGrid','on');
end
% Linkage
function Tree = Linkage(Matrix,Method,Metric)
% Number of elements in input vector
N = size(Matrix,1);
% Dissimilarity triangular matrix (lower triangle)
Md = inf(N,N);
for m = 1:N
    for n = m+1:N
        Md(n,m) = Metric(Matrix(n,:),Matrix(m,:));
    end
end
% Preallocations
Indices  = 1:N;
Weights  = ones(1,N);
Tree = zeros(N-1,3);
% Linkage
for t = 1:N-1
    
    % Indices of the minimal dissimilarity
    [n,m] = find(Md == min(min(Md)),1,'first');
    
    % New cluster corresponding to the minimum dissimilarity
    Tree(t,:) = [Indices(sort([n m])), Md(n,m)];
    
    % Preallocation
    I = numel(Indices);
    D = zeros(I,1);
    
    % Distance vector wrt new cluster
    for i = 1:I
        
        % Dissimilarity to the first member in lower triangle
        if i < n,  d1 = Md(n,i);
        else       d1 = Md(i,n);
        end
        
        % Dissimilarity to the second member in lower triangle
        if i <= m, d2 = Md(m,i);
        else       d2 = Md(i,m);
        end
        
        % Distance calculation
        switch upper(Method)
            
            % Weighted pair group method with averaging
            case 'WPGMA'
                D(i) = mean([d1 d2]);
                
            % Unweighted pair group method with averaging    
            case 'UPGMA'
                D(i) = (Weights(n)*d1 + Weights(m)*d2) / (Weights(n) + Weights(m)); 
                
        end
        
    end
    
    % Indices including the new cluster and excluding its members
    Indices(n) = [];  Indices(m) = [];    
    Indices = [N+t Indices]; %#ok<AGROW>
    
    % Weight vector    
    w = Weights(n) + Weights(m);
    Weights(n) = [];  Weights(m) = []; 
    Weights = [w Weights] ; %#ok<AGROW>
    
    % Dissimilarity matrix including the new cluster and excluding its members
    D(n) = [];        D(m) = [];         % Dissimilarity vector without new cluster members
    Md([n m],:) = []; Md(:,[n m]) = [];  % Dissimilarity matrix without new cluster members
    Md = [inf(1,I-1); D,Md]; %#ok<AGROW> % Dissimilarity matrix with new cluster
    
end
end
% Analysis of tree
function [Roots, Clusters, Number] = Analysis(Tree,Parameter,Value)
% Empty roots and nodes vectors
Roots = [];
Nodes = [];
% Vector length
N = max(max(Tree(:,1:2)))/2+1;
% Clustering parameter
switch lower(Parameter)
    
    % Number of clusters
    case 'number'
        Number = Value;
        Limit = Tree(end-Number+1+1,3);
        if Limit == 0
            return
        end
        
    % Dissimilarity limit
    case 'limit'
        Limit = Value;
        Number = N-find(Tree(:,3)>=Limit,1,'first')+1;
        
end
% Clusters
Clusters = cell(Number,1);
% Cluster number
Cluster = 0;
% Exploration of the node
ExplorationDown(Tree(end,1),Tree(end,3));
ExplorationDown(Tree(end,2),Tree(end,3));
    % Exploration of the subnodes of a node
    function ExplorationDown(Node,Dissimilarity)
        
        if ismember(Node,Nodes) || ismember(Node,Roots)
            return
        end
        
        % Adding of the current node in nodes list
        Nodes = [Nodes Node];        
        
        if Node <= N
            
            % Root
            Roots = [Roots Node];
           
            % Root whose distance is higher than limit
            [n,~]=find(Tree(:,1:2)==Node);            
            if Tree(n,3) >= Limit
                Cluster = Cluster+1;
            end
            
            % Cluster
            Clusters{Cluster} = [Clusters{Cluster} Node];
            
        else
            
            % Nodes
            Node = Node-N;
            
            % Subnodes
            N1 = Tree(Node,1);
            N2 = Tree(Node,2);
            
            % Cluster index increment
            if Tree(Node,3) < Limit && Dissimilarity >= Limit
                Cluster = Cluster+1;
            end
            
            % Dissimilarity of current node
            Dissimilarity = Tree(Node,3);
            
            if N1 <= N && N2 <= N
                
                % Roots subnodes
                if N1 < N2
                    ExplorationDown(N1,Dissimilarity);
                    ExplorationDown(N2,Dissimilarity);
                else
                    ExplorationDown(N2,Dissimilarity);
                    ExplorationDown(N1,Dissimilarity);
                end
                
            else
                
                % Exploration of the subnodes
                ExplorationDown(N1,Dissimilarity);
                ExplorationDown(N2,Dissimilarity);
                
            end
            
        end
        
    end
end
% Creation of the dendrogram
function Axes = Creation_dendrogram(Axes,Bars,Tree,Roots,Clusters,Matrix,Method,Metric,Number,Colormap)
% Axes
Axes(2) = subplot(4,1,1:3);
set(Axes(2),'Ygrid','on','TickLength',[0 0],'YMinorGrid','on');
title('Hierarchical clustering & dendrogram',...
      'Fontname',   'Times',...
      'Fontsize',   12,...
      'Fontweight', 'Demi',...
      'FontAngle',  'Normal');
ylabel(Axes(2),'Dissimilarity');
hold(Axes(2),'on');
% Number of elements
N = max(max(Tree(:,1:2)))/2+1;
% Matrix width
M = size(Matrix,2);
% Number of roots
R = numel(Roots);
% Number of Tree or nodes
T = size(Tree,1);
% Roots coordinates
AbscissaRoots  = 1:R;
OrdinatesRoots = zeros(1,R);
% Nodes coordinates
AbscissaNodes  = zeros(1,N);
OrdinatesNodes = zeros(1,N);
% Axes limits
Xlim = [0 N+1];
k = 5;
Ylim = [0 ceil(max(Tree(:,3))*k)/k];
set(Axes(2),...
    'Xlim',     Xlim,...
    'Ylim',     Ylim,...
    'Fontsize', 9);
% Indices
Indices = 1:N;
% Colormap
Colors = colormap(Colormap);
% Extrema
Min = min(arrayfun(@(n)norm(Matrix(n,:)),1:size(Matrix,1)));
Max = max(arrayfun(@(n)norm(Matrix(n,:)),1:size(Matrix,1)));
% Y axes limit
Ylim = get(Axes(2),'Ylim');
% Clustering limit
Limit = mean([Tree(end-Number+2,3) Tree(end-Number+1,3)]);
plot([0.5 N+0.5],Limit*[1 1],'r--');
% Text lines
Line1 = sprintf('Method  : %s',  Method);
Line2 = sprintf('Metric    : %s',char(Metric));
Line3 = sprintf('Clusters : %u',Number);
% Text
Text = ...
    text(0.5,Ylim(2),           {Line1,Line2,Line3},...
         'HorizontalAlignment', 'Left',...
         'VerticalAlignment',   'Top',...
         'Fontsize',            9,...
         'FontAngle',           'Normal',...
         'Color',               'b',...
         'EdgeColor',           'b',...
         'BackGroundColor',     'w');
 
% Figure
Figure = get(Axes(2),'Parent');
% Spinner
Spinner(Figure,'WIP: 00%','Start');
for t = 1:T
    
    % Nodes or members of the current cluster
    Node1 = Tree(t,1);
    Node2 = Tree(t,2);
    
    % First member of the current cluster
    if Node1 <= N
        
        % Index of the current node
        i = find(Roots==Node1);
        
        % Indices permutation
        Index = Indices(i);
        Indices(Indices==Node1) = Index;
        Indices(i) = Node1;
        
        % Coordinates of the current root
        x1 = AbscissaRoots(i);
        y1 = OrdinatesRoots(i);
        
    else
        
        % Coordinates of the current node
        x1 = AbscissaNodes(Node1-N);
        y1 = OrdinatesNodes(Node1-N);
        
    end
    
    % Second member of the cluster
    if Node2 <= N
        
        % Index of the current node
        i = find(Roots==Node2);
        
        % Indices permutation
        Index = Indices(i);
        Indices(Indices==Node2) = Index;
        Indices(i) = Node2;
        
        % Coordinates of the current root
        x2 = AbscissaRoots(i);
        y2 = OrdinatesRoots(i);
        
    else
        
        % Coordinates of the current node
        x2 = AbscissaNodes(Node2-N);
        y2 = OrdinatesNodes(Node2-N);
        
    end
    
    % Bars update
    switch M
        
        case 1
            
            % Update of scalar bar
            for i = 1:N
                v = Indices(i);
                Index = floor( (Matrix(v)-Min) / (Max-Min) *(64-1)+1);
                Color = Colors(Index,:);
                set(Bars(i),'Ydata',Matrix(v),'Facecolor',Color);
            end
            
        otherwise
            
            % Update of monocolor multiple bars
            for m = 1:M
                set(Bars(m),'Ydata',Matrix(Indices,m));
            end
            
    end
    
    % Coordinates of the current node
    y = Tree(t,3);
    AbscissaNodes(t) = mean([x1 x2]);
    OrdinatesNodes(t) = y;
    
    % Data label
    Xticklabel = arrayfun(@(n)sprintf('%u',Indices(n)),1:N,'UniformOutput',false);
    set(Axes,'Xtick',1:N,'XtickLabel',Xticklabel);
    
    % Creation of a new graphical node
    plot(Axes(2),[x1 x1 x2 x2],[y1 y y y2],'b');
    set(Axes(2),'Xlim',Xlim,'Ylim',Ylim);
    
    drawnow();
    
    % Spinner
    Spinner(Figure,sprintf('WIP: %02d%%',round(100*t/T)),'Update');
    
end
uistack(Text,'top');
% Clusters rectangles colors
Colors = {'y','m','c','g','r','b','w','k'};
% Clusters rectangles
for c = 1:Number
    
    if c == 1
        
        % First rectangle
        x1 = 0.5;
        xc1 = find(Indices==Clusters{c}(end));
        xc2 = find(Indices==Clusters{c+1}(1));
        x2 = (xc1+xc2)/2;
        
    elseif c < Number
        
        % Intermediate rectangles
        xc1 = find(Indices==Clusters{c-1}(end));
        xc2 = find(Indices==Clusters{c}(1));
        x1 = (xc1+xc2)/2;                
        xc1 = find(Indices==Clusters{c}(end));
        xc2 = find(Indices==Clusters{c+1}(1));
        x2 = (xc1+xc2)/2;
        
    else
        
        % Last rectangle
        xc1 = find(Indices==Clusters{c-1}(end));
        xc2 = find(Indices==Clusters{c}(1));
        x1 = (xc1+xc2)/2;  
        x2 = N+0.5;
        
    end
    
    % Cluster
    Rectangle = ...
        patch('Faces',     1:4,...
              'Vertices',  [x1 0; x2 0; x2 Limit ; x1 Limit],...
              'FaceColor', Colors{mod(c-1,numel(Colors))+1},...
              'EdgeColor', 'None',...
              'FaceAlpha', 0.5);
    uistack(Rectangle,'bottom');
    
end
% Spinner
Spinner(Figure,'100%','Stop');
    
end
% Spinner
function Spinner(Figure,Message,Action)
persistent jObj
Position = [10,10,80,80];
switch lower(Action)
    
    case 'start'
        
        if ~isempty(jObj)
            jObj.stop;
            jObj.setBusyText(Message);
            [~, hContainer] = javacomponent(jObj.getComponent, Position, Figure);
            delete(hContainer);
        end
        
        try
            
            % R2010a and newer
            iconsClassName = 'com.mathworks.widgets.BusyAffordance$AffordanceSize';
            iconsSizeEnums = javaMethod('values',iconsClassName);
            SIZE_32x32 = iconsSizeEnums(2);  % (1) = 16x16,  (2) = 32x32
            jObj = com.mathworks.widgets.BusyAffordance(SIZE_32x32,Message);
            jObj.getComponent.setBackground(java.awt.Color(1, 1, 1));
            jObj.setPaintsWhenStopped(0);
            
        catch
            
            % R2009b and earlier
            redColor   = java.awt.Color(1,0,0);
            blackColor = java.awt.Color(0,0,0);
            jObj = com.mathworks.widgets.BusyAffordance(redColor, blackColor);
            
        end
        
        jObj.setPaintsWhenStopped(true);  % default = false
        jObj.useWhiteDots(false);         % default = false (true is good for dark backgrounds)
        javacomponent(jObj.getComponent, Position, Figure);
        jObj.start;
        drawnow();
        
    case 'update'
        
        jObj.setBusyText(Message);
        
    case 'stop'
        
        if isempty(jObj)
            return
        end
        
        jObj.stop;
        jObj.setBusyText(Message);
        pause(0.5);
        [~, hContainer] = javacomponent(jObj.getComponent, Position, Figure);
        delete(hContainer);
        
    otherwise
        error('Unknown spinner action : "%s".',Action);
        
end
end


