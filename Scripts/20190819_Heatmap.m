%% Create heatmap
%
% Input:
%   'data'      X-by-Y cell matrix containing heatmap values
%                   col1: y-variable names
%                   row1: x-variable names
%

data = info;

yvar = data([2:end],1);
xvar = data(1,[2:end]);
data = rescale(data,1,1000);
% data = cell2mat(data([2:end],[2:end]));



h = heatmap(data);
% h.XDisplayLabels = xvar;
% h.YDisplayLabels = yvar;
h.GridVisible = false;
h.CellLabelColor = 'none';

set(gca,'Colormap',cmap)

clear data h values xvar yvar



