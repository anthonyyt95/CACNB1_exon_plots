%% Acquires datapoints associated with a list of genes & a list of celltypes
%
% Input
%   'cellList'      X-by-1 cell matrix of celltypes to include. Rows
%                   represent replicates of celltypes to be averaged;
%                   Columns represent different celltypes.
%   'geneList'      X-by-1 cell list of gene names (+ aliases) to be
%                   retrieved

data = ThSubsetsMouse.values;
genes = ThSubsetsMouse.genes;
cells = ThSubsetsMouse.conditions;
cellList = clist;
geneList = cacnList;

% Acquires celltypes (& performs averaging)
data1 = [];
output_cells = {};
for i = [1:length(cellList(1,:))]
    
    val = [];
    for j = [1:length(cellList(:,1))]
        cell = cellList{j,i};
        if isempty(cell)
            break
        end
        
        loc = find(strcmp(lower(cells), lower(cell)));
        if not(isempty(loc))
            val = [val data(:,loc)];
        end
    end
    val = mean(val,2);
    data1 = [data1 val];
    output_cells = [output_cells cell];
end

% Acquires gene information
output_data = [];
output_genes = {};
missing = {};
for i = [1:length(geneList(:,1))]
    aliases = strsplit(geneList{i,1}, ' ');
    
    match = 0;
    for j = [1:length(aliases)]
        gene = aliases{j};
        loc = find(strcmp(lower(genes), lower(gene)));
        if not(isempty(loc))
            output_data = [output_data; data1(loc,:)];
            output_genes = [output_genes; aliases{1}];
            match = 1;
            break
        end
    end
    
    if match == 0
        missing = [missing; aliases{1}];
    end
end
output.values = output_data;
output.cells = output_cells;
output.genes = output_genes;
output.missing = missing;

clear aliases cell cellList cells data data1 gene geneList genes i j loc
clear match missing output_cells output_data output_genes val 



%% Heatmap - scale within genes

data = output.values([1:26],:);
genes = output.genes([1:26],1);
cells = output.cells;
%%
data = data + 1;
data = log10(data);
data = rescale(data,1,10);
% for i = [1:length(genes(:,1))]
%     val = data(i,:);
%     val = rescale(val,1,100);
%     data(i,:) = val;
% end


h = heatmap(data);
h.CellLabelColor = 'none';
%h.XDisplayLabels = cells;
%h.YDisplayLabels = genes;
h.FontSize = 8;
h.Colormap = cmap;
%%
clear cells data do_average doAvg_doSort genes hdata

