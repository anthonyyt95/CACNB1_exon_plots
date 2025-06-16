%% Merge RNAseq files


directory = 'C:\Users\antho\Documents\NYU\NYU Langone\PhD\Feske Lab\Experiments\03.31.19_ICT Analysis\Data\2019.11.08_RandomEffectors\GSE131359';
flist = dir(directory);

oval = [];
ogenes = {};
ocond = {};
first = 1;
for i = [1:length(flist)]
    if flist(i).bytes == 0
        continue
    end
    
    folder = flist(i).folder;
    fname = flist(i).name;
    path = [folder, '\', fname];
    
    data = importdata(path);
    values = data.data;
    
    if first == 1
        genes = data.textdata([2:end],1);
        first = 0
    end
    
    oval = [oval values];
    ocond = [ocond fname];
end

heading = ['Genes', ocond];
output = [genes num2cell(oval)];
output = [heading; output];

xlswrite('temp.xlsx',output);

clear directory flist ogenes oval ocond heading i folder fname path data
clear values first genes

%% Organizes the downloaded data from GEO (GSE116347)

fname = 'C:\Users\antho\Documents\NYU\NYU Langone\PhD\Feske Lab\Experiments\03.31.19_ICT Analysis\Data\2019.11.08_RandomEffectors\GSE137057_tpm_geneIDd.xlsx';
data = importdata(fname);
acquireICT = 1;
list = ICT797;

tic
output = {};
output.values = data.data.Sheet1;
output.cells = data.textdata.Sheet1(1,[2:end]);
output.genes = data.textdata.Sheet1([2:end],1);
output.annot = data.textdata.Sheet2';
output.description = {};

data = [];
genes = {};
missing = {};
if acquireICT == 1
    for i = [1:length(list(:,1))]
        aliases = strsplit(list{i,1},' ');
        
        match = 0;
        for j = [1:length(aliases)]
            gene = aliases{j};
            loc = find(strcmp(lower(output.genes),lower(gene)));
            if not(isempty(loc))
                data = [data; output.values(loc(1),:)];
                genes = [genes; aliases{1}];
                match = 1;
                break
            end
        end
        if match == 0
            missing = [missing; aliases{1}];
        end
    end
    output.genes = genes;
    output.values = data;
    output.missing = missing;
end
toc

clear fname data acquireICT data genes missing i j 
clear aliases match loc list


%% Assigns ENSEMBLE IDs to gene names

fname = 'C:\Users\antho\Documents\NYU\NYU Langone\PhD\Feske Lab\Pan-Experiment Data\Gene Lengths for TPM\mmusculus.xlsx';
data = 'C:\Users\antho\Documents\NYU\NYU Langone\PhD\Feske Lab\Experiments\03.31.19_ICT Analysis\Data\2019.11.08_RandomEffectors\GSE116288.xlsx';
remove_ensembl_period = 1;

ref = table2cell(readtable(fname));
data = importdata(data);
genes = data.textdata.Sheet1([2:end],1);

tic
% Removes period from ENSEMBL ID if present
if remove_ensembl_period == 1
    for i = [1:length(genes(:,1))]
        val = genes{i,1};
        loc = strfind(val, '.');
        val = val([1:loc-1]);
        genes{i,1} = val;
    end
end

% Matches ENSEMBL ID with 'ref' and retrieves gene name from 'ref'
output = {};
ref_genes = lower(ref(:,1));
for i = [1:length(genes(:,1))]
    gene = lower(genes{i,1});
    loc = find(strcmp(ref_genes, gene));
    i
    if not(isempty(loc))
        output = [output; upper(ref{loc,3})];
    else
        output = [output; upper(gene)];
    end
end
toc

clear fname data ref ref_genes i gene loc

%% Calculates TPM

values = output.values;
genes = output.genes;
conditions = output.conditions;
ref_file = 'C:\Users\antho\Documents\NYU\NYU Langone\PhD\Feske Lab\Pan-Experiment Data\Gene Lengths for TPM\mmusculus.xlsx';

% Acquires gene lengths for genes ----------------------------------------
tic
% Proecss gene lengths file
glengths = importdata(ref_file);
rgenes = glengths.textdata([2:end],3);
rens = glengths.textdata([2:end],2);
rlengths = glengths.data;
ref = [rgenes rens num2cell(rlengths)];

% Convert to lower string
ref_genes = lower(ref(:,1));

% Creates a list of gene lengths associated with each item in 'genes'
exclude = [];
list = [];
for i = [1:length(genes(:,1))]
    gene = lower(string(genes{i}));
    
    loc = find(strcmp(ref_genes, gene));
    if isempty(loc)
        exclude = [exclude; i];
        glen = NaN;
    else
        glen = ref{loc,6};
    end
    
    list = [list; glen];
end
genes(exclude,:) = [];
values(exclude,:) = [];
list(exclude,:) = [];

% Performs TPM calculation
output_values = [];
for i = [1:length(values(1,:))]
    val = values(:,i);
    scaling_factor = sum(val)/(10^6);

    val = val./(list/(10^3));
    val = val/scaling_factor;
    output_values = [output_values val];
end

% Saves file
output = {};
output = [genes num2cell(output_values)];
delete 'output.xlsx'
xlswrite('output.xlsx', output);
toc


clear gene genes glen glengths i loc ref ref_genes rens
clear rgenes rlengths conditions ref_file exclude list output output_values
clear scaling_factor val values

%% Acquires ICT list

