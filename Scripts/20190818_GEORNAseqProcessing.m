%% Performs TPM normalization of gene counts
%
%
% Input:
%   'filename'          Path: xlsx file
%                           - rows: gene names
%                           - columns: celltypes (ie sample IDs)
%   'gene2len'          file path for gene names file
%                           - hsapiens: 'C:\Users\antho\Documents\NYU\NYU Langone\PhD\Feske Lab\Pan-Experiment Data\Gene Lengths for TPM\hsapiens.xlsx'
%                           - mmusculus: 'C:\Users\antho\Documents\NYU\NYU Langone\PhD\Feske Lab\Pan-Experiment Data\Gene Lengths for TPM\hsapiens.xlsx'
%
% Output:
%   'output'            Structure
%                           output.data:        X-by-Y double of TPM values
%                           output.genes:       X-by-1 cell list of genes
%                           output.celltypes:   1-by-Y cell list of celltypes/samples


filename = 'C:\Users\antho\Documents\NYU\NYU Langone\PhD\Feske Lab\Experiments\06.26.18_Misc\Review Figures_Pie chart\GEO data\GSE135803_norm.xlsx';
gene2len = 'C:\Users\antho\Documents\NYU\NYU Langone\PhD\Feske Lab\Pan-Experiment Data\Gene Lengths for TPM\hsapiens.xlsx';


tic
data = importdata(filename);
celltypes = data.textdata.Sheet1(1,[2:end]);
genes = data.textdata.Sheet1([2:end],1);
values = data.data.Sheet1;
annot = data.textdata.Sheet2';

gene2len = importdata(gene2len);
rgenes = gene2len.textdata([2:end],3);
rlen = gene2len.data(:,4);
rID = gene2len.textdata([2:end],1);



%% Acquire genes associated with ENS IDs

gref = [rID rgenes];
gref = sortrows(gref,1);

num_gref_id = {};
for i = [1:length(gref(:,1))]
    gref_id = gref{i,1};
    gref_id = gref_id([8:end]);
    gref_id = str2num(gref_id);
    
    num_gref_id{gref_id,1} = gref{i,1};
    num_gref_id{gref_id,2} = gref{i,2};    
end
    
gene_names = {};
exclude = [];
for i = [1:length(genes(:,1))]
    gene = genes{i,1};
    i
    % Removes period from ENS ID
    pdloc = strfind(gene,'.');
    if not(isempty(pdloc))
        gene = gene([1:pdloc-1]);
    end

    % Converts ENS ID into a form that can be indexed
    gene = str2num(gene([8:end]));
    gene_assoc = num_gref_id{gene,2};
    
    if not(isempty(gene_assoc))
        gene_names = [gene_names; gene_assoc];
    else
        exclude = [exclude; i];
    end
end
genes = gene_names;
values(exclude,:) = [];

output = [genes num2cell(values)];
xlswrite('temp.xlsx',output);

load handel
sound(y,Fs)


%% Assigns genes to transcripts

% Performs TPM normalization
ref = [lower(rgenes) num2cell(rlen)];
ref = sortrows(ref,1);
gannot = alph_search_annot(ref,1,0);

lengths = [];
exclude = [];
for i = [1:length(genes(:,1))]
    gene = genes{i,1};

    isperiod = strfind(gene,'.');
    islinc = strfind(gene,'LINC');
    isloc = strfind(gene,'LOC');
    if not(isempty(isperiod)) | not(isempty(islinc)) | not(isempty(isloc))
        exclude = [exclude; i];
        continue
    end
    i
    % Searches for gene using alphabetized annotation
    first_let = gene(1);
    start_ind = find(strcmp(gannot(:,1),lower(first_let)));
    if start_ind == length(gannot(:,1))
        search_start = gannot{start_ind,2};
        search_end = length(ref(:,1));
    else
        search_start = gannot{start_ind,2};
        search_end = gannot{start_ind+1,2};
    end

    temp_rgenes = ref([search_start:search_end],1);
    loc = find(strcmp(lower(gene), lower(temp_rgenes)));

    if not(isempty(loc)) & length(loc)==1
        lengths = [lengths; ref{loc,2}];
    else
        exclude = [exclude; i];
    end
end
genes(exclude,:) = [];
values(exclude,:) = [];


% Performs count to TPM conversion
for i = [1:length(values(:,1))]  
    len = lengths(i,1)/1000;
    val = values(i,:);
    val = val/len;

    values(i,:) = val;
end

for i = [1:length(values(1,:))]
    val = values(:,i);
    permil = sum(val)/1e6;

    val = val/permil;
    values(:,i) = val;
end


output = [genes num2cell(values)];
xlswrite('temp.xlsx', output);

load handel
sound(y,Fs)





