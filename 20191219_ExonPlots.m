%% Plots exon counts with error bars

filename = 'E:\Documents\NYU\NYU Langone\PhD\Feske Lab\Experiments\04.10.19_(Serap)CACN project\Misc\2020.08.16_ExonData5\CACNA1I-202.csv';
data = importdata(filename);
%cond1 = 'CD4_0H';
%cond2 = 'CD4_24H';
% cond3 = 'CD4_24H';
 cond1 = 'BRAIN';
 cond2 = 'T_CELLS';
%cond1 = 'BRAIN';
%cond2 = '_unstimulated';
bColor = [0 0 0.5;
    1 0.4 0.4,
    0 0.8 0];
bWidth = 1;

% Cleans up data
values = data.data(:,[2:end]);
values = log10(values+1);
heading = data.textdata(1,[3:end]);
exonNum = data.data(:,1);

% Acquires reference & T-cell indices
ind1 = [];
ind2 = [];
% ind3 = [];
for i = [1:length(heading(1,:))]
    loc = strfind(lower(heading{1,i}), lower(cond1));
    if not(isempty(loc))
        ind1 = [ind1 i];
    end
    
    loc = strfind(lower(heading{1,i}), lower(cond2));
    if not(isempty(loc))
        ind2 = [ind2 i];
    end
    
%     loc = strfind(lower(heading{1,i}),lower(cond3));
%     if not(isempty(loc))
%         ind3 = [ind3 i];
%     end
end

% Subsets values
cond1Data = values(:,ind1);
cond1Mean = mean(cond1Data,2);
cond1STD = std(cond1Data,0,2);

cond2Data = values(:,ind2);
cond2Mean = mean(cond2Data,2);
cond2STD = std(cond2Data,0,2);

% cond3Data = values(:,ind3);
% cond3Mean = mean(cond3Data,2);
% cond3STD = std(cond3Data,0,2);

exonNum = exonNum';

% Graphs bar figure
figure
hold
xval = 1;
xtic = [];
% barData = [cond1Mean,cond2Mean,cond3Mean];
% errData = [cond1STD,cond2STD,cond3STD];
for i = [1:length(barData(:,1))]
    yval1 = barData(i,1);
    yval2 = barData(i,2);
%     yval3 = barData(i,3);
    b1 = bar(xval,yval1);
    xval = xval + 1;
    val = xval;
    b2 = bar(xval,yval2);
%     xval = xval + 1;
%     b3 = bar(xval,yval3);
    
    xtic = [xtic, val];
    
    b1.EdgeColor = 'none';
    b2.EdgeColor = 'none';
%     b3.EdgeColor = 'none';
    b1.FaceColor = bColor(1,:);
    b2.FaceColor = bColor(2,:);
%     b3.FaceColor = bColor(3,:);
    b1.BarWidth = bWidth;
    b2.BarWidth = bWidth;
%     b3.BarWidth = bWidth;
    
    
    err1 = errData(i,1);
    err2 = errData(i,2);
%     err3 = errData(i,3);
    e1 = plot([xval-1 xval-1],[yval1,yval1+err1]);
    e2 = plot([xval xval],[yval2,yval2+err2]);
%     e3 = plot([xval,xval],[yval3,yval3+err3]);
    e1.Color = bColor(1,:);
    e2.Color = bColor(2,:);
%     e3.Color = bColor(3,:);

    xval = xval + 1.5;
end

ylim([0,3.5]);
set(gca,'TickDir','out');
set(gca,'FontSize',8);
set(gca,'XTick',xtic);
set(gca,'XTickLabel',num2cell(exonNum));
set(gca,'YTick',[0 1 2 3]);
set(gca,'YTickLabel',{'10^0','10^1','10^2','10^3'});
set(gca,'LineWidth',2);
%%
clear b1 b2 barData bColor bWidth cond1 cond1Data cond1Mean cond1STD
clear cond2 cond2Data cond2Mean cond2STD data e1 e2 err1 err2 errData exonNum
clear filename heading i ind1 ind2 loc values xtic xticVal xval yval1 yval2

%% Create exon diagram based on exon data

filename = 'E:\Documents\NYU\NYU Langone\PhD\Feske Lab\Experiments\04.10.19_(Serap)CACN project\Misc\2020.12.28_ExonData_CACNB1\Cacnb1_205_MusMusculus.csv';
yval = 100;

% Downloads & organizes data
data = importdata(filename);
temp = {};
header = strsplit(data{1,:},',','CollapseDelimiters',false);
for i = [3:length(data(:,1))-1]
    val = data{i,1};
    val = val([2:end-1]);
    val = strsplit(val,'","','CollapseDelimiters',false);
    val = val(1,[1:5]);
    
    if not(isequal(val{1,1},' '))
        temp = [temp; val];
    end
end

% Calculates coordinates
xval = [];
for i = [1:length(temp(:,1))]
%     startVal = temp{i,3};
%     endVal = temp{i,4};
    
     startVal = regexprep(temp{i,3},',','');
     endVal = regexprep(temp{i,4},',','');
     startVal = str2double(startVal);
     endVal = str2double(endVal);
    
    xval = [xval; startVal, endVal];
end
base = xval(1,1);
xval = xval - base;
xval = abs(xval);

% Creates graph
figure
hold
l = plot([xval(1,1),xval(end,2)],[yval/2,yval/2])
l.LineWidth = 2;
l.Color = [0 0 0];
for i = [1:length(xval(:,1))]
    leftVal = xval(i,1);
    rightVal = xval(i,2);
    fillX = [leftVal,leftVal,rightVal,rightVal];
    fillY = [0, yval, yval, 0];
    f = fill(fillX,fillY,[0 0 0]);
end
set(gca,'Color','none')
set(gca,'visible','off')

clear base bColor bWidth cond1 cond2 data endVal f filename fillX fillY
clear header i l leftVal list output rightVal startVal temp val xval yval

%% Acquires unique exon IDs

input = struc

strucFields = fields(input);
allExons = {};
for i = [1:length(strucFields)]
    eval(['data = input.',strucFields{i},';']);
    
    allExons = [allExons; data];
end

allExons = sortrows(allExons,2,'desc');
uniqueExons = unique(allExons(:,1), 'stable');
occurExons = [];
exonLoc = {};
for i = [1:length(uniqueExons(:,1))]
    exon = uniqueExons{i,1};
    loc = find(strcmp(allExons(:,1),exon));
    occurExons = [occurExons; length(loc)];
    exonLoc = [exonLoc; allExons(loc(1),[2,3])];
    
end

output = [uniqueExons, exonLoc, num2cell(occurExons)];

%% Acquires ENSME IDs for each transcript

data = input;
transcripts = {'Cacnb1-201','Cacnb1-202','Cacnb1-203','Cacnb1-204','Cacnb1-205'};

output = nan(20,length(transcripts));
output = num2cell(output);
for i = [1:length(transcripts)]
    transcript = transcripts{i};
    loc = find(strcmp(data(:,3),transcript));
    temp = data(loc,:);
    temp = sortrows(temp,2,'asc');
    temp = unique(temp(:,1),'stable');
    
    for j = [1:length(temp)]
        val = temp{j,1};
        val = strsplit(val,'.');
        temp{j,1} = val(1);
    end

    output([1:length(temp)],i) = temp(:,1);
end

%% Combines data according to a reference Exon structure

barData = input1;
errData = input2;
ref = input3;

uniqueID = ref(:,1);
for i = [1:length(uniqueID(:,1))]
    uniqueID{i,1} = num2str(uniqueID{i,1});
end
ref(:,1) = uniqueID;
uniqueID = unique(uniqueID,'stable');

barOut = [];
errOut = [];
for i = [1:length(uniqueID(:,1))]
    id = uniqueID{i,1};
    loc = find(strcmp(ref(:,1),id));
    
    barTemp = [];
    errTemp = []
    for j = [1:length(loc)]
        exon = ref{loc(j),2};
        loc2 = find(strcmp(barData(:,1),exon));
        
        barTemp = [barTemp; barData(loc2(1),[3:5])];
        errTemp = [errTemp; errData(loc2(1),[3:5])];
    end
    
    barTemp = cell2mat(barTemp);
    errTemp = cell2mat(errTemp);
    barOut = [barOut; mean(barTemp,1)];
    errOut = [errOut; mean(errTemp,1)];
end

%% Plots exon counts with error bars + individual points

filename = 'E:\Documents\NYU\NYU Langone\PhD\Feske Lab\Experiments\04.10.19_(Serap)CACN project\Misc\2020.06.19_ExonData3\Cacna1a_213_TCELLS_VS_BRAIN.csv';
data = importdata(filename);
%cond1 = 'CD4_0H';
%cond2 = 'CD4_24H';
% cond3 = 'CD4_24H';
%cond1 = 'BRAIN';
%cond2 = 'T_CELLS';
cond1 = 'BRAIN';
cond2 = 'CD4_0H';
bColor = [0 0 0.5;
    1 0.4 0.4,
    0 0.8 0];
bWidth = 1;

% Cleans up data
values = data.data(:,[2:end]);
values = log10(values+1);
heading = data.textdata(1,[3:end]);
exonNum = data.data(:,1);

% Acquires reference & T-cell indices
ind1 = [];
ind2 = [];
% ind3 = [];
for i = [1:length(heading(1,:))]
    loc = strfind(lower(heading{1,i}), lower(cond1));
    if not(isempty(loc))
        ind1 = [ind1 i];
    end
    
    loc = strfind(lower(heading{1,i}), lower(cond2));
    if not(isempty(loc))
        ind2 = [ind2 i];
    end
    
%     loc = strfind(lower(heading{1,i}),lower(cond3));
%     if not(isempty(loc))
%         ind3 = [ind3 i];
%     end
end

% Subsets values
cond1Data = values(:,ind1);
cond1Mean = mean(cond1Data,2);
cond1STD = std(cond1Data,0,2);

cond2Data = values(:,ind2);
cond2Mean = mean(cond2Data,2);
cond2STD = std(cond2Data,0,2);

exonNum = exonNum';

% Graphs bar figure
f = figure
hold
xval = 1;
xtic = [];
errData = [cond1STD,cond2STD];
for i = [1:length(cond1Data(:,1))]
    yval1 = cond1Mean(i,1);
    yval2 = cond2Mean(i,1);
    dots1 = cond1Data(i,:)';
    dots2 = cond2Data(i,:)';
    
    b1 = bar(xval,yval1);
    for j = [1:length(dots1)]
        s = scatter(xval,dots1(j),'o');
        s.MarkerEdgeColor = [0.5 0.5 0.5];
        s.MarkerFaceColor = [0.5 0.5 0.5];
        s.MarkerFaceAlpha = 0.5;
        s.SizeData = 5;
    end
    xval = xval + 1;
    val = xval;
    b2 = bar(xval,yval2);
    for j = [1:length(dots2)]
        s = scatter(xval,dots2(j),'o');
        s.MarkerEdgeColor = [0.5 0.5 0.5];
        s.MarkerFaceColor = [0.5 0.5 0.5];        
        s.MarkerFaceAlpha = 0.5;
        s.SizeData = 5;
    end
    
    xtic = [xtic, val];
    
    b1.EdgeColor = 'none';
    b2.EdgeColor = 'none';
    b1.FaceColor = bColor(1,:);
    b2.FaceColor = bColor(2,:);
    b1.BarWidth = bWidth;
    b2.BarWidth = bWidth;
   
    
    err1 = errData(i,1);
    err2 = errData(i,2);
    e1 = plot([xval-1 xval-1],[yval1,yval1+err1]);
    e2 = plot([xval xval],[yval2,yval2+err2]);
    e1.Color = bColor(1,:);
    e2.Color = bColor(2,:);

    xval = xval + 1.5;
end

ylim([0,3.5]);
set(gca,'TickDir','out');
set(gca,'FontSize',8);
set(gca,'XTick',xtic);
set(gca,'XTickLabel',num2cell(exonNum));
set(gca,'YTick',[0 1 2 3]);
set(gca,'YTickLabel',{'10^0','10^1','10^2','10^3'});
set(gca,'LineWidth',2);
        


    





