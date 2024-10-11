%grep KO GENE
clear
[filename1,pathname1]=uigetfile({'\*.txt'});
microb=readtable([pathname1,filename1]);
   T1=microb.Properties.VariableNames;
   ta=microb(:,1);
   T=microb(:,2:end);
 A=table2array(T);clear T;
%db=readcell('D:\temp\DBGHG\2KO\new3.txt');
[filename,pathname]=uigetfile({'\*.txt'});
db=readtable([pathname,filename]);

dbko=db(:,1);
dbgene=db(:,2);

ta=table2cell(ta);
dbko=table2cell(dbko);
alpha=ismember(ta,dbko);
beta=ismember(dbko,ta);

geneta=ta(alpha);
genenum=A(alpha,:);

genenum=num2cell(genenum);
x1=[geneta,genenum];
gene=cell2table(x1);
gene.Properties.VariableNames=T1;

write(gene,[pathname1,'geneKO.txt'])