function [M]=hiCtsv_to_MATLABcsv(Resolution,ChrSize,tsv_IN,csv_M_OUT)

%tsv_IN ... output from HiCExplorer --> convert to .csv for matlab to read
%M_tsv=load(sprintf("../hicMatrix/Chr%d/%d/hiCmatrix_GSM862724.corrected.tsv",ChrNr,Resolution));

tsv_IN=load(tsv_IN);
BINS=Bin_Map(ChrSize,Resolution);
row=tsv_IN(1,2:3);

if mod(row(1,2),2)~=mod(BINS(1,3),2)
    BINS(:,2:3)=BINS(:,2:3)-1;
end

%BINS(:,3)=BINS(:,3);
%BINS(:,2)=BINS(:,2)-1;
%BINS=BINS(First_NonZero_Bin:Last_NonZero_Bin,:);

format long
M=zeros(size(BINS,1),size(BINS,1));

for k=1:size(tsv_IN,1)
    row=tsv_IN(k,2:3);
    col=tsv_IN(k,5:6);

    BIN_r=BINS(row(1)==BINS(:,2)&row(2)==BINS(:,3),1);
    BIN_c=BINS(col(1)==BINS(:,2)&col(2)==BINS(:,3),1);  
        
    M(BIN_r,BIN_c)=tsv_IN(k,end);
    M(BIN_c,BIN_r)=tsv_IN(k,end);
    
end

csvwrite(csv_M_OUT,M)
M=load(csv_M_OUT);
