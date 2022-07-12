function matlabM_to_homer(BINs,ChrNr,HOME_TAG_DIR,HOMER_OUT,M)
format long
OUT_=fopen(HOMER_OUT,'w');
vector=[repmat(ChrNr,size(M(1,:)));(BINs(:,2))'];
vector=reshape(vector,1,2*size(vector,2));
string=sprintf(repmat('%d-%d\t',size(M(1,:))),vector);
fprintf(OUT_,'HiCMarix\t(directory=%s)\t%s\n',HOME_TAG_DIR,string);clear string
for row=1:size(BINs,1)
    ROW=BINs(row,1);
    string=sprintf(repmat('%f\t',size(M(1,:))),M(row,:));string=string(1:end-1);
    %string=sprintf(repmat('%s\t',size(M(1,:))),strrep(num2str(M(row,:)),'.',','));
    fprintf(OUT_,'%d-%d\t%d-%d\t%s\n',ChrNr,BINs(BINs(:,1)==ROW,2),ChrNr,BINs(BINs(:,1)==ROW,2),string);
end
fclose(OUT_);


end