function [PC_OUT,whichPC] = Norm_PC(PC_IN,BED_FILE,ChrNr,BINs,EXCLUDE_PCA)
a1=0;a2=0.5;
b1=0.5;b2=1;

FN=strrep(BED_FILE,'.bed','.mat');
%if exist(FN,'file')==0
    [H3K4me1_DENSITY] = H3K4me1Density(BINs,ChrNr,BED_FILE);H3K4me1_DENSITY=H3K4me1_DENSITY(:,2);
    H3K4me1_DENSITY= (H3K4me1_DENSITY-min(H3K4me1_DENSITY))/(max(H3K4me1_DENSITY) - min(H3K4me1_DENSITY));
    save(FN,'H3K4me1_DENSITY');
%end
load(FN);
if ~isempty(EXCLUDE_PCA)
            H3K4me1_DENSITY(EXCLUDE_PCA(1):EXCLUDE_PCA(2))=[];
end



Nr_consider=3;

CC=nan(1,Nr_consider);
PC_OUT=nan(size(PC_IN,1),Nr_consider);
for k=1:Nr_consider
    PC_IN_2=PC_IN(:,k); PC_IN_1=[];
    
    PC_IN_1=PC_IN(:,k);
    PC_IN_1=filloutliers(PC_IN_1,'nearest');
    % center around zero mean and unit var
    PC_IN_1 = (PC_IN_1-mean(PC_IN_1))/std(PC_IN_1);
    
    C=corrcoef(PC_IN_1,H3K4me1_DENSITY);
    CC(k)=C(1,2);
    
    if CC(k)<0
        PC_IN_1=-PC_IN_1;
    end
    
    EV=PC_IN_1;X=nan(size(EV));
    for n=1:length(EV)
        if EV(n)<0
            X(n)=(a2-a1)*(EV(n)-min(EV))./(max(EV)-min(EV))+a1;
        elseif EV(n)>=0
            X(n)=(b2-b1)*(EV(n)-min(EV))./(max(EV)-min(EV))+b1;
        end
    end
    PC_OUT(:,k)=X;clear X
    
    clear PC_IN_1
    
end
%  CC

whichPC=find(abs(CC(1:Nr_consider))==nanmax(abs(CC(1:Nr_consider))));

PC_OUT=PC_OUT(:,whichPC);

end


