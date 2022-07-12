function [r,pval,AUC,ACC,F1] = perf_metrics(PRED,TRUE)

[r, pval]=corrcoef(PRED',TRUE);
r=r(1,2);pval=pval(1,2);

AB_DDA=zeros(length(PRED),1);AB_hic=zeros(length(TRUE),1);
AB_DDA(PRED>0.5)=1;
AB_hic(TRUE>0.5)=1;

[~,~,~,AUC] = perfcurve(AB_hic,PRED',1);

TP=0;FP=0;TN=0;FN=0;
for i=1:length(AB_DDA)
    if(AB_hic(i)==1 && AB_DDA(i)==1)
        TP=TP+1;
    elseif(AB_hic(i)==0 && AB_DDA(i)==1)
        FP=FP+1;
    elseif(AB_hic(i)==0 && AB_DDA(i)==0)
        TN=TN+1;
    else
        FN=FN+1;
    end
end

ACC=(TP+TN)/(TP+TN+FP+FN);
F1=TP/(TP+.5*(FP+FN));

end

