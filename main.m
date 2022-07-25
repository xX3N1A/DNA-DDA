%% load variables
load('vars.mat')
FIG_PATH='Figures';

EXCLUDE(EXCLUDE==0)=[];EXCLUDE=unique(EXCLUDE);
if EXCLUDE~=0
    HiCM(EXCLUDE,:)=nan;HiCM(:,EXCLUDE)=nan;
    HiCP(EXCLUDE,:)=nan;HiCP(:,EXCLUDE)=nan;
end

HiCP=HiCP(BINs(:,1),BINs(:,1));
HiCM=HiCM(BINs(:,1),BINs(:,1));

if ~isempty(EXCLUDE_PCA)
    HiCP(EXCLUDE_PCA(1):EXCLUDE_PCA(2),:)=[];
    HiCP(:,EXCLUDE_PCA(1):EXCLUDE_PCA(2))=[];
end

%% generate **1D DNA walk** and save to file <FN_ASCII>

makeTS=0;
ASCII_DIR=('..');
FN_ASCII=(sprintf('%s/RandomWalk_Chr22.ascii',ASCII_DIR));
% vergleiche ascii file
if makeTS==1
    DATA=load(sprintf('../../../../data/PERL_SEQ_ALL/Chr%d_BINALL.perl.ascii',ChrNr));
    unix(sprintf('touch %s',FN_ASCII))
    TS=nan(Resolution,size(BINs,1));
    for w_x=1:size(BINs,1)
        X=DATA(BINs(w_x,2):(BINs(w_x,3)-1));
        if sum(X)~=0
            TS(1,w_x)=X(1);
            for n=2:length(X)
                if X(n-1)==0
                    TS(n,w_x)=TS(n-1,w_x);
                elseif X(n-1)==1 || X(n-1)==3
                    TS(n,w_x)=TS(n-1,w_x)-1;
                elseif X(n-1)==2 || X(n-1)==4
                    TS(n,w_x)=TS(n-1,w_x)+1;
                end
            end
        end
    end
    unix(['touch ',FN_ASCII]); % create file
    fid = fopen(FN_ASCII,'w'); % open file
    string=repmat('%15.10f ',1,size(TS,2));
    fprintf(fid,[string ' \n'],TS');% write to file...format in spalten! also eine zeitserie ist die erste ZEILE (tsx length)
    fclose(fid);
    
    clear X*;clear TS;
    
    
    X=load(FN_ASCII);X=X(:,BINs(:,1));
    h=plot(X(:,:),'linewidth',1);
    xlabel('steps [nucleotides]','fontsize',20,'interpreter','latex')
    ylabel('$x(s)$','fontsize',20,'interpreter','latex')
    title(sprintf('Random Walk DNA chr%d:%d-%d',ChrNr,BINs(1,2),BINs(end,3)),'fontsize',20,'interpreter','latex')
end

%% generate bash script to run DDA
OD_DDA=sprintf('DDA_OUT',ChrNr,Resolution);

BIN_VEC=(1:size(BINs,1));
BIN_VEC=BIN_VEC(randperm(length(BIN_VEC)));

TAU_FN='DELAY_FILE';
TAU_LIST=load(TAU_FN);
N_TAU=size(TAU_LIST,1);

% Window lenght and window shift
WL=Resolution-max(TAU_LIST(:))-2*4; % one calculation per 100kbp window (-2*4 implies data points needed for numerical derivative)
WS=WL;
NrBins=size(BINs,1);ST_BIN=BINs(1,1);
FID=sprintf('%s_runDDAbash.sh',date);
unix(sprintf('touch %s',FID));
file=fopen(FID,'w');
fprintf(file,'#!/bin/bash\n\n');
fclose(file);
for w=BIN_VEC
    
    FN_DDA=sprintf('%s/Chr%d_%d-%d__%d.ascii',OD_DDA,ChrNr,BINs(w,2),BINs(w,3),BINs(w,1));
    
    LIST=[ones(1,NrBins)*BINs(w,1); BINs(1,1):BINs(end,1)];
    LIST=LIST(:,w:end);
    LIST=reshape(LIST,1,2*(NrBins-w+1));
    LIST=LIST(3:end);
    
    if exist([FN_DDA,'_CT'],'file')==0

        CMD='./run_DDA_ASCII_DNA ';
        CMD=sprintf('%s -TAU_file %s',CMD,TAU_FN); % delay file
        CMD=sprintf('%s -DATA_FN %s -OUT_FN %s',CMD,FN_ASCII,FN_DDA); % time series ascii file and output file name
        CMD=sprintf('%s -WL %d -WS %d',CMD,WL,WS); % Window lenght and window shift
        if BINs(w,1)==ST_BIN
            CMD=sprintf('%s -SELECT 1 1 0 0',CMD); % generate ST and CT DDA outputs only for first bin
        else
            CMD=sprintf('%s -SELECT 0 1 0 0',CMD); % generate CT DDA outputs only for all other bins
        end
        CMD=sprintf('%s -CT_CH_list %s',CMD,sprintf('%d ',LIST'));% which bins to CT
        disp(CMD)
        file=fopen(FID,'a');
        fprintf(file,'FN=%s\n\n',FN_DDA);
        fprintf(file,'if [[ ! -f "$FN"_CT ]]; then\n');
        fprintf(file,'touch "$FN"_CT\n');
        fprintf(file,'\t%s\n',CMD);
        fprintf(file,'fi\n\n');
        fclose(file);
    end
end
unix(sprintf('chmod u+x %s',FID));
%% generate DNA-DDA contact matrix from DDA outputs

FN_DDA_MAT=sprintf('DDA_OUT/ERGODICITY_GM12878.mat');
NrBINs=size(BINs_all,1);

if exist(FN_DDA_MAT,'file')==0
    unix(sprintf('touch %s',FN_DDA_MAT))    
    DNA_DDA=nan(NrBINs,NrBINs,N_TAU);
    
    Q_ST=load(sprintf('DDA_OUT/Chr%d_%d-%d__%d.ascii_ST',ChrNr,BINs(1,2),BINs(1,3),BINs(1,1)));
    T=Q_ST(:,1:2);WN=size(T,1);
    Q_ST=Q_ST(:,3:end);
    N_SUB_OUTST=size(Q_ST,2)/N_TAU/4;
    for w=1:size(BINs,1)-1
        Q_CT=load(sprintf('DDA_OUT/Chr%d_%d-%d__%d.ascii_CT',ChrNr,BINs(w,2),BINs(w,3),BINs(w,1)));
        T=Q_CT(:,1:2);WN=size(T,1); Q_CT=Q_CT(:,3:end);
        N_SUB_OUT=size(Q_CT,2)/N_TAU/4;

        for tau_f=1:N_TAU
            for S=w+1:size(BINs,1)
                q_ST=(reshape(Q_ST,WN,4,N_SUB_OUTST,N_TAU));
                q_ST=q_ST(1:WN,4,[w S],tau_f);
                q_ST=(nanmean(q_ST,length(size(q_ST))));%mean over subj
                
                q_CT=(reshape(Q_CT,WN,4,N_SUB_OUT,N_TAU));
                q_CT=squeeze(q_CT(1:WN,4,S-w,tau_f));
                q_CT=nanmean(q_CT,2);
                
                erg=nanmean(abs(q_ST/q_CT-1));
                DNA_DDA(BINs(w,1),BINs(S,1),tau_f)=erg;
                DNA_DDA(BINs(S,1),BINs(w,1),tau_f)=erg;
            end
       end
        w
    end
    save(FN_DDA_MAT,'DNA_DDA','-v7.3')
end

load(FN_DDA_MAT,'DNA_DDA')

tau=[3 6];% load DNA-DDA matrix for GM12878 cell line
tau_f=find(TAU_LIST(:,1)==tau(1)&TAU_LIST(:,2)==tau(2));
DNA_DDA=DNA_DDA(:,:,tau_f);

%% matrix post processing

DNA_DDA(EXCLUDE,:)=nan;
DNA_DDA(:,EXCLUDE)=nan;

%map high to low values
[~,idx_D] = sort(DNA_DDA(:), 'descend');
[~,idx_A] = sort(DNA_DDA(:), 'ascend');
idx_NAN=find(isnan(DNA_DDA));
DNA_DDA(idx_D(~ismember(idx_D,idx_NAN)))=DNA_DDA(idx_A(~ismember(idx_A,idx_NAN)));

%fill in diagonal with neighboring values
DNA_DDA(1,1)=DNA_DDA(1,2);
for d=2:size(DNA_DDA,1)
    DNA_DDA(d,d)=DNA_DDA(d-1,d);
end

%take log
DNA_DDA(DNA_DDA==0)=nan;
DNA_DDA=(log(DNA_DDA));
%% Use HiCExplorer to normalize and generate DNA-DDA-Pearson Maps
% https://hicexplorer.readthedocs.io/en/latest/index.html

DNA_DDA(isnan(DNA_DDA))=0;

BED_FILE='GM12878_GSM733772.V2.GRCh38.bed';%ChIP-Seq
BED_FILE_GENE='Homo_sapiens.GRCh38.103.bed';
HOME_TAG_DIR='HiCExplorer/TAG';
HiCExplorer_OUT=sprintf('HiCExplorer/DNA-DDA.%d.homer',Resolution);
matlabM_to_homer(BINs_all,ChrNr,HOME_TAG_DIR,HiCExplorer_OUT,DNA_DDA)

IN_h5=strrep(HiCExplorer_OUT,'homer','h5');
unix(sprintf('hicConvertFormat --matrices %s --outFileName %s --inputFormat homer --outputFormat h5'...
    ,HiCExplorer_OUT,IN_h5));
IN_norm_h5=strrep(IN_h5,'DDA','DDA.norm');
unix(sprintf('hicNormalize -m %s --normalize norm_range -o %s',...
    IN_h5,IN_norm_h5));
Pearson_h5=strrep(IN_norm_h5,'norm','norm.pearson');
unix(sprintf('hicPCA --matrix %s --whichEigenvectors "1" --method "dist_norm" --ignoreMaskedBins --extraTrack %s --pearsonMatrix %s --format bigwig --outputFileName HiCExplorer/pca.bw',...
    IN_norm_h5,BED_FILE_GENE,Pearson_h5));
unix(sprintf('hicConvertFormat --matrices %s --inputFormat h5 --outFileName %s --outputFormat ginteractions',...
    IN_norm_h5,strrep(IN_norm_h5,'.h5','')));
unix(sprintf('hicConvertFormat --matrices %s --inputFormat h5 --outFileName %s --outputFormat ginteractions',...
    Pearson_h5,strrep(Pearson_h5,'.h5','')));

tsv_IN=sprintf('%s',strrep(IN_norm_h5,'.h5','.tsv'));
csv_M_OUT=sprintf('%s',strrep(IN_norm_h5,'.h5','.csv'));
[DNA_DDA]=hiCtsv_to_MATLABcsv(Resolution,ChrSize,tsv_IN,csv_M_OUT);
DNA_DDA=load(csv_M_OUT);
DNA_DDA=DNA_DDA(BINs(:,1),BINs(:,1));

figure('units','normalized','outerposition',[0 0 1 1]);
subplot(1,2,1);imagesc(log(HiCM));axis square;colormap jet;colorbar;title('HiC contact map Chr22 GM12878');
subplot(1,2,2);imagesc(DNA_DDA);axis square;colormap jet;colorbar;title('DNA-DDA contact map Chr22 GM12878');
saveas(gcf,sprintf('%s/ContactMaps.svg',FIG_PATH))


tsv_IN=sprintf('%s',strrep(Pearson_h5,'.h5','.tsv'));
csv_M_OUT=sprintf('%s',strrep(Pearson_h5,'.h5','.csv'));
%if exist(csv_M_OUT,'file')==0,
    [DNA_DDA_P]=hiCtsv_to_MATLABcsv(Resolution,ChrSize,tsv_IN,csv_M_OUT);%end
DNA_DDA_P=load(csv_M_OUT);
DNA_DDA_P=DNA_DDA_P(BINs(:,1),BINs(:,1));

if ~isempty(EXCLUDE_PCA) % EXCLUDE_PCA 
    DNA_DDA_P(EXCLUDE_PCA(1):EXCLUDE_PCA(2),:)=[];
    DNA_DDA_P(:,EXCLUDE_PCA(1):EXCLUDE_PCA(2))=[];
end

s=unique(sort(DNA_DDA_P(:)));
figure('units','normalized','outerposition',[0 0 1 1]);clf;
subplot(1,2,1);b=imagesc(HiCP);axis square;colormap jet;colorbar;title('HiC Pearson correlation matrix GM12878');
subplot(1,2,2);b=imagesc(DNA_DDA_P,[s(2) 1]);axis square;colormap jet;colorbar;title('DNA-DDA Pearson correlation map GM12878');set(b,'AlphaData',~isnan(DNA_DDA_P))
saveas(gcf,sprintf('%s/Pearson_Matrices.svg',FIG_PATH))
%% calling A/B compartments 
HiCP(isnan(HiCP))=0;  
DNA_DDA_P(isnan(DNA_DDA_P))=0;  

EV_HiC=pca(HiCP);
EV_DDA=pca(DNA_DDA_P);
EV_DDA=movmean(EV_DDA,5);

[EV_HiC,whichPC_HiC]=Norm_PC(EV_HiC,BED_FILE,ChrNr,BINs,EXCLUDE_PCA);
[EV_DDA,whichPC_DDA]=Norm_PC(EV_DDA,BED_FILE,ChrNr,BINs,EXCLUDE_PCA);

EXCLUDE_PCA
[C,~,AUC_AB,ACC,F1] = perf_metrics(EV_DDA,EV_HiC);

sprintf('Chr%d r=%4.2f AUC=%4.2f ACC=%4.2f F1=%4.2f PC_hic=%d PC_dda=%d',ChrNr,C,AUC_AB,ACC,F1,whichPC_HiC,whichPC_DDA)

figure('units','normalized','outerposition',[-0.98 0.69 0.93 0.23]);clf;
plot(EV_HiC,'k','linewidth',1.5);axis tight;hold on;plot(EV_DDA,'m','linewidth',.8);set(gca,'FontSize',20)
xticks([1 size(HiCP,1)/2 size(HiCP,1)]);
saveas(gcf,sprintf('%s/PCs.svg',FIG_PATH))


close all


