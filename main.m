%% load variables
load('vars.mat')

%% generate **1D DNA walk** and save to file <FN_ASCII>

makeTS=0;
ASCII_DIR=('.');
FN_ASCII=(sprintf('%s/RandomWalk_Chr1.ascii',ASCII_DIR));
% vergleiche ascii file
if makeTS==1
    DATA=load(sprintf('../../../data/PERL_SEQ_ALL/Chr%d_BINALL.perl.ascii',ChrNr));
    unix(sprintf('touch %s',FN_ASCII))
    TS=nan(Resolution,50);
    for w_x=1:size(BINs_ALL,1)
        X=DATA(BINs_ALL(w_x,2):(BINs_ALL(w_x,3)-1));
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
    title(sprintf('Random Walk DNA chr1:%d-%d',BINs(1,2),BINs(end,3)),'fontsize',20,'interpreter','latex')
end

%% generate bash script to run DDA
OD_DDA=sprintf('DDA_OUT',ChrNr,Resolution);

BIN_VEC=(1:size(BINs,1));
BIN_VEC=BIN_VEC(randperm(length(BIN_VEC)));

TAU_FN='DELAY_FILE';
TAU_LIST=load(TAU_FN);N_TAU=size(TAU_LIST,1);

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
        fprintf(file,'\ttouch "$FN"_CT\n');
        fprintf(file,'\t%s\n',CMD);
        fprintf(file,'fi\n\n');
        fclose(file);
    end
end
unix(sprintf('chmod u+x %s',FID));
%% generate DNA-DDA contact matrix from DDA outputs

FN_DDA_MAT=sprintf('DDA_OUT/ERGODICITY.mat');

if exist(FN_DDA_MAT,'file')==0
    unix(sprintf('touch %s',FN_DDA_MAT))
    DNA_DDA=nan(NrBins,NrBins,N_TAU);
    
    Q_ST=load(sprintf('DDA_OUT/Chr%d_%d-%d__%d.ascii_ST',ChrNr,BINs(1,2),BINs(1,3),BINs(1,1)));
    T=Q_ST(:,1:2);WN=size(T,1);
    Q_ST=Q_ST(:,3:end);
    N_SUB_OUTST=size(Q_ST,2)/N_TAU/4;
    for w=1:size(BINs,1)-1
        Q_CT=load(sprintf('DDA_OUT/Chr%d_%d-%d__%d.ascii_CT',ChrNr,BINs(w,2),BINs(w,3),BINs(w,1)));
        T=Q_CT(:,1:2);WN=size(T,1); Q_CT=Q_CT(:,3:end);
        N_SUB_OUT=size(Q_CT,2)/N_TAU/4;
        
        for tt=1:N_TAU
            for S=w+1:size(BINs,1)
                q_ST=(reshape(Q_ST,WN,4,N_SUB_OUTST,N_TAU));
                q_ST=q_ST(1:WN,4,[w S],tt);
                q_ST=(nanmean(q_ST,length(size(q_ST))));%mean over subj
                
                q_CT=(reshape(Q_CT,WN,4,N_SUB_OUT,N_TAU));
                q_CT=squeeze(q_CT(1:WN,4,S-w,tt));
                q_CT=nanmean(q_CT,2);
                
                erg=nanmean(abs(q_ST/q_CT-1));
                DNA_DDA(BINs(w,1),BINs(S,1),tt)=erg;
                DNA_DDA(BINs(S,1),BINs(w,1),tt)=erg;
            end
        end
        w
    end
    save(FN_DDA_MAT,'DNA_DDA','-v7.3')
end

load(FN_DDA_MAT,'ERGODICITY')
DNA_DDA=ERGODICITY(First_NonZero_Bin:end,First_NonZero_Bin:end,:);
tau=[3 6];% load DNA-DDA matrix for GM12878 cell line
tau=find(TAU_LIST(:,1)==tau(1)&TAU_LIST(:,2)==tau(2));

DNA_DDA=DNA_DDA(:,:,tau);

%% matrix post processing

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
DNA_DDA(isnan(DNA_DDA))=0;


figure;subplot(1,2,1);imagesc(DNA_DDA);axis square;colormap jet;colorbar;title('DNA-DDA contact map GM12878')
subplot(1,2,2);imagesc(log(HiCM(BINs(1,1):BINs(end,1),BINs(1,1):BINs(end,1))));axis square;colormap jet;colorbar;title('HiC contact map GM12878')



