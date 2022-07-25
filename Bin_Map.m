function [BINS] = Bin_Map(ChrSize,Resolution)
% Resolution in kB

Nr_Bins=ceil(ChrSize/(Resolution));

%BINS=1:Nr_Bins;
BINS=nan(Nr_Bins,2);
BINS(1,:)=[0, Resolution];
for k=2:Nr_Bins-1
    BINS(k,:)=BINS(k-1,:)+Resolution;
end
BINS(end,:)=[BINS(end-1,2)+1,ChrSize];

BINS=[transpose(1:Nr_Bins),BINS];
BINS(:,2)=BINS(:,2)+1;
BINS(:,3)=BINS(:,3)+1;

% start=1:Resolution*1e3:ChrSize;start(start>ChrSize)=[];
% ende=[Resolution*1e3:Resolution*1e3:ChrSize];
% ende(end)=ChrSize;
% BINS=[BINS',start',ende'];



end

