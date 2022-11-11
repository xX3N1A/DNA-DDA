function [DENSITY] = H3K4me1Density(BINs,Chr,BED_FILE)
   
DENSITY=zeros(size(BINs,1),2);
DENSITY(:,1)=BINs(:,1);
fid = fopen(BED_FILE);
tline = fgetl(fid);
while ischar(tline)
  %  disp(tline)
    if contains(tline,sprintf('chr%d',Chr))
       B = regexp(tline,'(\t)(\d*)(\t)(\d*)(\t)','tokens');
       B=B{1};
       A=str2num(B{[2]});
       E=str2num(B{[4]});clear B
       f=find(BINs(:,2)<=A&BINs(:,3)>=E);
       if ~isempty(f)
           DENSITY(f,2)=DENSITY(f,2)+1;
       end
    end
    
    tline = fgetl(fid);
end  

end