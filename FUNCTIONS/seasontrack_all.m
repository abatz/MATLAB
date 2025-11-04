function [season]=seasontrack(tmx);

% Wotton and Flannigan 1993 say you need 3 conseq days of tmx>12C to start
% season, and 3 conseq days of tmx<5C to end season
season=single(zeros(size(tmx)));lastseason=zeros(size(tmx,1),1);
on=ones(size(tmx,1),1);days12=on*0;days5=on*0;
for yr=1:size(tmx,3);
    for doy=1:size(tmx,2);
        foff=find(lastseason==0);fon=find(lastseason==1);
        finc_hot=find(tmx(foff,doy,yr)>=12);days12(foff(finc_hot))=days12(foff(finc_hot))+1;
        finc_cold=find(tmx(fon,doy,yr)<=5);days5(fon(finc_cold))=days5(fon(finc_cold))+1;
        f1=find(days12==3);season(f1,doy,yr)=1;days12(f1)=0;on(f1)=1;
        f1=find(days5==3);season(f1,doy,yr)=0;days5(f1)=0;on(f1)=0;
        season(:,doy,yr)=on;lastseason=season(:,doy,yr);
    end
end
%season=ones(size(season));

