
Iso = 'DT';

path = [getenv('SamakPath'),'inputs/FSD/'];

f1  = [path,'FSD_KNM2_',Iso,'.txt'];			
f2  = [path,'FSD_KNM2_',Iso,'_0p5eV.txt'];			
f3  = [path,'FSD_KNM2_',Iso,'_Blinding.txt'];

d_orig = importdata(f1);
d_rebin = importdata(f2);
d_blind = importdata(f3);

GetFigure;
plot(d_orig(:,1),d_orig(:,2));
hold on;
plot(d_rebin(:,1),d_rebin(:,2));

%%
E = d_orig(:,1);
Prob = d_orig(:,2);

skipLog=0;

for t = 1:5
    for i=1:numel(E)-1
        if  E(i+1)-E(i)<0.05 && skipLog==0
            ENew(i)     = 0.5.*(E(i)+E(i+1));
            ProbNew(i)  = Prob(i)+Prob(i+1);
            skipLog = 1; % skip next one -> already merged
        elseif skipLog==1
            ENew(i)    = NaN;
            ProbNew(i) = NaN;
            skipLog = 0;
        else skipLog==0
            ENew(i)     = E(i);
            ProbNew(i)  = Prob(i);
        end
    end
    ENew(isnan(ENew))=[];
    ProbNew(isnan(ProbNew))=[];
    
    E = ENew;
    Prob = ProbNew;
    
    ENew    = NaN.*zeros(numel(E),1);
    ProbNew = NaN.*zeros(numel(E),1);
end



GetFigure;
plot(d_orig(:,1),d_orig(:,2));
hold on;
plot(E,Prob);

%% save
Write2Txt('filename',[path,'FSD_KNM2_',Iso,'_0p1eV'],'nCol',2,'variable',[E';Prob'])
%%
d = importdata(['../../inputs/FSD/FSD_KNM2_',Iso,'_0p1eV.txt']);
close all
plot(d(:,1),d(:,2))