
Iso = 'T2';
%BinStep
path = [getenv('SamakPath'),'inputs/FSD/'];

f1  = [path,'FSD_KNM2_',Iso,'.txt'];			
f2  = [path,'FSD_KNM2_',Iso,'_0p1eV.txt'];			
%f3  = [path,'FSD_KNM2_',Iso,'_Blinding.txt'];

d_orig = importdata(f1);
d_rebin = importdata(f2);
%d_blind = importdata(f3);

% GetFigure;
% plot(d_orig(:,1),d_orig(:,2));
% hold on;
% plot(d_rebin(:,1),d_rebin(:,2));
%%
E    = d_orig(:,1);
Prob = d_orig(:,2);
% clean bins with too small probability
E = E(Prob>1e-08);
Prob = Prob(Prob>1e-08);

skipLog=0;

for t = 1:5
    for i=1:numel(E)-1
        if  E(i+1)-E(i)<0.05 && skipLog==0
            % binning is finer than required -> merge bin with next one
            ENew(i)     = 0.5.*(E(i)+E(i+1));
            ProbNew(i)  = Prob(i)+Prob(i+1);
            if E(i)>5
                a = 1;
            end
            skipLog = 1; % skip next one -> already merged
        elseif skipLog==1
            % bin is skipped (becasue already merged)
            ENew(i)    = NaN;
            ProbNew(i) = NaN;
            skipLog = 0;
        elseif skipLog==0
            % binning is already larger or equal than required -> take bin as it is
            ENew(i)     = E(i);
            ProbNew(i)  = Prob(i);
            %fprintf('E = %.01f eV \n',ENew(i))
        end
    end
    ENew(isnan(ENew))=[];
    ProbNew(isnan(ProbNew))=[];
    
    E = ENew;       
    Prob = ProbNew;
    
    ENew    = NaN.*zeros(numel(E),1);
    ProbNew = NaN.*zeros(numel(E),1);
    
    fprintf('Iteration %.0f -> %.0f bins (%.0f) \n',t,numel(E),numel(d_orig(:,1)))
end



GetFigure;
plot(d_orig(:,1),d_orig(:,2));
hold on;
plot(E,Prob);

%% save
%Write2Txt('filename',[path,'FSD_KNM2_',Iso,'_0p1eV'],'nCol',2,'variable',[E';Prob'])
%%
% d = importdata(['../../inputs/FSD/FSD_KNM2_',Iso,'_0p1eV.txt']);
% close all
% plot(d(:,1),d(:,2))