function PixList = GetPixList(DataSet) 
% function to return the standard list of good pixels 
% Lisa: April 2019
PixList=1:148;
switch DataSet
    case 'FirstTritium.katrin'
%         % FPD
%         PixExclFPD = [97, 98, 110, 111, 121, 122]+1;
%         PixList=PixList(~ismember(PixList,PixExclFPD));
%         % FBM
%         PixExclFBM =[100, 112, 123]+1;
%         PixList=PixList(~ismember(PixList,PixExclFBM));
%         % Alignement
%         PixExclAli=[103:105, 113:117, 124:130, 135, 136:142, 147]+1;
%         PixList=PixList(~ismember(PixList,PixExclAli));

         PixExclFTpaper = [100,112,123:147]+1;
         PixList=PixList(~ismember(PixList,PixExclFTpaper));
    case 'Knm1' % Katrin Neutrino Mass 1
        % Save stacked pixels excluding two outer rings and pixels 98, 110, 111, 122
        % PixelList=[1:98,100:110,113:122,124];
        % https://rocketchat.kaas.kit.edu/channel/analysis?msg=kWLqnfAkm3TYLFW83
        % EGun Pixel
        PixExclEGun = []; %3+1; % (EGun)
        PixList=PixList(~ismember(PixList,PixExclEGun));
        % FPD
        PixExclFPD = [97, 98, 110, 111, 121, 122]+1; % (FPD)
        PixList=PixList(~ismember(PixList,PixExclFPD));
        % FBM
        PixExclFBM = [99, 100, 112, 123]+1;% (FBM)
        PixList=PixList(~ismember(PixList,PixExclFBM));
        % Alignement
        PixExclAli = [124   125   126   127   128   129   130   134   135   136   137   138   139   140   141   142   143   144   145   146   147]+1; % (alignment)
        PixList=PixList(~ismember(PixList,PixExclAli));  

        %%% WARNING : TECHNICAL PATH - TEMPORARY - START
        % PATCH TO REMOVE LAST PSEUDO-RING
        %PixExclOuter = 100:148;
        %PixList=PixList(~ismember(PixList,PixExclOuter));  
        % PATCH TO REMOVE 20 random pixels
        % PixList=sort(datasample(PixList,97,'Replace',false));
        %%% WARNING : TECHNICAL PATH - TEMPORARY - STOP
    case 'Knm2'
        % FPD efficiency
         PixExclFPD = [97, 98, 110, 111, 121, 122]+1; % (FPD)
         PixList=PixList(~ismember(PixList,PixExclFPD));
         % FBM shaddow
        PixExclFBM = [100]+1;% (FBM)
        PixList=PixList(~ismember(PixList,PixExclFBM)); 
       % Alignement
        PixExclAli = [112,113,123:130,134:147]+1; % (alignment)
        PixList=PixList(~ismember(PixList,PixExclAli));  
end