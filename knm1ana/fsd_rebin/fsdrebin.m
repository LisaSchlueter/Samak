%
% Rebin FSD Files
%
extra = '_FullRange';
dataset = 'KNM1';
isotopologue='T2';

file = importdata([getenv('SamakPath'),sprintf('/inputs/FSD/FSD_%s_%s_Doppler%s.txt',dataset,isotopologue,extra)]);

% clear
clear ee; clear eb;

% Read Files
nbin = numel(file(:,1));

% Loop on original bins
counter = 1;
for i = 1:1:nbin-1
    %fprintf(2,'Bin Original %.0f : %.3f %.3g \n',i,file(i,1),file(i,2));
    ee(counter) = 0;
    ep(counter) = 0;
    if i<70
        for j = i:1:i
            ee(counter) = ee(counter) + file(j,1);
            ep(counter) = ep(counter) + file(j,2);
        end
        ee(counter) = ee(counter);
        ep(counter) = ep(counter);
        counter     = counter+1;
         
    else
        if ~rem(i,5)
            disp(i)
            for j = i:1:i+4
                %fprintf('Bin Original %.0f : %.3f %.3g \n',j,file(j,1),file(j,2));
                ee(counter) = ee(counter) + file(j,1);
                ep(counter) = ep(counter) + file(j,2);
            end
            ee(counter) = ee(counter)/5;
            ep(counter) = ep(counter);
            counter = counter+1;
            fprintf('Bin New %.0f : %.3f %.3g \n',counter-1,ee(counter-1),ep(counter-1));
        end        
    end
    
end

figure(1)
plot(file(:,1),file(:,2));
hold on
plot(ee,ep);
hold off

variable = [ee' ep'];
  
switch isotopologue 
    case 'T2'
save([getenv('SamakPath'),'/inputs/FSD/FSD_KNM1_T2_Doppler0p5eV',extra,'.txt'],'variable','-ascii');
    case 'HT'
save([getenv('SamakPath'),'/inputs/FSD/FSD_KNM1_HT_Doppler0p5eV',extra,'.txt'],'variable','-ascii');
    case 'DT'
save([getenv('SamakPath'),'/inputs/FSD/FSD_KNM1_DT_Doppler0p5eV',extra,'.txt'],'variable','-ascii');
end
