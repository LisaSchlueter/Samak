SanityPlot = 'OFF';
ZoomPlot = 'ON';
RecomputeFlag = 'ON';
E0 = 18574;
mnu = 0;

Sigma = sqrt(0.0935.^2 + 0.0001.^2);

FSDdir = [getenv('SamakPath'),'inputs/FSD/'];
ttfsdfilename = [FSDdir,'FSD_KNM2_T2_Blinding.txt'];
%Rebinning
ttfsdfile_temp=importdata(ttfsdfilename);
ttfsdfile = ttfsdfile_temp(ttfsdfile_temp(:,2)>(10^-8),[1 2]);
%
[obj.TTexE, TTexE_index] = sort(ttfsdfile(:,1)); %sort from small to large excitation energies           
obj.TTexE = obj.TTexE';
obj.TTexP = (ttfsdfile(TTexE_index,2))';

plot(obj.TTexE,obj.TTexP);
hold on;

[obj.TTexE,obj.TTexP] = FSD_Convfun_relic(obj.TTexE,obj.TTexP,...
                    squeeze(Sigma(1,:,:)),...
                    'SanityPlot',SanityPlot,'ZoomPlot',ZoomPlot,...
                    'filename',ttfsdfilename,'RecomputeFlag',RecomputeFlag);
                
plot(obj.TTexE,obj.TTexP);
hold off;