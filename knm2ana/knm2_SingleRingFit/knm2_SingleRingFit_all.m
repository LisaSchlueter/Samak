%  script to fit all single rings for
% - different ranges
% - different RW periods
% with and without neutrino mass
RunLists = {'KNM2_RW1','KNM2_RW2','KNM2_RW3','KNM2_Prompt'};
freePars = {'E0 Bkg Norm','mNu E0 Bkg Norm'};
Ranges = [40,90];
ROIFlag = '14keV';
MosCorrFlag = 'ON';
RecomputeFlag = 'ON';

InputArg = [...
    {'RunList',RunLists{1},'freePar',freePars{1},'Range',Ranges(1)};... % mnu fix - 40 eV
    {'RunList',RunLists{2},'freePar',freePars{1},'Range',Ranges(1)};...
    {'RunList',RunLists{3},'freePar',freePars{1},'Range',Ranges(1)};...
    ...
    {'RunList',RunLists{1},'freePar',freePars{1},'Range',Ranges(2)};... % mnu fix - 90 eV
    {'RunList',RunLists{2},'freePar',freePars{1},'Range',Ranges(2)};...
    {'RunList',RunLists{3},'freePar',freePars{1},'Range',Ranges(2)};...
    ...
    {'RunList',RunLists{1},'freePar',freePars{2},'Range',Ranges(1)};...% mnu free - 40 eV
    {'RunList',RunLists{2},'freePar',freePars{2},'Range',Ranges(1)};...
    {'RunList',RunLists{3},'freePar',freePars{2},'Range',Ranges(1)};...
    ...
    {'RunList',RunLists{1},'freePar',freePars{2},'Range',Ranges(2)};... % mnu free - 90 eV
    {'RunList',RunLists{2},'freePar',freePars{2},'Range',Ranges(2)};...
    {'RunList',RunLists{3},'freePar',freePars{2},'Range',Ranges(2)};...
    ...
    ];

parfor i=1:size(InputArg,1)
    knm2_SingleRingFit(InputArg{i,:},'ROIFlag',ROIFlag,'MosCorrFlag',MosCorrFlag,RecomputeFlag,RecomputeFlag)
end