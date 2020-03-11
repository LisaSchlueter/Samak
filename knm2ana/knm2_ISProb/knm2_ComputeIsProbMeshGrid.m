


p = inputParser;
p.addParameter('WGTS_','',@(x)isfloat(x));
p.addParameter('MACE_Bmax_T','',@(x)isfloat(x));
p.addParameter('WGTS_B_T','',@(x)isfloat(x));

MACE_Bmax_T = p.Results.MACE_Bmax_T;
WGTS_B_T    = p.Results.WGTS_B_T;
%%
A = ref_FakeRun_KNM2_CD84_2hours;

%%
Eiscs = 18575+(-1000:10:1000);
myTheta =  asin(sqrt(A.WGTS_B_T./AMACE_Bmax_T));

Theta = @(bmax)  asin(sqrt(A.WGTS_B_T./bmax));

Bmax_Min = A.MACE_Bmax_T.*0.90;
Bmax_Max = A.MACE_Bmax_T.*1.10;

BmaxSamples = linspace(Bmax_Min,Bmax_Max,20);
for i=1:numel(BmaxSamples)
    progressbar(i/numel(BmaxSamples));
    A.MACE_Bmax_T = BmaxSamples(i);
    A.ComputeISProb('Energy',reshape(Eiscs,[1,1,numel(Eiscs)]));
end

