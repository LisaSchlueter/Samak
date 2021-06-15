
Bconst = 0.221; %cps

alpha = 3*1e-06;

qUfrac= [0.00263112088588541,0.00263224146547905,0.00263149441241661,0.00263522967772878,0.00264307373488432,0.00263336204507270,0.00263186793894784,0.00263747083691607,0.00263522967772878,0.00263448262466635,0.00452788861140019,0.00452377981955681,0.00451892397465099,0.00452452687261923,0.00452415334608802,0.00452602097874410,0.00452004455424465,0.00546133141290893,0.00707907481960535,0.00938111883148949,0.0126128703795702,0.0173308839953593,0.0242056398023895,0.0341862687164807,0.0579201445099445,0.0676295931623727,0.0773323183372393,0.0870391523039491,0.0678914352607550,0.0482073341187350,0.0482032253268917,0.0482047194330165,0.0482062135391415,0.0481740902574570,0.0482099488044536,0.0482133105432345,0.0482054664860791,0.0482039723799543];
qUfrac = qUfrac(26:29);
TimeSec = 7416;
TimeSecSubrun = qUfrac.*TimeSec;
TimeSecSubrun_CumSum = cumsum(TimeSecSubrun);

BkgPtFunc = @(t) alpha.*t+Bconst;
Time_global = 1:1:max(TimeSecSubrun_CumSum);
Rate_global = zeros(numel(Time_global),1);

Time_local = 0;
SubRun_local = 1;
for i=1:numel(Time_global)
    
    Rate_global(i) = BkgPtFunc(Time_local);
    
    if TimeSecSubrun_CumSum(SubRun_local)<Time_global(i)
        Time_local = 0;
        SubRun_local = SubRun_local+1;
    else
        Time_local = Time_local+1;
    end  
end

%%
GetFigure;
plot(Time_global/60,Rate_global.*1e3,'LineWidth',1.5)
xlabel('Measurement time (min)');
ylabel('Background rate (mcps)');
xlim([0,max(Time_global)/60]);
PrettyFigureFormat;

pltname = [getenv('SamakPath'),'knm2ana/knm2_PngBkg/plots/PtBkgDummyPlot.png'];
print(pltname,'-dpng','-r300');


