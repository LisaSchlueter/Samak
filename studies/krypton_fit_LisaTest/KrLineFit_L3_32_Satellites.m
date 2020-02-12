function [par err chi2min ndof] = KrLineFit_L3_32_Satellites(varargin)
%
%             FIT INTEGRAL SPECTRUM
%         Fit a 83mKr Lines + Satellites
%
%          Th. Lasserre - CEA Saclay
%                August 2017
%
%              L. Schlueter - TUM
%                December 2017   

% Initialization
clear par ; clear parnomix;
addpath(genpath('../../../Samak2.0'));

% Parser
p = inputParser;
p.addParameter('fign',1,@(x)isfloat(x) && x>0);
p.addParameter('pub','YES',@(x)ismember(x,{'ON','OFF','YES'}));
p.addParameter('display','ON',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('Mode','Data',@(x)ismember(x,{'Data', 'Sim'}));
p.addParameter('CPS','OFF',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('FPD_Segmentation','OFF',@(x)ismember(x,{'OFF','RING','PIXEL'}));
p.addParameter('Pixel',1,@(x)isfloat(x) && x>0);
p.addParameter('nfit',1,@(x)isfloat(x) && x>=0);
p.addParameter('Fitter','Minuit',@(x)ismember(x,{'Minuit','Matlab'}));
p.addParameter('TD','KrL3_32_Satellites',@(x)ismember(x,{'KrL3_32','KrL3_32_HS','KrL3_32_Satellites'}));
p.addParameter('DopplerEffectFlag','ON',@(x)ismember(x,{'OFF','ON'}));
p.addParameter('SystCorr','OFF',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('HVRipples','OFF',@(x)ismember(x,{'OFF','ON'}));
p.addParameter('HVRipplesP2PV',0.52,@(x)isfloat(x) && x>0); % V or eV
p.addParameter('ConvFlag','DopplerEffect', @(x)ismember(x,{'DopplerEffect','Voigt','OFF'}))
p.addParameter('MultiPeaksFlag','ON',@(x)ismember(x,{'OFF','ON'}));

p.parse(varargin{:});

fign                   =    p.Results.fign;
pub                    =    p.Results.pub;
display                =    p.Results.display;
Mode                   =    p.Results.Mode;
CPS                    =    p.Results.CPS;
FPD_Segmentation       =    p.Results.FPD_Segmentation;
Pixel                  =    p.Results.Pixel;
nfit                   =    p.Results.nfit;
Fitter                 =    p.Results.Fitter;
TD                     =    p.Results.TD;
DopplerEffectFlag      =    p.Results.DopplerEffectFlag;
SystCorr               =    p.Results.SystCorr;
HVRipples              =    p.Results.HVRipples;
HVRipplesP2PV          =    p.Results.HVRipplesP2PV;
ConvFlag               =    p.Results.ConvFlag;
MultiPeaksFlag         =    p.Results.MultiPeaksFlag;

myOpt = {'TD', TD,'HVRipples',HVRipples,'HVRipplesP2PV',HVRipplesP2PV ...
         'FPD_Segmentation', FPD_Segmentation,'FPD_Pixel',Pixel,...
         'DopplerEffectFlag', DopplerEffectFlag, 'ConvFlag',ConvFlag,...
         'CPS',CPS, 'MultiPeaksFlag', MultiPeaksFlag};

% Parametrization: Initial (Model) Values
global A ;
A=InitKrKATRIN_LisaFit(myOpt{:});
A.ComputeKrDS(); A.ComputeKrIS();

%Initial (Model) Values 
LineE_i    = A.L3_32_E_i; 
LineW_i    = A.L3_32_W_i;
LineE_s3_i = A.L3_32_E_s3_i; 


% Init: Fit Parameters
i_E        = 0; i_W        = 0;
i_Bkg      = 0; i_Phi0     = 0;
i_N        = 0;
npar =5;      %Nmb of Fit Parameters
npar_sat = 0; %Nmb of add. Parameters for Satellite
switch MultiPeaksFlag
    case 'ON'
        %Init: Additional Fit Parameters
        i_E_s3 = 0; 
        i_Phi0_s3 = 0;
        npar_sat=2;
end
npar = 5 + npar_sat;
fit_p      = ones(npar,nfit);
fit_xi2min = ones(1,nfit);


% Read Data from File
[qU krdatamap] = readKrData('TD',TD);
qU=flip(qU);
% Background / Amplitude Initialization (Used later in KrLineModel7par_sat.m while Fitting!)
%Background
A.BKG_RateSec_i = krdatamap(Pixel,1,end);

%Amplitude
switch TD
    case {'KrL3_32_HS'} 
        qumin=100;
    case {'KrL3_32_Satellites', 'KrL3_32'}
        qumin=1;
        qu_temp=find(qU>LineE_i-5);
        qumid=qu_temp(1); %position of qU=LineE_i-5 in array 
        
end
tmp = krdatamap(Pixel,1,qumin); %take #Counts for smallest qU
tmp_mid=krdatamap(Pixel,1,qumid); %counts between main and satellite line

%%%%%%%%%%%%%%%%%%%% WORK IN PROGRESS %%%%%%%%%%%%%%%%%%%%%%%%%%
tmp_time_index_main= find(A.qU>LineE_i); %index of qU in main region
time_index_main=tmp_time_index_main(1);
tmp_time_index_sat= find(A.qU>qU(1)); %index of qU in sat region
time_index_sat=tmp_time_index_sat(1);

A.L3_32_Phi0_i  =(tmp_mid-A.BKG_RateSec_i)/(A.qUfrac(time_index_main)*A.TimeSec*0.4739); %Phi0 = Voigt_amplitude / 0.4739
A.L3_32_Phi0_s3_i = (tmp-tmp_mid)/(A.qUfrac(time_index_sat)*A.TimeSec*0.4739);%(A.p_sat3/A.p_main)*A.L3_32_Phi0_i; %  %%tmp*(A.p_sat3/A.p_main);
%fprintf('L3_32_Phi0_s3_i = %.2f \n',A.L3_32_Phi0_s3_i);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LinePhi0_i = A.L3_32_Phi0_i;
LinePhi0_s3_i = A.L3_32_Phi0_s3_i;
fprintf('Init: L3_32_Phi0_i = %.2f \n',A.L3_32_Phi0_i);
fprintf('Init: L3_32_Phi0_i_sat = %.2f \n',A.L3_32_Phi0_s3_i);
fprintf('Compare with probabilities:  Phi0_sat= %.2f \n', (A.p_sat3/A.p_main)*A.L3_32_Phi0_i);


switch Fitter
    case 'Minuit'
    case 'Matlab'
        %             tbdisf = @(mb,qb,bb,nb,qu) A.ComputeTBDISf(qu,mb,qb,bb,nb,0,0);
        %             tmp = @(mb,qb,bb,nb,e) interp1(A.qU,tbdisf(mb,qb,bb,nb,A.qU),e);
        %             myf = fittype(@(mb,qb,bb,nb,qu) tmp(mb,qb,bb,nb,qu),...
        %                 'independent','qu','coefficients', {'mb','qb','bb','nb'});
        %             opts = fitoptions('StartPoint',[0,0,0,0],...
        %                 'Method', 'NonlinearLeastSquares','Display','off',...
        %                 'Weights',1./A.TBDIS);
end


% Number of Fits
tic
progressbar('Krypton Line Fits');

for f=1:1:nfit
    progressbar(f/nfit);
    switch Mode
        case 'Data'
            Count    = (krdatamap(Pixel,1,qumin:end));
            CountErr = (krdatamap(Pixel,2,qumin:end));
            Data = [A.qU(:),Count(:),CountErr(:)];
            %disp(Data);
        case 'Sim'
            A.ComputeKrDS(); A.ComputeKrIS();
            AddStatFluctKrIS(A);
            Data = [A.qU,A.KrIS,A.KrISE];
    end
    
    
    switch Fitter
        case 'Minuit'
            fprintf('Minuit fitter\n');
            % Initializing Fit
            switch MultiPeaksFlag
                case 'ON'
                    fprintf('MultiPeaksFlag ON\n');
                    ParIni = [i_E i_W i_Phi0 i_Bkg i_N i_E_s3 i_Phi0_s3];
                    parnames = ['E W Phi0 Bkg N E_s3 Phi0_s3'];
                case 'OFF'
                    fprintf('MultiPeaksFlag OFF\n');
                    ParIni = [i_E i_W i_Phi0 i_Bkg i_N];
                    parnames = ['E W Phi0 Bkg N'];
            end
            
            tmparg = sprintf(['set pri -10 ; fix  5'...
                'set now; min ; imp' ]);
            Args = {ParIni, Data, '-c',tmparg};
            
            switch MultiPeaksFlag
                case 'OFF'
                    [par, err, chi2min, errmat] = fminuit('KrLineChi2',Args{:});
                case 'ON'
                    [par, err, chi2min, errmat] = fminuit('KrLineChi2_sat',Args{:});
            end
            % ndof = size(A.qU)-npar; %new Lisa
            
            fit_p(:,f)=par;  fit_xi2min(f) = chi2min;
            
        case 'Matlab'
            %                 [fit1,gof,fitinfo] = fit(A.qU,A.TBDIS,myf,opts);
            %                 par = [fit1.mb fit1.qb fit1.bb fit1.nb 0 0];
            %                 fit_p(:,f)=par;
            %                 err = [0 0 0 0 0 0];
            %                 chi2min = gof.sse; fit_xi2min(f) = chi2min;
    end
    
    
    switch display
        case 'ON'
            fprintf(2,'--------------------------------------------------------------\n');
            fprintf('  Processing \t= %g %% \n',f/nfit*100);
            fprintf(2,'--------------------------------------------------------------\n');
            fprintf(2,'  E\t\t= %.2f \t +- \t %.2f \t eV \n', (LineE_i+par(1)),err(1));
            fprintf(2,'  W\t\t= %.2f \t +- \t %.2f \t meV \n', (LineW_i+par(2))*1e3,err(2)*1e3);
            fprintf(2,'  Phi0\t\t= %.2f \t +- \t %.2f \t    \n', (A.L3_32_Phi0_i+par(3)),err(3));
            fprintf(2,'  Bkg\t\t= %.2f \t +- \t %.2f \t cps \n', (A.BKG_RateSec_i+par(4)),err(4));
            fprintf(2,'  N\t\t= %.2f \t +- \t %.2f \t    \n', par(5),err(5))
            switch MultiPeaksFlag
                case 'ON'
                    fprintf(2,'  E(sat)\t= %.2f \t +- \t %.2f \t eV \n', (LineE_s3_i+par(6)), err(6));
                    fprintf(2,'  Phi0(sat)\t= %.2f \t +- \t %.2f \t eV \n', (A.L3_32_Phi0_s3_i+par(7)), err(7));
            end
            ndof = A.nqU-4;
            fprintf(2,'  Chi2\t\t= %g / %g dof \n',chi2min,ndof);
            fprintf(2,'--------------------------------------------------------------\n');
    end
  toc
end
