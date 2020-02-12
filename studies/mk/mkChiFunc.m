%%Class that defines some functions for the senstivity study
%
% Marc Korzeczek - June, 2017


classdef mkChiFunc
    
    %% Minimizer functions and models. - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    methods(Static)
        %Model - Minuit
        function model = calcModel(p)
            global Theory;
            Theory.ComputeTBDDS('m_bias',p(1),'Q_bias',p(2),'B_bias',p(3),'mnu4Sq_Bias',p(5),...
                'sin2T4_Bias',p(6)); Theory.ComputeTBDIS();
            model = ((1+p(4)).*Theory.TBDIS);
        end
        
        %Minimizer - Minuit
        function chi2 = calcChi2(p,X)
            % p = vector of nuisance parameters
            % X = ['vector of bin centers';'vector of data';'vector of uncertainties';'info']
            nn = X(:,2)~=0; 
            Data = X(nn,2); 
            Sig = X(nn,3);
            Model = mkChiFunc.calcModel(p);
            chi2 = real( sum( (Data - Model).^2 ./ Sig.^2 ) );
            
            %Additional chi2-term from external CeSOX simulation
            global CeSox;
            if(~isempty(CeSox))
                if(~strcmp(CeSox.Type,'OFF'))
                    sinSqTh4 = abs(p(6)); sinSq2Th4 = abs(1-(1-2*sinSqTh4)^2);  mnu4Sq_eV = abs(p(5));
                    chi2 = chi2 + CeSox.chi2f(sinSq2Th4,mnu4Sq_eV);%+(mnu4Sq_eV>100)*1e3 +(sinSqTh4>.5)*1e3;
                end
            end
        end
        
        %Minimizer - Minuit w. Covariance Matrix
        function chi2 = calcChi2CovMat(p,X)
            nn = X(:,2)~=0;  Data = X(nn,2);  Model = mkChiFunc.calcModel(p);
            
            global Theory;
            chi2 = (Data - Model)' * ( Theory.fitcov \ (Data - Model) );
            
            global CeSox;
            if(~isempty(CeSox))
                if(~strcmp(CeSox.Type,'OFF'))
                    sinSqTh4 = abs(p(6)); sinSq2Th4 = abs(1-(1-2*sinSqTh4)^2);  mnu4Sq_eV = abs(p(5));
                    chi2 = chi2 + CeSox.chi2f(sinSq2Th4,mnu4Sq_eV);%+(mnu4Sq_eV>100)*1e3 +(sinSqTh4>.5)*1e3;
                end
            end
        end
    end
    
    
    %% General mimization function for Matlab and Minuit fitter
    methods(Static)
        function [chi2,par,err,TBDData] = minimizeTheoryToData(DataOpt, TheoryOpt, varargin )
            %SpecOpt = {'TD','DR30','TimeYear',5,'nTeBinningFactor',20, ...
            %           'BKG_RateAllFPDSec',10e-3,'TTFSD','SAENZ'};
            %DataOpt =[Specopt,{...}];
            %TheoryOpt =[Specopt,{...}];
            %varargin: (Property,Value)-pairs (s. below)
            
            %Parse funciton specific options
            p = inputParser;
            p.addParameter('Fitter','Minuit',@(x)ismember(x,{'Minuit','Matlab'}));
            p.addParameter('CovMat','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('CeSoxType','OFF',@(x)ismember(x,{'OFF','5','50','100','500'}));
            p.addParameter('PlotsShow','OFF',@(x)ismember(x,{'ON','OFF'}));
            
            %                            mSqB QB    RbkgB NormB m4SqB sinSqTh4B
            p.addParameter('FitParInit',[0    0     0     0     0     0]   ,@(x)isa(x,'numeric')&& length(x)==6);
            p.addParameter('FitParFix', [true false false false true  true],@(x)isa(x,'logical')&& length(x)==6);
            
            p.addParameter('FitParUseLimits','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('FitParLower',[-10  -5     -10e-3 -1   -500    0]   ,@(x)isa(x,'numeric')&& length(x)==6);
            p.addParameter('FitParUpper',[10   5      1      1    500   .5]   ,@(x)isa(x,'numeric')&& length(x)==6);
            p.parse(varargin{:});
            
            Fitter=p.Results.Fitter;
            CovMat=p.Results.CovMat;
            CeSoxType = p.Results.CeSoxType;
            PlotsShow = p.Results.PlotsShow;
            FitParInit= p.Results.FitParInit;
            FitParFix= p.Results.FitParFix;
            FitParUseLimits= p.Results.FitParUseLimits;
            FitParLower= p.Results.FitParLower;
            FitParUpper= p.Results.FitParUpper;

            
            addpath('fminuit'); addpath('tools');  global Theory;
            
            
            global CeSox;
            CeSox.Type = CeSoxType;
            
            % Read CeSox data
            if(~strcmp(CeSox.Type,'OFF'))
                CeSox.fname = ['data/CeSox/CeSOXsmall_100x100_chi2map', CeSox.Type ,'days.mat'];
                if exist(CeSox.fname,'file')==2
                    t=load(CeSox.fname,'dm2','sin22theta','chi2map');
                    CeSox.c2map = t.chi2map;   CeSox.s22 = t.sin22theta;   CeSox.dm2 = t.dm2;
                    CeSox.chi2f = @(x,y) interp2(t.sin22theta,t.dm2,t.chi2map(:,:,1),x,y,'linear',0);
                    % returns 0 if (x,y) out of bounds (dm(1)<x<dm(end), sin22theta(1)<y<sin22theta(end))
                    clear t;
                end
            end
            
            
            %define theoretic spectrum
            Theory=InitKATRIN(TheoryOpt{:});
            
            switch Fitter
                case 'Minuit'
                    %Build covariance matrix
                    switch CovMat
                        case 'ON'
                            nCMTrials = 100;
                            fName = sprintf('./data/CovMat/WGTSMACE_CovMat_%gTrials_%s_Year_%g.mat',...
                                nCMTrials,Theory.TD,Theory.TimeYear);
                            if exist(fName, 'file') == 2
                                CM=importdata(fName);
                            else
                                Theory.ComputeTBDISCovarianceMatrix('nTrials',nCMTrials);
                                CM=importdata(fName);
                            end
                    end
                case 'Matlab'
                    %define theory
                    tmp = @(mSqb,qb,bb,nb,m4Sqb,s2T4b,qu) Theory.ComputeTBDIS4f(qu,mSqb,qb,bb,nb,0,0,m4Sqb,s2T4b);
                    TBDISLIf = @(mSqb,qb,bb,nb,m4Sqb,s2T4b,e) interp1(Theory.qU,tmp(mSqb,qb,bb,nb,m4Sqb,s2T4b,Theory.qU),e);
                    
                    Theoryf = fittype(@(mSqf,qf,bf,nf,m4Sqf,sSqT4f,qu) TBDISLIf(mSqf,qf,bf,nf,m4Sqf,sSqT4f,qu),...
                        'independent','qu','coefficients', {'mSqf','qf','bf','nf','m4Sqf','sSqT4f'});
            end
            
            % Init parameters
            npar = 6;  par = ones(1,npar);  err = zeros(1,npar);   chi2 = 0;
            
            %define measured data
            TBDData=InitKATRIN(DataOpt{:});
            Data = [TBDData.qU,TBDData.TBDIS,TBDData.TBDISE];
            % Set Fitter dependent fit procedure
            switch Fitter
                case 'Minuit'
                    %Fit Options
                    %                mSqB   QB    RbkgB  NormB  m4SqB  sinSqTh4B
                    %ParIni =   [      0      0      0      0      0        0     ];
                    %ParFix = 'fix     1                          5        6     ;';
                    %ParLim = 'lim 1 %f %f                                        '
                    ParNum = 1:npar; 
                    ParFix=['fix ',num2str(ParNum(FitParFix))];
                    ParLim=''; if(strcmp(FitParUseLimits,'ON')); ParLim = reshape([repmat('; set lim ',npar,1),num2str([ParNum;FitParLower;FitParUpper]')]',1,[]); end
                    Args = {FitParInit, Data, '-c ',[ParFix,ParLim,'; set pri -10 ; set now; min']};
                    %Start fit
                    switch CovMat
                        case 'OFF'
                            [ fpar, ferr, chi2min, ~ ] = fminuit('mkChiFunc.calcChi2',Args{:});
                        case 'ON'
                            switch SystCorr
                                case 'OFF'
                                    Theory.fitcov = diag(diag(CM.WGTSMACE_CovMat)) + diag(Theory.TritSpecData.TBDIS);
                                case 'ON'
                                    Theory.fitcov = (CM.WGTSMACE_CovMat) + diag(Theory.TritSpecData.TBDIS);
                            end
                            [ fpar, ferr, chi2min, ~ ] = fminuit('mkChiFunc.calcChi2CovMat',Args{:});
                    end
                    par(:)=fpar; err(:)=ferr; chi2=chi2min;
                    
                case 'Matlab'
                    %Fit options                    
                    ParIni = {'StartPoint',FitParInit};
                    FitParLower(FitParFix)=0; FitParUpper(FitParFix)=0; %fix variables; ParFix defined by ParLim
                    ParLim = {'lower',FitParLower,'upper',FitParUpper}; %assign (non-loose) limits
                    
                    opts = fitoptions('Display','off','Robust','Off',ParIni{:},ParLim{:},...
                        'Method', 'NonlinearLeastSquares','Weights',1./Data(:,2));
                    %Start fitting
                    [fobj,gof,output] = fit(Data(:,1),Data(:,2),Theoryf,opts);
                    par(:)=[fobj.mSqf fobj.qf fobj.bf fobj.nf fobj.m4Sqf fobj.sSqT4f];
                    chi2=gof.sse;  out=output;
            end
            %print result
            fprintf(2,'--------------------------------------------------------------\n');
            fprintf(2,'  dmnuSq   \t= %g ±\t%g\tmeV^2 \n', par(1)*1e6, err(1)*1e6);
            fprintf(2,'  dQ       \t= %g ±\t%g\tmeV \n', par(2)*1e3, err(2)*1e3);
            fprintf(2,'  dB       \t= %g ±\t%g\tmcps \n',par(3)*1e3, err(3)*1e3);
            fprintf(2,'  dN       \t= %g ±\t%g\n',       par(4), err(4));
            fprintf(2,'  dmnu4Sq  \t= %g ±\t%g\tmeV^2\n',  par(5)*1e6, err(5)*1e6);
            fprintf(2,'  dsinSqTh4\t= %g ±\t%g\n',       par(6), err(6));
            fprintf(2,'  Chi2 \t= %g / %g dof \n',       chi2,TBDData.nqU-4);
            fprintf(2,'--------------------------------------------------------------\n');
            
            %% Standard plots
            if(strcmp(PlotsShow,'ON'))
                mkChiFunc.PlotTheoVsData(DataOpt,TheoryOpt,par,9);
                %TBDData.PlotTD('fign',10);
                %TBDData.PlotTBDDS('fign',11);
                %TBDData.PlotTBDIS('fign',12);
                %TBDData.PlotKTF('fign',13);
                %TBDData.PlotFSD('fign',14,'TT','ON');
                %TBDData.PlotPIS('fign',15);
            end
        end
        
        
        
        function      PlotTheoVsData(DataOpt,TheoryOpt,par,fign)
            if(~exist('fign','var')); fign=20; end
            
            TBDData=InitKATRIN(DataOpt{:});
            x=TBDData.qU/1e3; yData=TBDData.TBDIS; err=TBDData.TBDISE;
            
            global Theory; Theory=InitKATRIN(TheoryOpt{:});
            yTheo=mkChiFunc.calcModel(par(:));
            
            %switch Fitter
            %    case 'Minuit'
            %        gTheory = plot(x,mkChiFunc.calcModel(par(:))); hold off;
            %    case 'Matlab'
            %        plot(x,Theoryf(par(:),Data(:,3))); hold off;
            %end
            
            %% Plot last theory vs data
            figure(fign);
            subplot(2,1,1);
            errorbar(x,yData,err*100,'k+'); hold on;
            plot(x,yTheo,'b:','LineWidth',2); hold off;
            
            title( num2str( par(:)' ) );
            ax=gca; set(ax,'YScale','log');  legend('data (error x100)','theory','location','southwest');
            ylabel('integral spectrum'); %xlabel('electron kinetic energy in keV');
            
            subplot(2,1,2);
            errorbar(x,yData./yData,err./yData,'k+'); hold on;
            plot(x,yTheo./yData,'b:','LineWidth',2); hold off;
            
            ax=gca; set(ax,'YScale','log');  legend('data/data','theory/data','location','southwest');
            xlabel('electron kinetic energy in keV'); ylabel('ratio');
            
            FigureFormat;
            
        end
        
    end
    %% General function reimplementations to avoid using Matlabs statistics and machine learning
    %toolbox. - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    methods(Static)
        
        function out = myGaus(x,mean,sig)
            if ~isa(x,'double'); error('1st argument should be of type double.'); end
            if ~isa(mean,'double'); error('2nd argument should be of type double.'); end
            if ~isa(sig,'double'); error('3rd argument should be of type double.'); end
            if sig<=0; error('Choose Sigma >0'); end
            
            out = 1/sig/sqrt(2*pi)*exp(-(x-mean).^2./sig^2/2);
        end
        
        function out = myGausCdf(x,mean,sig)
            if ~isa(mean,'double'); error('2nd argument should be of type double.'); end
            if ~isa(sig,'double'); error('3rd argument should be of type double.'); end
            if sig<=0; error('Choose Sigma >0'); end
            
            xmin=mean-10*sig;
            f=@(y) mkChiFunc.myGaus(y,mean,sig);
            ff=@(xmax) integral(f,xmin,xmax);
            out = arrayfun(ff,x);
        end
        
        function out = myChi2(x, df)
            addpath('fminuit');
            out=1/2.^(df/2)./gammac(df/2).*x.^(df/2-1).*exp(-x/2);
        end
        
        function out = myChi2Cdf(x, df)
            addpath('fminuit');
            f=@(y) mkChiFunc.myChi2(y,df);
            ff=@(y) integral(f,0,y);
            out = arrayfun(ff,x);
        end
        
        function out = myInvChi2CdfRec(prob, df)
            %'Recursive' implementation for chi2inv() from matlabs statistics
            %and machine learning toolbox. Rather slow but safe.
            %out: inverse chi2 cdf at p (0..inf)
            %prob: probability of invesre chi2 cdf (0..1)
            %df: degrees of freedom (0,1,2,..)
            addpath('fminuit');
            
            func=@(A) calcInvChi2Cdf(A,df);
            out=arrayfun(func,prob);
            
            %Actual recursive function
            function out = calcInvChi2Cdf(probVal,df)
                if(probVal>1); out=nan; return; end
                if(probVal>1-1e-12); out=inf; return; end
                if(probVal<1e-12); out=0; return; end
                outt=0; pt=0; step=5; res=0; acc=1e-5*probVal*(1-probVal);
                out = calcInvChi2CdfRec(probVal,df);
                
                function out = calcInvChi2CdfRec(p,df)
                    diff=p-pt;
                    if( diff > acc )
                        outt=outt+step;
                        pt=mkChiFunc.myChi2Cdf(outt,df); %caching calculated cdf
                        calcInvChi2CdfRec(p,df);
                    elseif( diff < -acc )
                        step=step/2; outt=outt-step;
                        pt=mkChiFunc.myChi2Cdf(outt,df); %caching calculated cdf
                        calcInvChi2CdfRec(p,df);
                    else
                        res=outt;
                    end
                    out=res;
                end
            end
        end
        
        function out = myInvChi2Cdf(prob, df)
            %Quick and dirty implementation of  chi2inv().
            %out: inverse chi2 cdf at p
            %prob: chi2 cdf at out
            %df: degrees of freedom
            addpath('fminuit');
            persistent InvC2cdf; persistent scale; persistent dftmp;
            fprintf(2,'nmkChiFunc.myInvChi2Cdf - might give wrong values in certain scenarios.\n');
            if( isempty(dftmp) ); dftmp=df; end
            if( isempty(scale) || isempty(InvC2cdf) || (df~=dftmp) )
                dftmp=df;
                scale=logspace(-5,2,100);
                InvC2cdf=mkChiFunc.myChi2Cdf(scale,df); %caching calculated cdf
            end
            out = interp1(InvC2cdf, scale , prob,'PCHIP');
        end
        
        function out = myChi2Inv(x, df)
            %Returns the inverse chi2 pdf for given position x and degrees of freedom df.
            addpath('fminuit');
            out=2^(-df/2)/gammac(df/2).*x.^(-df/2-1).*exp(-1/2./x);
        end
        
        function out = myChi2InvCdf(x, df)
            %Returns the cumulative of the inverse chi2 pdf for given position x and degrees of
            %freedom df.
            if ~isa(x,'double'); error('2nd argument should be of type double.'); end
            if ~isa(df,'double'); error('3rd argument should be of type double.'); end
            
            fun=@(y) mkChiFunc.myChi2Inv(y,df);
            funfun=@(xmax) integral(fun,0,xmax);
            out = arrayfun(funfun,x);
        end
    end
    
    
    
    %% Functions for sensitivity studies. - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    methods(Static)
        
        function out = GetChiCL(prob,df)
            %out = chi2inv( prob, df ); %#Matlab-Toolbox
            
            %out = mkChiFunc.myInvChi2Cdf( prob, df ); %might be unsafe in certain scenarios
            %out = mkChiFunc.myInvChi2CdfRec( prob, df ); %rather slow
            
            %fct = memoize(@mkChiFunc.myInvChi2Cdf);    out=fct( prob, df );
            fct = memoize(@mkChiFunc.myInvChi2CdfRec); out=fct(prob,df);
            
        end
        
        function out = GetMassStatCL(chi2Arr,massArr,prob,df)
            % ... calculating the statistical effective mass - chi2 - equivalent
            %   prob    = 0.9      -> look for 90% CL chi2-meff-equivalent
            %   df      = 3        -> assuming 3 nuisance parameters
            %   massArr = [0,1,..] -> array of tested neutrino masses
            %   chi2Arr = [0,1,..] -> array of corresponding (chi2=chi2(mass)) chi2-values
            out=interp1( chi2Arr-min(chi2Arr), massArr, mkChiFunc.GetChiCL(prob,df),'linear' );
        end
        
        function out = ProbNorm(x)
            % ... calculating the symmetrical coverage probability around the average value of a
            % normal distribution.
            %out = normcdf(x,0,1)-normcdf(-x,0,1); %#Matlab-Toolbox
            
            %out = mkChiFunc.myGausCdf(x,0,1)-mkChiFunc.myGausCdf(-x,0,1);
            fct = memoize(@mkChiFunc.myGausCdf); out= fct(x,0,1)-fct(-x,0,1);
        end
        
        function out = CL2Sig(probCL)
            % ... converting CL values to multiples of sigma
            out = interp1(mkChiFunc.ProbNorm( linspace(0,5,1000) ) , linspace(0,5,1000) , probCL );
        end
        
        function out = GetMassTotCL(chi2Arr,massArr,prob,df)
            % ... calculating the total effective neutrino mass - chi2 - equivalent
            CL_m2stat=mkChiFunc.GetMassStatCL( chi2Arr,massArr,prob,df ).^2;
            CL2SigFac=mkChiFunc.CL2Sig(prob);
            
            sigma_m2syste = 0.017; %eV^2
            out = sqrt( ...
                sqrt( (CL_m2stat./CL2SigFac).^2 + sigma_m2syste^2 ) .*CL2SigFac ...
                ...%sqrt( (CL_m2stat).^2 + sigma_m2syste^2 .* CL2SigFac.^2 ) ...
                );
        end
        
        function     print(x)
            % ... printing a matlab arrays
            % x = [1,2,3,4,5], y = [6,7,8,9,10]
            % print(x)     -> |  1.0000  |  2.0000  |  3.0000  |  4.0000  |  5.0000  |
            % print([x;y]) -> |  1.0000  |  6.0000  |  2.0000  |  7.0000  |  3.0000  |
            %               > |  8.0000  |  4.0000  |  9.0000  |  5.0000  |  10.0000  |
            if(isfloat(x))
                fprintf(['|  ', repmat('%5.4f  |  ',1,length(x)), '\n'], x);
            end
        end
        
    end
    
end