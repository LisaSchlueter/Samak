%Debug class for relic neutrino analysis

classdef RelicNuDebug < handle
    
   properties
       Params;      %parameters of modelobj
       R;           %TBD obj
       ToggleES;    %whether to use excited states in the neutrino capture spectrum
       etaSensitivity;
   end
   
   methods %constructor
       function obj = RelicNuDebug(varargin)
           p = inputParser;
           p.addParameter('R','',@(x)isa(x,'TBD') || isempty(x));
           p.addParameter('Params','TDR',@(x)ismember(x,{'TDR','KNM1','KNM2','Formaggio'}));
           p.addParameter('ToggleES','OFF',@(x)ismember(x,{'ON','OFF'}));
           p.addParameter('etaSensitivity','',@(x)isfloat(x) || isempty(x));
           p.parse(varargin{:})
           obj.R              = p.Results.R;
           obj.Params         = p.Results.Params;
           obj.ToggleES       = p.Results.ToggleES;
           obj.etaSensitivity = p.Results.etaSensitivity;
           
           if isempty(obj.R)
                if strcmp(obj.Params,'TDR')
                    obj.R = ref_RelicNuBkg_DesignReport('ToggleES',obj.ToggleES);
                elseif strcmp(obj.Params,'KNM1')
                    obj.R = ref_RelicNuBkg_KNM1('ToggleES',obj.ToggleES);
                elseif strcmp(obj.Params,'Formaggio')
                    obj.R = ref_RelicNuBkg_Formaggio('ToggleES',obj.ToggleES);
                else
                    sprintf('Params input not known!');
                end
           end
            
           obj.R.ComputeTBDDS;
           obj.R.ComputeTBDIS;
       end
   end
   methods %Modify settings
       function SetRhoDmnu2TimeEta(obj,varargin)                 %set variable parameters
           p = inputParser;
           p.addParameter('RhoD',5e17,@(x)isfloat(x) && x>0);   %column density [molecules/cm^2]
           p.addParameter('mnu2',0,@(x)isfloat(x) && x>=0);      %neutrino mass [eV]
           p.addParameter('Time',1,@(x)isfloat(x) && x>0);      %measurement time [s]
           p.addParameter('Eta',1,@(x)isfloat(x) && x>0);       %relic neutrino overdensity
           p.addParameter('RhoDMode','CD',@(x)ismember(x,{'CD','Mass','EffMass'}));     %whether to input column density directly, a total tritium mass or an effective tritium mass
           p.parse(varargin{:});
           RhoD = p.Results.RhoD;
           mnu2  = p.Results.mnu2;
           Time = p.Results.Time;
           Eta  = p.Results.Eta;
           RhoDMode = p.Results.RhoDMode;
           
           switch RhoDMode
               case 'Mass'
                   RhoD = RhoD/(pi*obj.R.WGTS_FTR_cm^2*obj.R.M*1e3*obj.R.WGTS_epsT);
               case 'EffMass'
                   RhoD = RhoD./((pi*obj.R.WGTS_FTR_cm^2*obj.R.M*1e3*obj.R.WGTS_epsT)...
                       .*0.5*(1-cos(asin(sqrt(obj.R.WGTS_B_T./obj.R.MACE_Bmax_T)))) ...    %angle of acceptance
                        .*(obj.R.FPD_MeanEff*obj.R.FPD_Coverage)...                         %detector efficiency and coverage
                        .*(numel(obj.R.FPD_PixList)/148)...                                 %number of pixels
                        .*0.3983);                                                          %zero scattering prob.
           end
           
           if strcmp(obj.Params,'TDR')
               obj.R = ref_RelicNuBkg_DesignReport('mnuSq_i',mnu2,'eta_i',Eta,'TimeSec',Time,'ToggleES',obj.ToggleES,'WGTS_CD_MolPerCm2',RhoD);
           elseif strcmp(obj.Params,'KNM1')
               obj.R = ref_RelicNuBkg_TDR('mnuSq_i',mnu2,'eta_i',Eta,'TimeSec',Time,'ToggleES',obj.ToggleES,'WGTS_CD_MolPerCm2',RhoD);
           elseif strcmp(obj.Params,'Formaggio')
               obj.R = ref_RelicNuBkg_Formaggio('mnuSq_i',mnu2,'eta_i',Eta,'TimeSec',Time,'ToggleES',obj.ToggleES,'WGTS_CD_MolPerCm2',RhoD);
           end
           obj.R.ComputeNormFactorTBDDS;
           obj.R.ComputeTBDDS;
           obj.R.ComputeTBDIS;
       end
   end
   
   methods %Diff and Int Spec.
       function TotalCountsStudy(obj,varargin)
           p = inputParser;
           p.addParameter('Plot','ON',@(x)ismember(x,{'ON','OFF'}));
           p.parse(varargin{:});
           Plot = p.Results.Plot;
           
           NormFactor = 0.5*(1-cos(asin(sqrt(obj.R.WGTS_B_T./obj.R.MACE_Bmax_T)))) ...      %angle of acceptance
                        .*(obj.R.FPD_MeanEff*obj.R.FPD_Coverage)...                         %detector efficiency and coverage
                        .*numel(obj.R.FPD_PixList)/148;                                     %number of pixels
           sprintf(' m_nu: %.1i eV \n Measurement time: %.2s s (%.2i yr) \n Total rho_D: %0.2g molecules/cm^2 \n Total N_T: %0.2g tritium molecules \n Total N_mol: %.2s mol tritium molecules \n Total m_T: %.2s microgram \n Effective rho_D: %0.2g molecules/cm^2 \n Effective N_T: %0.2g tritium molecules \n Effective N_mol: %0.2g mol tritium molecules \n Effective m_T: %0.2g microgram \n Total capture rate: %0.2g 1/s (%0.2g 1/yr) \n Total number of capture events: %0.2g \n Effective number of capture events: %0.2g \n In ground state: %0.2g',...
               sqrt(obj.R.mnuSq),...                                                                        %neutrino mass [eV]
               obj.R.TimeSec,...                                                                            %measurement time [s]
               obj.R.TimeSec./(365.242*24*3600),...                                                         %measurement time [yr]
               obj.R.WGTS_CD_MolPerCm2.*obj.R.WGTS_epsT,...                                                 %column density
               obj.R.WGTS_CD_MolPerCm2.*obj.R.WGTS_epsT.*pi*obj.R.WGTS_FTR_cm^2,...                         %number of tritium molecules
               obj.R.WGTS_CD_MolPerCm2.*obj.R.WGTS_epsT.*pi*obj.R.WGTS_FTR_cm^2./obj.R.NA,...               %number of tritium molecules [mol]
               pi*obj.R.WGTS_FTR_cm^2*obj.R.WGTS_CD_MolPerCm2*obj.R.M*1e9*obj.R.WGTS_epsT,...               %tritium mass [µg]
               obj.R.WGTS_CD_MolPerCm2.*obj.R.WGTS_epsT.*NormFactor,...                                     %effective column density
               obj.R.WGTS_CD_MolPerCm2.*obj.R.WGTS_epsT.*pi*obj.R.WGTS_FTR_cm^2.*NormFactor,...             %effective number of tritium molecules
               obj.R.WGTS_CD_MolPerCm2.*obj.R.WGTS_epsT.*pi*obj.R.WGTS_FTR_cm^2./obj.R.NA.*NormFactor,...   %effective number of tritium molecules [mol]
               pi*obj.R.WGTS_FTR_cm^2*obj.R.WGTS_CD_MolPerCm2*obj.R.M*1e9*obj.R.WGTS_epsT.*NormFactor,...   %effective tritium mass [µg]
               obj.R.NormFactorTBDDS_R./NormFactor,...                                                      %total neutrino capture rate [1/s]
               obj.R.NormFactorTBDDS_R./NormFactor.*365.242.*24.*3600,...                                   %total neutrino capture rate [1/yr]
               obj.R.NormFactorTBDDS_R./NormFactor.*obj.R.TimeSec,...                                       %total number of capture events
               obj.R.NormFactorTBDDS_R.*obj.R.TimeSec,...                                                   %effective number of capture events
               obj.R.NormFactorTBDDS_R.*obj.R.TimeSec.*obj.R.TTNormGS)                                      %effective number of capture events in ground state
           
           if strcmp(Plot,'ON')
               hR = semilogy((obj.R.Te-obj.R.Q),obj.R.TBDDS.*obj.R.TimeSec,'LineWidth',2,'Color','Black','LineStyle','-');
               hold on;
               hB = semilogy((obj.R.Te-obj.R.Q),(obj.R.TBDDS-obj.R.TBDDS_R).*obj.R.TimeSec,'LineWidth',2,'Color','Blue','LineStyle',':');
               hP = semilogy((obj.R.Te-obj.R.Q),obj.R.TBDDS_R.*obj.R.TimeSec,'LineWidth',2,'Color','Red','LineStyle','--');
               grid on;
               ylim([1e-20,inf]);
               xlabel('E-E_0 (eV)','FontSize',12);
               str = sprintf('Count per energy bin');
               ylabel(str,'FontSize',14);
               strT1 = sprintf('Tritium Beta Decay Spectrum');
               strT2 = sprintf('Beta Decay + E-capture');
               str1 = sprintf('E-capture: %.3g evts',obj.R.NormFactorTBDDS_R.*obj.R.TimeSec);
               lh1 = legend([hR hB hP],strT2,strT1,str1);
               %legend(lh1,'box','off');
               set(lh1,'FontSize',12);
               %axis([-1 2 1e-2 10000])
               title(sprintf('Relic Neutrinos Capture rate - %.1f years',obj.R.TimeSec/(365.242*24*3600)),'FontSize',14);
               PrettyFigureFormat;
               hold off;
           end
       end
       function IntegralSpectrum(obj,varargin)
           p = inputParser;
           p.addParameter('PlotIsolatedRelicSpectrum','OFF',@(x)ismember(x,{'ON','OFF'}));
           p.parse(varargin{:});
           PlotIsolatedRelicSpectrum = p.Results.PlotIsolatedRelicSpectrum;
           
           TBDIS_R = obj.R.TBDIS./obj.R.qUfrac./obj.R.TimeSec;
           TBDIS_B = (obj.R.TBDIS-obj.R.TBDIS_R)./obj.R.qUfrac./obj.R.TimeSec;
           
           if strcmp(PlotIsolatedRelicSpectrum,'ON')
            TBDIS_R_R = obj.R.TBDIS_R./obj.R.qUfrac./obj.R.TimeSec;
           end
           
           e1 = obj.R.qU(obj.R.qU(:,1)>(obj.R.Q-310),:)-18575;
           tmpis1 = TBDIS_R(obj.R.qU(:,1)>(obj.R.Q-310),:);
           tmpis2 = TBDIS_B(obj.R.qU(:,1)>(obj.R.Q-310),:);
           
           if strcmp(PlotIsolatedRelicSpectrum,'ON')
            ris = TBDIS_R_R(obj.R.qU(:,1)>(obj.R.Q-310),:);
           end
           
           hT1 = semilogy((e1),tmpis1,'-s','LineWidth',2,'Color','Black','LineStyle','-');
           hold on;
           hT2 = semilogy((e1),tmpis2,'-s','LineWidth',2,'Color','Blue','LineStyle','-');
           
           if strcmp(PlotIsolatedRelicSpectrum,'ON')
            hR = semilogy((e1),ris,'-s','LineWidth',2,'Color','Red','LineStyle',':');
           end
           
           grid on;
           %ylim([1e-20,1e2]);
           xlabel('E-E_0 (eV)','FontSize',12);
           str = sprintf('Event Rate (Hz)');
           ylabel(str,'FontSize',14);
           strT1 = sprintf('Beta Decay');
           strT2 = sprintf('Beta Decay + E-capture');
           switch PlotIsolatedRelicSpectrum
               case 'ON'
                   strR  = sprintf('Capture Spectrum');
                   lh1 = legend([hT1 hT2 hR],strT1,strT2,strR);
               case 'OFF'
                   lh1 = legend([hT1 hT2],strT1,strT2);
           end
           set(lh1,'FontSize',12);
           PrettyFigureFormat;
           hold off;
       end
   end
   
   methods %Chi2 scans
        function Chi2Scan_Fake(obj,varargin)
            p=inputParser;
            p.addParameter('Recompute',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('range',30,@(x)isfloat(x));
            p.addParameter('RunNr',1,@(x)isfloat(x));                          % 1 for default setings from initfile, 10 if you want to modify settings
            p.addParameter('Init_Opt','',@(x)iscell(x) || isempty(x));         % cell array containing options to change in initfile (only relevant if RunNr==10)
            p.addParameter('fitPar','mNu E0 Norm Bck',@(x)ischar(x));
            p.addParameter('Syst','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('Netabins',10,@(x)isfloat(x));
            p.addParameter('etarange',10,@(x)isfloat(x));
            p.addParameter('etafactor',1.5,@(x)isfloat(x));                    % max(eta)=etafactor*10^etarange
            p.addParameter('mode','SCAN',@(x)ismember(x,{'SCAN','SEARCH'}));
            %% =========== SEARCH mode settings =============
            p.addParameter('etalower',0,@(x)isfloat(x));
            p.addParameter('etaupper',1.5e10,@(x)isfloat(x));                  % initial upper and lower search bounds
            p.addParameter('delta',0.1e9,@(x)isfloat(x));                      % amount by which to shift eta if fit fails
            p.parse(varargin{:});
            Recompute = p.Results.Recompute;
            RunNr     = p.Results.RunNr;
            range     = p.Results.range;
            Init_Opt  = p.Results.Init_Opt;
            fitPar    = p.Results.fitPar;
            Syst      = p.Results.Syst;
            Netabins  = p.Results.Netabins;
            etarange  = p.Results.etarange;
            etafactor = p.Results.etafactor;
            mode      = p.Results.mode;
            etalower  = p.Results.etalower;
            etaupper  = p.Results.etaupper;
            delta     = p.Results.delta;
            
            if strcmp(obj.Params,'TDR')
                initfile=@ref_RelicNuBkg_DesignReport;
                TwinBias_Q=18575;
                RingList=1:14;
            elseif strcmp(obj.Params,'Formaggio')
                initfile=@ref_RelicNuBkg_Formaggio;
                TwinBias_Q=18575;
                RingList=1:14;
            elseif strcmp(obj.Params,'KNM1')
                initfile=@ref_RelicNuBkg_KNM1;
                TwinBias_Q=18573.73;
                RingList=1:12;
            end
            
            if strcmp(Syst,'ON')
                Chi2opt='chi2CMShape';
            else
                Chi2opt='chi2Stat';
            end

            if RunNr==1
                U = RunAnalysis('RunNr',RunNr,...             
                    'FakeInitFile',initfile,...
                    'chi2',Chi2opt,...                 % uncertainties: statistical or stat + systematic uncertainties
                    'DataType','Fake',...                 % can be 'Real' or 'Twin' -> Monte Carlo
                    'RingList',RingList,...
                    'TwinBias_Q',TwinBias_Q,...
                    'fixPar',fitPar,...                   % free Parameter!!
                    'NonPoissonScaleFactor',1,...         % background uncertainty are enhanced
                    'minuitOpt','min ; minos',...         % technical fitting options (minuit)
                    'FSDFlag','SibilleFull',...          % final state distribution                        !!check ob initfile hier überschrieben wird
                    'ELossFlag','KatrinT2',...            % energy loss function
                    'SysBudget',22,...                    % defines syst. uncertainties -> in GetSysErr.m;
                    'DopplerEffectFlag','FSD',...
                    'SynchrotronFlag','OFF',...
                    'AngularTFFlag','OFF');
            elseif RunNr==10
                U = RunAnalysis('RunNr',RunNr,...             
                    'FakeInitFile',initfile,...
                    'Init_Opt',Init_Opt,...
                    'RecomputeFakeRun','ON',...
                    'chi2',Chi2opt,...                 % uncertainties: statistical or stat + systematic uncertainties
                    'DataType','Fake',...                 % can be 'Real' or 'Twin' -> Monte Carlo
                    'RingList',RingList,...
                    'TwinBias_Q',TwinBias_Q,...
                    'fixPar',fitPar,...                   % free Parameter!!
                    'NonPoissonScaleFactor',1,...         % background uncertainty are enhanced
                    'minuitOpt','min ; minos',...         % technical fitting options (minuit)
                    'FSDFlag','SibilleFull',...          % final state distribution                        !!check ob initfile hier überschrieben wird
                    'ELossFlag','KatrinT2',...            % energy loss function
                    'SysBudget',22,...                    % defines syst. uncertainties -> in GetSysErr.m;
                    'DopplerEffectFlag','FSD',...
                    'SynchrotronFlag','OFF',...
                    'AngularTFFlag','OFF');
            end

            U.exclDataStart = U.GetexclDataStart(range); % set region of interest
            
            if RunNr==10
                U.InitModelObj_Norm_BKG('Recompute','ON');
            end
            
            if strcmp(mode,'SCAN')
                Chi2      = 1:Netabins;
                mnuSq     = 1:Netabins;
                mnuSq_err = 1:Netabins;
                E0        = 1:Netabins;
                E0_err    = 1:Netabins;
                Bkg       = 1:Netabins;
                Bkg_err   = 1:Netabins;
                Norm      = 1:Netabins;
                Norm_err  = 1:Netabins;
                matFilePath = [getenv('SamakPath'),sprintf('RelicNuBkg/Chi2Scans/')];
                if RunNr==1
                    savename=[matFilePath,sprintf('RelicChi2Scan_Fake_Syst%s_range%g_%s_[0 %g]_%s.mat',Syst,range,obj.Params,etafactor*10^etarange,fitPar)];
                else
                    SaveStr='';
                    for i=1:numel(Init_Opt)
                        if ischar(Init_Opt{i})
                            SaveStr=[SaveStr,sprintf('_%s',Init_Opt{i})];
                        elseif isfloat(Init_Opt{i})
                            SaveStr=[SaveStr,sprintf('_%f',Init_Opt{i})];
                        end
                    end
                savename=[matFilePath,sprintf('RelicChi2Scan_Fake_Syst%s_range%g_%s_[0 %g]%s_%s.mat',Syst,range,obj.Params,etafactor*10^etarange,SaveStr,fitPar)];
                end
                
                if exist(savename,'file') && strcmp(Recompute,'OFF')
                    plotchi2scan(savename);
                else
                    for i=1:Netabins
                       U.ModelObj.eta_i = (i-1)*((etafactor*10^(etarange))/(Netabins-1));
                       U.ModelObj.eta   = (i-1)*((etafactor*10^(etarange))/(Netabins-1));
                       U.ModelObj.ComputeNormFactorTBDDS;
                       U.ModelObj.ComputeTBDDS;
                       U.ModelObj.ComputeTBDIS;
                       %% Fit
                       U.Fit;
                       Chi2(i)      = U.FitResult.chi2min;
                       mnuSq(i)     = U.ModelObj.mnuSq_i+U.FitResult.par(1);
                       mnuSq_err(i) = U.FitResult.err(1);
                       E0(i)        = U.ModelObj.Q_i+U.FitResult.par(2);
                       E0_err(i)    = U.FitResult.err(2);
                       Bkg(i)       = U.ModelObj.BKG_RateSec_i+U.FitResult.par(3);
                       Bkg_err(i)   = U.FitResult.err(3);
                       Norm(i)      = U.FitResult.par(3+U.ModelObj.nPixels:3+2*U.ModelObj.nPixels-1) + 1;
                       Norm_err(i)  = U.FitResult.err(3+U.ModelObj.nPixels:3+2*U.ModelObj.nPixels-1);
                    end

                    save(savename,'Chi2','Netabins','etafactor','etarange','mnuSq','mnuSq_err','E0','E0_err','Bkg','Bkg_err','Norm','Norm_err');
                    plotchi2scan(savename);
                end
            end

            if strcmp(mode,'SEARCH')
                
                matFilePath = [getenv('SamakPath'),sprintf('RelicNuBkg/UpperLimits/')];
                if RunNr==1
                    savename=[matFilePath,sprintf('RelicLimit_Fake_Syst%s_range%g_%s_%s.mat',Syst,range,obj.Params,fitPar)];
                else
                    SaveStr='';
                    for i=1:numel(Init_Opt)
                        if ischar(Init_Opt{i})
                            SaveStr=[SaveStr,sprintf('_%s',Init_Opt{i})];
                        elseif isfloat(Init_Opt{i})
                            SaveStr=[SaveStr,sprintf('_%f',Init_Opt{i})];
                        end
                    end
                    savename=[matFilePath,sprintf('RelicLimit_Fake_Syst%s_range%g_%s%s_%s.mat',Syst,range,obj.Params,SaveStr,fitPar)];
                end
                
                if exist(savename,'file') && strcmp(Recompute,'OFF')
                    load(savename,'eta');
                    switch RunNr
                        case 1
                            relic_global('eta',eta,'Params',obj.Params,'Syst',Syst);
                        case 10
                            relic_global('eta',eta,'Params',obj.Params,'Init_Opt',Init_Opt,'Syst',Syst);
                    end
                    sprintf('Final Result: eta = %g',eta)
                    obj.etaSensitivity = eta;
                else
                    if RunNr==1
                        F = RunAnalysis('RunNr',RunNr,...             
                            'FakeInitFile',initfile,...
                            'chi2',Chi2opt,...                 % uncertainties: statistical or stat + systematic uncertainties
                            'DataType','Fake',...                 % can be 'Real' or 'Twin' -> Monte Carlo
                            'RingList',RingList,...
                            'TwinBias_Q',TwinBias_Q,...
                            'fixPar',fitPar,...                   % free Parameter!!
                            'NonPoissonScaleFactor',1,...         % background uncertainty are enhanced
                            'minuitOpt','min ; minos',...         % technical fitting options (minuit)
                            'FSDFlag','SibilleFull',...          % final state distribution                        !!check ob initfile hier überschrieben wird
                            'ELossFlag','KatrinT2',...            % energy loss function
                            'SysBudget',22,...                    % defines syst. uncertainties -> in GetSysErr.m;
                            'DopplerEffectFlag','FSD',...
                            'SynchrotronFlag','OFF',...
                            'AngularTFFlag','OFF');
                    elseif RunNr==10
                        F = RunAnalysis('RunNr',RunNr,...             
                            'FakeInitFile',initfile,...
                            'Init_Opt',Init_Opt,...
                            'RecomputeFakeRun','ON',...
                            'chi2',Chi2opt,...                 % uncertainties: statistical or stat + systematic uncertainties
                            'DataType','Fake',...                 % can be 'Real' or 'Twin' -> Monte Carlo
                            'RingList',RingList,...
                            'TwinBias_Q',TwinBias_Q,...
                            'fixPar',fitPar,...                   % free Parameter!!
                            'NonPoissonScaleFactor',1,...         % background uncertainty are enhanced
                            'minuitOpt','min ; minos',...         % technical fitting options (minuit)
                            'FSDFlag','SibilleFull',...          % final state distribution                        !!check ob initfile hier überschrieben wird
                            'ELossFlag','KatrinT2',...            % energy loss function
                            'SysBudget',22,...                    % defines syst. uncertainties -> in GetSysErr.m;
                            'DopplerEffectFlag','FSD',...
                            'SynchrotronFlag','OFF',...
                            'AngularTFFlag','OFF');
                    end

                    F.exclDataStart = F.GetexclDataStart(range); % set region of interest
                    
                    if RunNr==10
                        U.InitModelObj_Norm_BKG('Recompute','ON');
                    end

                    eta=etaupper;
                    F.ModelObj.eta_i=etalower;
                    F.ModelObj.eta=etalower;
                    F.ModelObj.ComputeNormFactorTBDDS;
                    F.ModelObj.ComputeTBDDS;
                    F.ModelObj.ComputeTBDIS;
                    U.ModelObj.eta_i=etaupper;
                    U.ModelObj.eta=etaupper;
                    U.ModelObj.ComputeNormFactorTBDDS;
                    U.ModelObj.ComputeTBDDS;
                    U.ModelObj.ComputeTBDIS;
                    U.Fit;
                    F.Fit;
                    chi2lower = F.FitResult.chi2min;
                    chi2upper = U.FitResult.chi2min;
                    chi2=chi2upper;
                    if chi2upper<2.71
                        sprintf('chi^2 too small! Set larger initial eta.')
                    elseif chi2upper>50
                        sprintf('Minos failed! Vary initial eta.')
                    else
                        
                        while abs(chi2-2.71)>0.01
                            if chi2>2.71
                                etaupper=eta;
                                chi2upper=chi2;
                                eta=((etaupper-etalower)./(chi2upper-chi2lower)).*2.71-chi2upper.*((etaupper-etalower)./(chi2upper-chi2lower))+etaupper;
                                U.ModelObj.eta_i=eta;
                                U.ModelObj.eta=eta;
                                U.ModelObj.ComputeNormFactorTBDDS;
                                U.ModelObj.ComputeTBDDS;
                                U.ModelObj.ComputeTBDIS;
                                U.Fit;
                                if U.FitResult.chi2min>50
                                    sprintf('Fit failed! Retrying...')
                                    etaupper = etaupper + delta;
                                    eta=etaupper;
                                else
                                    chi2=U.FitResult.chi2min;
                                end
                            else
                                etalower=eta;
                                chi2lower=chi2;
                                eta=((etaupper-etalower)./(chi2upper-chi2lower)).*2.71-chi2upper.*((etaupper-etalower)./(chi2upper-chi2lower))+etaupper;
                                F.ModelObj.eta_i=eta;
                                F.ModelObj.eta=eta;
                                F.ModelObj.ComputeNormFactorTBDDS;
                                F.ModelObj.ComputeTBDDS;
                                F.ModelObj.ComputeTBDIS;
                                F.Fit;
                                if F.FitResult.chi2min>50
                                    sprintf('Fit failed! Retrying...')
                                    etalower = etalower - delta;
                                    eta = etalower;
                                else
                                    chi2=F.FitResult.chi2min;
                                end
                            end
                        end
                        save(savename,'eta');
                        switch RunNr
                            case 1
                                relic_global('eta',eta,'Params',obj.Params,'fitPar',fitPar,'E0',TwinBias_Q,'Syst',Syst);
                            case 10
                                relic_global('eta',eta,'Params',obj.Params,'fitPar',fitPar,'Init_Opt',Init_Opt,'E0',TwinBias_Q,'Syst',Syst);
                        end
                        sprintf('Final Result: eta = %g',eta)
                        obj.etaSensitivity = eta;
                    end
                end
            end
        end
       
        function Chi2Scan_Twin(obj,varargin)
            p=inputParser;
            p.addParameter('Recompute',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('RunList','KNM1',@(x)ischar(x));                          % KNM1 or KNM2
            p.addParameter('fitPar','mNu E0 Norm Bck',@(x)ischar(x));
            p.addParameter('Syst','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('SystBudget',24,@(x)isfloat(x));
            p.addParameter('TwinBias_mnuSq',0,@(x)isfloat(x));
            p.addParameter('range',30,@(x)isfloat(x));
            p.addParameter('Netabins',10,@(x)isfloat(x));
            p.addParameter('etarange',10,@(x)isfloat(x));
            p.addParameter('etafactor',1.5,@(x)isfloat(x));                    % max(eta)=etafactor*10^etarange
            p.addParameter('mode','SCAN',@(x)ismember(x,{'SCAN','SEARCH'}));
            %% =========== SEARCH mode settings =============
            p.addParameter('etalower',0,@(x)isfloat(x));
            p.addParameter('etaupper',1.5e10,@(x)isfloat(x));                  % initial upper and lower search bounds
            p.addParameter('delta',0.1e9,@(x)isfloat(x));                      % amount by which to shift eta if fit fails
            p.parse(varargin{:});
            Recompute      = p.Results.Recompute;
            RunList        = p.Results.RunList;
            fitPar         = p.Results.fitPar;
            Syst           = p.Results.Syst;
            SystBudget     = p.Results.SystBudget;
            TwinBias_mnuSq = p.Results.TwinBias_mnuSq;
            range          = p.Results.range;
            Netabins       = p.Results.Netabins;
            etarange       = p.Results.etarange;
            etafactor      = p.Results.etafactor;
            mode           = p.Results.mode;
            etalower       = p.Results.etalower;
            etaupper       = p.Results.etaupper;
            delta          = p.Results.delta;
            
            if strcmp(Syst,'ON')
                Chi2opt='chi2CMShape';
            else
                Chi2opt='chi2Stat';
            end
            

            if strcmp(RunList,'KNM1')
                U = MultiRunAnalysis('RunList',RunList,... % runlist defines which runs are analysed -> set MultiRunAnalysis.m -> function: GetRunList()
                    'chi2',Chi2opt,...                 % uncertainties: statistical or stat + systematic uncertainties
                    'DataType','Twin',...                 % can be 'Real' or 'Twin' -> Monte Carlo
                    'fixPar',fitPar,...                   % free Parameter!!
                    'RadiativeFlag','ON',...              % theoretical radiative corrections applied in model
                    'NonPoissonScaleFactor',1.064,...     % background uncertainty are enhanced
                    'minuitOpt','min ; minos',...         % technical fitting options (minuit)
                    'FSDFlag','SibilleFull',...          % final state distribution
                    'ELossFlag','KatrinT2',...            % energy loss function
                    'SysBudget',SystBudget,...                    % defines syst. uncertainties -> in GetSysErr.m;
                    'DopplerEffectFlag','FSD',...
                    'Twin_SameCDFlag','OFF',...
                    'Twin_SameIsotopFlag','OFF',...
                    'SynchrotronFlag','ON',...
                    'AngularTFFlag','OFF',...
                    'TwinBias_Q',18573.73,...
                    'TwinBias_mnuSq',TwinBias_mnuSq);
            elseif strcmp(RunList,'KNM2')
                U = MultiRunAnalysis('RunList',RunList,... % runlist defines which runs are analysed -> set MultiRunAnalysis.m -> function: GetRunList()
                    'chi2',Chi2opt,...                 % uncertainties: statistical or stat + systematic uncertainties
                    'DataType','Twin',...                 % can be 'Real' or 'Twin' -> Monte Carlo
                    'fixPar',fitPar,...                   % free Parameter!!
                    'RadiativeFlag','ON',...              % theoretical radiative corrections applied in model
                    'NonPoissonScaleFactor',1.064,...     % background uncertainty are enhanced
                    'minuitOpt','min ; minos',...         % technical fitting options (minuit)
                    'FSDFlag','SibilleFull',...          % final state distribution
                    'ELossFlag','KatrinT2',...            % energy loss function
                    'SysBudget',SystBudget,...                    % defines syst. uncertainties -> in GetSysErr.m;
                    'DopplerEffectFlag','FSD',...
                    'Twin_SameCDFlag','OFF',...
                    'Twin_SameIsotopFlag','OFF',...
                    'SynchrotronFlag','ON',...
                    'AngularTFFlag','OFF',...
                    'TwinBias_Q',18573.73,...
                    'TwinBias_mnuSq',TwinBias_mnuSq);
            end

            U.exclDataStart = U.GetexclDataStart(range); % set region of interest
            U.InitModelObj_Norm_BKG('Recompute','ON');
            
            if strcmp(mode,'SCAN')
                Chi2      = 1:Netabins;
                mnuSq     = 1:Netabins;
                mnuSq_err = 1:Netabins;
                E0        = 1:Netabins;
                E0_err    = 1:Netabins;
                Bkg       = 1:Netabins;
                Bkg_err   = 1:Netabins;
                Norm      = 1:Netabins;
                Norm_err  = 1:Netabins;
                
                matFilePath = [getenv('SamakPath'),sprintf('RelicNuBkg/Chi2Scans/')];
                savename=[matFilePath,sprintf('RelicChi2Scan_Twin_BiasmnuSq%g_Syst%s_range%g_%s_[0 %g]_%s.mat',TwinBias_mnuSq,Syst,range,obj.Params,etafactor*10^etarange,fitPar)];
                
                if exist(savename,'file') && strcmp(Recompute,'OFF')
                    plotchi2scan(savename);
                else
                    for i=1:Netabins
                       U.ModelObj.eta_i = (i-1)*((etafactor*10^(etarange))/(Netabins-1));
                       U.ModelObj.eta   = (i-1)*((etafactor*10^(etarange))/(Netabins-1));
                       U.ModelObj.ComputeNormFactorTBDDS;
                       U.ModelObj.ComputeTBDDS;
                       U.ModelObj.ComputeTBDIS;
                       %% Fit
                       U.Fit;
                       Chi2(i)      = U.FitResult.chi2min;
                       mnuSq(i)     = U.ModelObj.mnuSq_i+U.FitResult.par(1);
                       mnuSq_err(i) = U.FitResult.err(1);
                       E0(i)        = U.ModelObj.Q_i+U.FitResult.par(2);
                       E0_err(i)    = U.FitResult.err(2);
                       Bkg(i)       = U.ModelObj.BKG_RateSec_i+U.FitResult.par(3);
                       Bkg_err(i)   = U.FitResult.err(3);
                       Norm(i)      = U.FitResult.par(3+U.ModelObj.nPixels:3+2*U.ModelObj.nPixels-1) + 1;
                       Norm_err(i)  = U.FitResult.err(3+U.ModelObj.nPixels:3+2*U.ModelObj.nPixels-1);
                    end

                    save(savename,'Chi2','Netabins','etafactor','etarange','mnuSq','mnuSq_err','E0','E0_err','Bkg','Bkg_err','Norm','Norm_err');
                    plotchi2scan(savename);
                end
            end

            if strcmp(mode,'SEARCH')
                
                matFilePath = [getenv('SamakPath'),sprintf('RelicNuBkg/UpperLimits/')];
                savename=[matFilePath,sprintf('RelicLimit_Twin_BiasmnuSq%g_Syst%s_range%g_%s_%s.mat',TwinBias_mnuSq,Syst,range,obj.Params,fitPar)];
                
                if exist(savename,'file') && strcmp(Recompute,'OFF')
                    load(savename,'eta');
                    relic_global_twin('eta',eta,'RunList',obj.Params,'fitPar',fitPar,'E0',U.TwinBias_Q,'mnuSq',U.TwinBias_mnuSq,'Syst',Syst);
                    sprintf('Final Result: eta = %g',eta)
                    obj.etaSensitivity = eta;
                else
                    if strcmp(RunList,'KNM1')
                        F = MultiRunAnalysis('RunList',RunList,... % runlist defines which runs are analysed -> set MultiRunAnalysis.m -> function: GetRunList()
                             'chi2',Chi2opt,...                 % uncertainties: statistical or stat + systematic uncertainties
                             'DataType','Twin',...                 % can be 'Real' or 'Twin' -> Monte Carlo
                             'fixPar',fitPar,...                   % free Parameter!!
                             'RadiativeFlag','ON',...              % theoretical radiative corrections applied in model
                             'NonPoissonScaleFactor',1.064,...     % background uncertainty are enhanced
                             'minuitOpt','min ; minos',...         % technical fitting options (minuit)
                             'FSDFlag','SibilleFull',...          % final state distribution
                             'ELossFlag','KatrinT2',...            % energy loss function
                             'SysBudget',SystBudget,...                    % defines syst. uncertainties -> in GetSysErr.m;
                             'DopplerEffectFlag','FSD',...
                             'Twin_SameCDFlag','OFF',...
                             'Twin_SameIsotopFlag','OFF',...
                             'SynchrotronFlag','ON',...
                             'AngularTFFlag','OFF',...
                             'TwinBias_Q',18573.73,...
                             'TwinBias_mnuSq',TwinBias_mnuSq);
                    elseif strcmp(RunList,'KNM2')
                        F = MultiRunAnalysis('RunList',RunList,... % runlist defines which runs are analysed -> set MultiRunAnalysis.m -> function: GetRunList()
                             'chi2',Chi2opt,...                 % uncertainties: statistical or stat + systematic uncertainties
                             'DataType','Twin',...                 % can be 'Real' or 'Twin' -> Monte Carlo
                             'fixPar',fitPar,...                   % free Parameter!!
                             'RadiativeFlag','ON',...              % theoretical radiative corrections applied in model
                             'NonPoissonScaleFactor',1.064,...     % background uncertainty are enhanced
                             'minuitOpt','min ; minos',...         % technical fitting options (minuit)
                             'FSDFlag','SibilleFull',...          % final state distribution
                             'ELossFlag','KatrinT2',...            % energy loss function
                             'SysBudget',SystBudget,...                    % defines syst. uncertainties -> in GetSysErr.m;
                             'DopplerEffectFlag','FSD',...
                             'Twin_SameCDFlag','OFF',...
                             'Twin_SameIsotopFlag','OFF',...
                             'SynchrotronFlag','ON',...
                             'AngularTFFlag','OFF',...
                             'TwinBias_Q',18573.73,...
                             'TwinBias_mnuSq',TwinBias_mnuSq);
                    end

                    F.exclDataStart = F.GetexclDataStart(range); % set region of interest
                    F.InitModelObj_Norm_BKG('Recompute','ON');

                    eta=etaupper;
                    F.ModelObj.eta_i=etalower;
                    F.ModelObj.eta=etalower;
                    F.ModelObj.ComputeNormFactorTBDDS;
                    F.ModelObj.ComputeTBDDS;
                    F.ModelObj.ComputeTBDIS;
                    U.ModelObj.eta_i=etaupper;
                    U.ModelObj.eta=etaupper;
                    U.ModelObj.ComputeNormFactorTBDDS;
                    U.ModelObj.ComputeTBDDS;
                    U.ModelObj.ComputeTBDIS;
                    U.Fit;
                    F.Fit;
                    chi2upper = U.FitResult.chi2min;
                    chi2lower = F.FitResult.chi2min;
                    chi2 = chi2upper;
                    if chi2<2.71
                        sprintf('chi^2 too small! Set larger initial eta.')
                    elseif chi2>50
                        sprintf('Minos failed! Vary initial eta.')
                    else
                        while abs(chi2-2.71)>0.01
                            if chi2>2.71
                                etaupper=eta;
                                chi2upper=chi2;
                                eta=((etaupper-etalower)./(chi2upper-chi2lower)).*2.71-chi2upper.*((etaupper-etalower)./(chi2upper-chi2lower))+etaupper;
                                U.ModelObj.eta_i=eta;
                                U.ModelObj.eta=eta;
                                U.ModelObj.ComputeNormFactorTBDDS;
                                U.ModelObj.ComputeTBDDS;
                                U.ModelObj.ComputeTBDIS;
                                U.Fit;
                                if U.FitResult.chi2min>50
                                    sprintf('Fit failed! Retrying...')
                                    etaupper = etaupper + delta;
                                    eta=etaupper;
                                else
                                    chi2=U.FitResult.chi2min;
                                end
                            else
                                etalower=eta;
                                chi2lower=chi2;
                                eta=((etaupper-etalower)./(chi2upper-chi2lower)).*2.71-chi2upper.*((etaupper-etalower)./(chi2upper-chi2lower))+etaupper;
                                F.ModelObj.eta_i=eta;
                                F.ModelObj.eta=eta;
                                F.ModelObj.ComputeNormFactorTBDDS;
                                F.ModelObj.ComputeTBDDS;
                                F.ModelObj.ComputeTBDIS;
                                F.Fit;
                                if F.FitResult.chi2min>50
                                    sprintf('Fit failed! Retrying...')
                                    etalower = etalower - delta;
                                    eta = etalower;
                                else
                                    chi2=F.FitResult.chi2min;
                                end
                            end
                        end
                        save(savename,'eta');
                        relic_global_twin('eta',eta,'RunList',obj.Params,'fitPar',fitPar,'E0',U.TwinBias_Q,'mnuSq',U.TwinBias_mnuSq,'Syst',Syst);
                        sprintf('Final Result: eta = %g',eta)
                        obj.etaSensitivity = eta;
                    end
                end
            end
        end
   end
   
   methods
       function Chi2Twin(obj,varargin)
            p=inputParser;
            p.addParameter('range',40,@(x)isfloat(x));
            p.addParameter('fitPar','mNu E0 Norm Bkg',@(x)ischar(x));
            p.addParameter('Syst','ON',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('SystBudget',24,@(x)isfloat(x));
            p.addParameter('TwinBias_mnuSq',1,@(x)isfloat(x));
            p.addParameter('NetaBins',10,@(x)isfloat(x));
            p.addParameter('etarange',11,@(x)isfloat(x));
            p.addParameter('etafactor',5,@(x)isfloat(x));
            p.addParameter('Recompute','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.parse(varargin{:});

            range          = p.Results.range;
            fitPar         = p.Results.fitPar;
            Syst           = p.Results.Syst;
            SystBudget     = p.Results.SystBudget;
            TwinBias_mnuSq = p.Results.TwinBias_mnuSq;
            NetaBins       = p.Results.NetaBins;
            etarange       = p.Results.etarange;
            etafactor      = p.Results.etafactor;
            Recompute      = p.Results.Recompute;

            
            obj.Chi2Scan_Twin('Recompute',Recompute,...
                'range',range,...
                'RunList',obj.Params,...
                'fitPar',fitPar,...
                'Syst',Syst,...
                'SystBudget',SystBudget,...
                'TwinBias_mnuSq',TwinBias_mnuSq,...
                'Netabins',NetaBins,...
                'etarange',etarange,...
                'etafactor',etafactor,...
                'mode','SCAN');

            matFilePath = [getenv('SamakPath'),sprintf('RelicNuBkg/Chi2Scans/')];
            savename=[matFilePath,sprintf('RelicChi2Scan_Twin_BiasmnuSq%g_Syst%s_range%g_%s_[0 %g]_%s.mat',TwinBias_mnuSq,Syst,range,obj.Params,etafactor*10^etarange,fitPar)];

            load(savename);

            if Chi2(end)<2.71
                sprintf('Increase etafactor or etarange')
            else
                etavalues=(0:(Netabins-1))*((etafactor*10^(etarange))/(Netabins-1));
                m=find(Chi2<2.71);
                n=find(Chi2>2.71);
                etalower=etavalues(m(end));
                etaupper=etavalues(n(1));

                obj.Chi2Scan_Twin('Recompute',Recompute,...
                    'range',range,...
                    'RunList',obj.Params,...
                    'fitPar',fitPar,...
                    'Syst',Syst,...
                    'SystBudget',SystBudget,...
                    'TwinBias_mnuSq',TwinBias_mnuSq,...
                    'etalower',etalower,...
                    'etaupper',etaupper,...
                    'mode','SEARCH');
            end
        end
   end
    
end