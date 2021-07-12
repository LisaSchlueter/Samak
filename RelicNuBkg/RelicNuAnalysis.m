%Debug class for relic neutrino analysis

classdef RelicNuAnalysis < handle
    
   properties
       Params;      %parameters of modelobj
       R;           %TBD obj
       M;           %(Multi)RunAnalysis object
       ToggleES;    %whether to use excited states in the neutrino capture spectrum
       etaSensitivity;
   end
   
   methods %constructor
       function obj = RelicNuAnalysis(varargin)
           p = inputParser;
           p.addParameter('R','',@(x)isa(x,'TBD') || isempty(x));
           p.addParameter('M','',@(x)isa(x,'MultiRunAnalysis') || isa(x,'RunAnalysis') || isempty(x));
           p.addParameter('Params','TDR',@(x)ismember(x,{'TDR','KNM1','KNM2_Prompt','Formaggio'}));
           p.addParameter('ToggleES','OFF',@(x)ismember(x,{'ON','OFF'}));
           p.addParameter('etaSensitivity','',@(x)isfloat(x) || isempty(x));
           p.parse(varargin{:})
           obj.R              = p.Results.R;
           obj.M              = p.Results.M;
           obj.Params         = p.Results.Params;
           obj.ToggleES       = p.Results.ToggleES;
           obj.etaSensitivity = p.Results.etaSensitivity;
           
           if isempty(obj.R)
                if strcmp(obj.Params,'TDR')
                    obj.R = ref_RelicNuBkg_DesignReport('ToggleES',obj.ToggleES);
                elseif strcmp(obj.Params,'KNM1')
                    obj.R = ref_RelicNuBkg_KNM1('ToggleES',obj.ToggleES);
                elseif strcmp(obj.Params,'KNM2_Prompt')
                    obj.R = ref_RelicNuBkg_KNM2('ToggleES',obj.ToggleES);
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
       function SetParam(obj,varargin)
          p=inputParser;
          p.addParameter('Parameter','',@(x)ismember(x,{'mNu','E0','Norm','Bkg','eta'}));
          p.addParameter('value',0,@(x)isfloat(x));
          p.addParameter('INIT',0,@(x)isfloat(x));
          p.addParameter('SystSelect','OFF',@(x)ismember(x,{'RF','TASR','Stack','FSD','TC','FPDeff','LongPlasma','BkgCM','BkgPtCM','NP','None','OFF'}));
          p.addParameter('fitPar','mNu E0 Norm Bkg',@(x)ischar(x));
          p.parse(varargin{:});
          SystSelect = p.Results.SystSelect;
          Parameter  = p.Results.Parameter;
          value      = p.Results.value;
          INIT       = p.Results.INIT;
          fitPar     = p.Results.fitPar;
          
          fitPar = erase(fitPar,Parameter);
          
          if INIT==1
              if ~(strcmp(SystSelect,'NP') || strcmp(SystSelect,'OFF'))
                  obj.M.NonPoissonScaleFactor = 1;
              end
              if strcmp(SystSelect,'NP')
                  if strcmp(obj.Params,'KNM1')
                      obj.M.NonPoissonScaleFactor = 1.064;
                  elseif strcmp(obj.Params,'KNM2_Prompt')
                      obj.M.NonPoissonScaleFactor = 1.112;
                  elseif strcmp(obj.Params,'TDR')
                      obj.M.NonPoissonScaleFactor = 1;
                  end
              end
              if ~contains(fitPar,'mNu') && ~strcmp(Parameter,'mNu')
                  mNufix=obj.M.ModelObj.mnuSq_i;
              end
              range = obj.M.exclDataStart;
              obj.M = MultiRunAnalysis('RunList',obj.Params,... % runlist defines which runs are analysed -> set MultiRunAnalysis.m -> function: GetRunList()
                    'chi2',obj.M.chi2,...                 % uncertainties: statistical or stat + systematic uncertainties
                    'DataType',obj.M.DataType,...                 % can be 'Real' or 'Twin' -> Monte Carlo
                    'fixPar',fitPar,...                   % free Parameter!!
                    'RadiativeFlag','ON',...              % theoretical radiative corrections applied in model
                    'NonPoissonScaleFactor',obj.M.NonPoissonScaleFactor,...     % background uncertainty are enhanced
                    'fitter',obj.M.fitter,...
                    'minuitOpt','min ; minos',...         % technical fitting options (minuit)
                    'pullFlag',obj.M.pullFlag,...
                    'FSDFlag',obj.M.FSDFlag,...          % final state distribution
                    'ELossFlag',obj.M.ELossFlag,...            % energy loss function
                    'SysBudget',obj.M.SysBudget,...                    % defines syst. uncertainties -> in GetSysErr.m;
                    'AnaFlag','StackPixel',...
                    'DopplerEffectFlag','FSD',...
                    'FSD_Sigma',obj.M.FSD_Sigma,...
                    'Twin_SameCDFlag','OFF',...
                    'Twin_SameIsotopFlag','OFF',...
                    'SynchrotronFlag','ON',...
                    'AngularTFFlag',obj.M.AngularTFFlag,...
                    'TwinBias_Q',obj.M.TwinBias_Q,...
                    'TwinBias_mnuSq',obj.M.TwinBias_mnuSq,...
                    'BKG_PtSlope',obj.M.BKG_PtSlope,...
                    'TwinBias_BKG_PtSlope',obj.M.TwinBias_BKG_PtSlope,...
                    'TwinBias_FSDSigma',obj.M.TwinBias_FSDSigma);
                obj.M.exclDataStart = range;
          end
          
          if strcmp(Parameter,'eta') && INIT==1 && strcmp(obj.M.chi2,'chi2Stat')
              obj.M.ModelObj.mnuSq_i = obj.M.TwinBias_mnuSq;
              obj.M.InitModelObj_Norm_BKG('Recompute','ON');
          end
          
          if ~contains(fitPar,'mNu') && ~strcmp(Parameter,'mNu')
              obj.M.ModelObj.mnuSq_i=mNufix;
          end
          if strcmp(Parameter,'mNu')
              obj.M.ModelObj.mnuSq_i = value;
          elseif strcmp(Parameter,'E0')
              obj.M.ModelObj.Q_i = value;
          elseif strcmp(Parameter,'Bkg')
              obj.M.ModelObj.BKG_RateSec_i = value;
          elseif strcmp(Parameter,'eta')
              obj.M.ModelObj.eta_i = value;
          end    
          
          if ~strcmp(SystSelect,'OFF') && strcmp(obj.M.chi2,'chi2CMShape')
              if ~(strcmp(SystSelect,'BkgCM') || strcmp(SystSelect,'None'))
                  obj.M.ComputeCM('SysEffects',struct(SystSelect,'ON'),'BkgCM','OFF','BkgPtCM','OFF');
              elseif strcmp(SystSelect,'BkgCM')
                  obj.M.ComputeCM('SysEffects',struct(),'BkgCM','ON','BkgPtCM','OFF');
              elseif strcmp(SystSelect,'None')
                  obj.M.NonPoissonScaleFactor = 1;
                  obj.M.ComputeCM('SysEffects',struct(),'BkgCM','OFF','BkgPtCM','OFF');
              elseif strcmp(SystSelect,'NP')
                  obj.M.ComputeCM('SysEffects',struct(),'BkgCM','OFF','BkgPtCM','OFF');
              elseif strcmp(SystSelect,'BkgPtCM')
                  obj.M.ComputeCM('SysEffects',struct(),'BkgCM','OFF','BkgPtCM','ON');
              end
          end
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
                        .*numel(obj.R.FPD_PixList)/148 ...                                  %number of pixels
                        .*0.3983;                                                           %zero scattering prob.
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
               obj.R.NormFactorTBDDS_R./NormFactor.*0.3983,...                                              %total neutrino capture rate [1/s]
               obj.R.NormFactorTBDDS_R./NormFactor.*365.242.*24.*3600.*0.3983,...                           %total neutrino capture rate [1/yr]
               obj.R.NormFactorTBDDS_R./NormFactor.*obj.R.TimeSec.*0.3983,...                               %total number of capture events
               obj.R.NormFactorTBDDS_R.*obj.R.TimeSec.*0.3983,...                                           %effective number of capture events
               obj.R.NormFactorTBDDS_R.*obj.R.TimeSec.*obj.R.TTNormGS.*0.3983)                              %effective number of capture events in ground state
           
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
            p.addParameter('fitPar','mNu E0 Norm Bkg',@(x)ischar(x));
            p.addParameter('Syst','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('pullFlag',3);
            p.addParameter('Netabins',10,@(x)isfloat(x));
            p.addParameter('etarange',10,@(x)isfloat(x));
            p.addParameter('etafactor',1.5,@(x)isfloat(x));                    % max(eta)=etafactor*10^etarange
            p.addParameter('mode','SCAN',@(x)ismember(x,{'SCAN','SEARCH'}));
            p.addParameter('Plot','ON',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('mnuSqfix',@(x)isfloat(x) && x>0);
            %% =========== SEARCH mode settings =============
            p.addParameter('etalower',0,@(x)isfloat(x));
            p.addParameter('etaupper',1.5e10,@(x)isfloat(x));                  % initial upper and lower search bounds
            p.addParameter('delta',0.1e9,@(x)isfloat(x));                      % amount by which to shift eta if fit fails
            p.addParameter('DeltaChi2',2.71,@(x)isfloat(x));
            p.parse(varargin{:});
            Recompute = p.Results.Recompute;
            RunNr     = p.Results.RunNr;
            range     = p.Results.range;
            Init_Opt  = p.Results.Init_Opt;
            fitPar    = p.Results.fitPar;
            Syst      = p.Results.Syst;
            pullFlag  = p.Results.pullFlag;
            Netabins  = p.Results.Netabins;
            etarange  = p.Results.etarange;
            etafactor = p.Results.etafactor;
            mode      = p.Results.mode;
            Plot      = p.Results.Plot;
            etalower  = p.Results.etalower;
            etaupper  = p.Results.etaupper;
            delta     = p.Results.delta;
            DeltaChi2 = p.Results.DeltaChi2;
            mnuSqfix  = p.Results.mnuSqfix;
            
            if strcmp(obj.Params,'TDR')
                initfile=@ref_RelicNuBkg_DesignReport;
                TwinBias_Q=18575;
                RingList=1:14;
                SysBudget=67;
                NonPoissonScaleFactor=1;
            elseif strcmp(obj.Params,'Formaggio')
                initfile=@ref_RelicNuBkg_Formaggio;
                TwinBias_Q=18575;
                RingList=1:14;
            elseif strcmp(obj.Params,'KNM1')
                initfile=@ref_RelicNuBkg_KNM1;
                TwinBias_Q=18573.73;
                RingList=1:12;
                SysBudget=24;
                NonPoissonScaleFactor=1.064;
            elseif strcmp(obj.Params,'KNM2_Prompt')
                initfile=@ref_RelicNuBkg_KNM2;
                TwinBias_Q=18573.7;
                RingList=1:12;
                SysBudget=40;
                NonPoissonScaleFactor=1.1120;
            end
            
            if strcmp(Syst,'ON')
                Chi2opt='chi2CMShape';
            else
                Chi2opt='chi2Stat';
                NonPoissonScaleFactor=1;
            end

            if RunNr==1
                U = RunAnalysis('RunNr',RunNr,...             
                    'FakeInitFile',initfile,...
                    'chi2',Chi2opt,...                 % uncertainties: statistical or stat + systematic uncertainties
                    'NonPoissonScaleFactor',NonPoissonScaleFactor,...     % background uncertainty are enhanced
                    'DataType','Fake',...                 % can be 'Real' or 'Twin' -> Monte Carlo
                    'RingList',RingList,...
                    'TwinBias_Q',TwinBias_Q,...
                    'fixPar',fitPar,...                   % free Parameter!!
                    'minuitOpt','min ; minos',...         % technical fitting options (minuit)
                    'pullFlag',pullFlag,...
                    'FSDFlag','SibilleFull',...          % final state distribution                        !!check ob initfile hier überschrieben wird
                    'ELossFlag','KatrinT2',...            % energy loss function
                    'SysBudget',SysBudget,...                    % defines syst. uncertainties -> in GetSysErr.m;
                    'DopplerEffectFlag','FSD',...
                    'SynchrotronFlag','OFF',...
                    'AngularTFFlag','OFF');
            elseif RunNr==10
                U = RunAnalysis('RunNr',RunNr,...             
                    'FakeInitFile',initfile,...
                    'Init_Opt',Init_Opt,...
                    'RecomputeFakeRun','ON',...
                    'chi2',Chi2opt,...                 % uncertainties: statistical or stat + systematic uncertainties
                    'NonPoissonScaleFactor',NonPoissonScaleFactor,...     % background uncertainty are enhanced
                    'DataType','Fake',...                 % can be 'Real' or 'Twin' -> Monte Carlo
                    'RingList',RingList,...
                    'TwinBias_Q',TwinBias_Q,...
                    'fixPar',fitPar,...                   % free Parameter!!
                    'minuitOpt','min ; minos',...         % technical fitting options (minuit)
                    'pullFlag',pullFlag,...
                    'FSDFlag','SibilleFull',...          % final state distribution                        !!check ob initfile hier überschrieben wird
                    'ELossFlag','KatrinT2',...            % energy loss function
                    'SysBudget',SysBudget,...                    % defines syst. uncertainties -> in GetSysErr.m;
                    'DopplerEffectFlag','FSD',...
                    'SynchrotronFlag','OFF',...
                    'AngularTFFlag','OFF');
            end

            U.exclDataStart = U.GetexclDataStart(range); % set region of interest
            obj.M = U;
            
            if ~contains(fitPar,'mNu')
                U.ModelObj.mnuSq_i=mnuSqfix;
            end
            
            if RunNr==10 && strcmp(U.chi2,'chi2Stat')
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
                    savename=[matFilePath,sprintf('RelicChi2Scan_Fake_Syst%s_range%g_%s_[0 %g]_%s',Syst,range,obj.Params,etafactor*10^etarange,fitPar)];
                else
                    SaveStr='';
                    for i=1:numel(Init_Opt)
                        if ischar(Init_Opt{i})
                            SaveStr=[SaveStr,sprintf('_%s',Init_Opt{i})];
                        elseif isfloat(Init_Opt{i})
                            SaveStr=[SaveStr,sprintf('_%f',Init_Opt{i})];
                        end
                    end
                savename=[matFilePath,sprintf('RelicChi2Scan_Fake_Syst%s_range%g_%s_[0 %g]%s_%s',Syst,range,obj.Params,etafactor*10^etarange,SaveStr,fitPar)];
                end
                if any(ismember(pullFlag,1))
                    savename = [savename,'_mnuSqPull'];
                end
                savename = [savename,'.mat'];
                
                if exist(savename,'file') && strcmp(Recompute,'OFF') && strcmp(Plot,'ON')
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
                    if strcmp(Plot,'ON')
                        plotchi2scan(savename);
                    end
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
                    savename=[matFilePath,sprintf('RelicLimit_Fake_Syst%s_range%g_%s%s_%s',Syst,range,obj.Params,SaveStr,fitPar)];
                    if any(ismember(pullFlag,1))
                        savename = [savename,'_mnuSqPull'];
                    end
                    savename = [savename,'.mat'];
                end
                
                if exist(savename,'file') && strcmp(Recompute,'OFF')
                    load(savename,'eta');
                    if strcmp(Plot,'ON')
                        switch RunNr
                            case 1
                                relic_global('eta',eta,'Params',obj.Params,'Syst',Syst);
                            case 10
                                relic_global('eta',eta,'Params',obj.Params,'Init_Opt',Init_Opt,'Syst',Syst);
                        end
                    end
                    sprintf('Final Result: eta = %g',eta)
                    obj.etaSensitivity = eta;
                else
                    if RunNr==1
                        F = RunAnalysis('RunNr',RunNr,...             
                            'FakeInitFile',initfile,...
                            'chi2',Chi2opt,...                 % uncertainties: statistical or stat + systematic uncertainties
                            'NonPoissonScaleFactor',NonPoissonScaleFactor,...     % background uncertainty are enhanced
                            'DataType','Fake',...                 % can be 'Real' or 'Twin' -> Monte Carlo
                            'RingList',RingList,...
                            'TwinBias_Q',TwinBias_Q,...
                            'fixPar',fitPar,...                   % free Parameter!!
                            'minuitOpt','min ; minos',...         % technical fitting options (minuit)
                            'pullFlag',pullFlag,...
                            'FSDFlag','SibilleFull',...          % final state distribution                        !!check ob initfile hier überschrieben wird
                            'ELossFlag','KatrinT2',...            % energy loss function
                            'SysBudget',SysBudget,...                    % defines syst. uncertainties -> in GetSysErr.m;
                            'DopplerEffectFlag','FSD',...
                            'SynchrotronFlag','OFF',...
                            'AngularTFFlag','OFF');
                    elseif RunNr==10
                        F = RunAnalysis('RunNr',RunNr,...             
                            'FakeInitFile',initfile,...
                            'Init_Opt',Init_Opt,...
                            'RecomputeFakeRun','ON',...
                            'chi2',Chi2opt,...                 % uncertainties: statistical or stat + systematic uncertainties
                            'NonPoissonScaleFactor',NonPoissonScaleFactor,...     % background uncertainty are enhanced
                            'DataType','Fake',...                 % can be 'Real' or 'Twin' -> Monte Carlo
                            'RingList',RingList,...
                            'TwinBias_Q',TwinBias_Q,...
                            'fixPar',fitPar,...                   % free Parameter!!
                            'minuitOpt','min ; minos',...         % technical fitting options (minuit)
                            'pullFlag',pullFlag,...
                            'FSDFlag','SibilleFull',...          % final state distribution                        !!check ob initfile hier überschrieben wird
                            'ELossFlag','KatrinT2',...            % energy loss function
                            'SysBudget',SysBudget,...                    % defines syst. uncertainties -> in GetSysErr.m;
                            'DopplerEffectFlag','FSD',...
                            'SynchrotronFlag','OFF',...
                            'AngularTFFlag','OFF');
                    end

                    F.exclDataStart = F.GetexclDataStart(range); % set region of interest
                    
                    if ~contains(fitPar,'mNu')
                        F.ModelObj.mnuSq_i=mnuSqfix;
                    end
                    
                    if RunNr==10 && strcmp(U.chi2,'chi2Stat')
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
                    mnuSq     = U.ModelObj.mnuSq_i+U.FitResult.par(1);
                    mnuSq_err = U.FitResult.err(1);
                    E0        = U.ModelObj.Q_i+U.FitResult.par(2);
                    E0_err    = U.FitResult.err(2);
                    Bkg       = U.ModelObj.BKG_RateSec_i+U.FitResult.par(3);
                    Bkg_err   = U.FitResult.err(3);
                    Norm      = U.FitResult.par(3+U.ModelObj.nPixels:3+2*U.ModelObj.nPixels-1) + 1;
                    Norm_err  = U.FitResult.err(3+U.ModelObj.nPixels:3+2*U.ModelObj.nPixels-1);
                    if chi2upper<DeltaChi2
                        sprintf('chi^2 too small! Set larger initial eta.')
%                     elseif chi2upper>50
%                         sprintf('Minos failed! Vary initial eta.')
                    else
                        
                        while abs(chi2-DeltaChi2)>0.02
                            if chi2>DeltaChi2
                                etaupper=eta;
                                chi2upper=chi2;
                                eta=((etaupper-etalower)./(chi2upper-chi2lower)).*DeltaChi2-chi2upper.*((etaupper-etalower)./(chi2upper-chi2lower))+etaupper;
                                U.ModelObj.eta_i=eta;
                                U.ModelObj.eta=eta;
                                U.ModelObj.ComputeNormFactorTBDDS;
                                U.ModelObj.ComputeTBDDS;
                                U.ModelObj.ComputeTBDIS;
                                U.Fit;
                                mnuSq     = U.ModelObj.mnuSq_i+U.FitResult.par(1);
                                mnuSq_err = U.FitResult.err(1);
                                E0        = U.ModelObj.Q_i+U.FitResult.par(2);
                                E0_err    = U.FitResult.err(2);
                                Bkg       = U.ModelObj.BKG_RateSec_i+U.FitResult.par(3);
                                Bkg_err   = U.FitResult.err(3);
                                Norm      = U.FitResult.par(3+U.ModelObj.nPixels:3+2*U.ModelObj.nPixels-1) + 1;
                                Norm_err  = U.FitResult.err(3+U.ModelObj.nPixels:3+2*U.ModelObj.nPixels-1);
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
                                eta=((etaupper-etalower)./(chi2upper-chi2lower)).*DeltaChi2-chi2upper.*((etaupper-etalower)./(chi2upper-chi2lower))+etaupper;
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
                        save(savename,'eta','mnuSq','mnuSq_err','E0','E0_err','Bkg','Bkg_err','Norm','Norm_err');
                        if strcmp(Plot,'ON')
                            switch RunNr
                                case 1
                                    relic_global('eta',eta,'Params',obj.Params,'fitPar',fitPar,'E0',TwinBias_Q,'Syst',Syst);
                                case 10
                                    relic_global('eta',eta,'Params',obj.Params,'fitPar',fitPar,'Init_Opt',Init_Opt,'E0',TwinBias_Q,'Syst',Syst);
                            end
                        end
                        sprintf('Final Result: eta = %g',eta)
                        obj.etaSensitivity = eta;
                    end
                end
            end
        end
       
        function Chi2Scan_Twin(obj,varargin)
            p=inputParser;
            p.addParameter('Recompute','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('RunList','KNM1',@(x)ischar(x));                          % KNM1 or KNM2
            p.addParameter('fitPar','mNu E0 Norm Bkg',@(x)ischar(x));
            p.addParameter('Syst','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('SystBudget',24,@(x)isfloat(x));
            p.addParameter('TBDISBias','',@(x)isfloat(x) || isempty(x));
            p.addParameter('pullFlag',3);
            p.addParameter('TwinBias_mnuSq',0,@(x)isfloat(x));
            p.addParameter('range',30,@(x)isfloat(x));
            p.addParameter('Netabins',10,@(x)isfloat(x));
            p.addParameter('etarange',10,@(x)isfloat(x));
            p.addParameter('etafactor',1.5,@(x)isfloat(x));                    % max(eta)=etafactor*10^etarange
            p.addParameter('mode','SCAN',@(x)ismember(x,{'SCAN','SEARCH'}));
            p.addParameter('Plot','ON',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('CheckErrors','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('mnuSqfix','',@(x)(isfloat(x) && x>=0) || isempty(x));
            %% =========== SEARCH mode settings =============
            p.addParameter('etalower',0,@(x)isfloat(x));
            p.addParameter('etaupper',1.5e10,@(x)isfloat(x));                  % initial upper and lower search bounds
            p.addParameter('minchi2',0,@(x)isfloat(x));
            p.addParameter('delta',0.1e9,@(x)isfloat(x));                      % amount by which to shift eta if fit fails
            p.addParameter('DeltaChi2',2.71,@(x)isfloat(x));
            p.parse(varargin{:});
            Recompute      = p.Results.Recompute;
            RunList        = p.Results.RunList;
            fitPar         = p.Results.fitPar;
            Syst           = p.Results.Syst;
            SystBudget     = p.Results.SystBudget;
            TBDISBias      = p.Results.TBDISBias;
            pullFlag       = p.Results.pullFlag;
            TwinBias_mnuSq = p.Results.TwinBias_mnuSq;
            range          = p.Results.range;
            Netabins       = p.Results.Netabins;
            etarange       = p.Results.etarange;
            etafactor      = p.Results.etafactor;
            mode           = p.Results.mode;
            Plot           = p.Results.Plot;
            CheckErrors    = p.Results.CheckErrors;
            mnuSqfix       = p.Results.mnuSqfix;
            etalower       = p.Results.etalower;
            etaupper       = p.Results.etaupper;
            minchi2        = p.Results.minchi2;
            delta          = p.Results.delta;
            DeltaChi2      = p.Results.DeltaChi2;
            
            if strcmp(Syst,'ON')
                Chi2opt='chi2CMShape';
                if strcmp(obj.Params,'KNM1')
                    NonPoissonScaleFactor=1.064;
                elseif strcmp(obj.Params,'KNM2_Prompt')
                    NonPoissonScaleFactor=1.112;
                end
            else
                Chi2opt='chi2Stat';
                NonPoissonScaleFactor=1;
            end
            
            fitter = 'minuit';

            if strcmp(RunList,'KNM1')
                U = MultiRunAnalysis('RunList',RunList,... % runlist defines which runs are analysed -> set MultiRunAnalysis.m -> function: GetRunList()
                    'chi2',Chi2opt,...                 % uncertainties: statistical or stat + systematic uncertainties
                    'DataType','Twin',...                 % can be 'Real' or 'Twin' -> Monte Carlo
                    'fixPar',fitPar,...                   % free Parameter!!
                    'RadiativeFlag','ON',...              % theoretical radiative corrections applied in model
                    'NonPoissonScaleFactor',NonPoissonScaleFactor,...     % background uncertainty are enhanced
                    'fitter',fitter,...
                    'minuitOpt','min ; minos',...         % technical fitting options (minuit)
                    'pullFlag',pullFlag,...
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
            elseif strcmp(RunList,'KNM2_Prompt')
                U = MultiRunAnalysis('RunList',RunList,... % runlist defines which runs are analysed -> set MultiRunAnalysis.m -> function: GetRunList()
                    'chi2',Chi2opt,...                 % uncertainties: statistical or stat + systematic uncertainties
                    'DataType','Twin',...                 % can be 'Real' or 'Twin' -> Monte Carlo
                    'fixPar',fitPar,...                   % free Parameter!!
                    'RadiativeFlag','ON',...              % theoretical radiative corrections applied in model
                    'NonPoissonScaleFactor',NonPoissonScaleFactor,...     % background uncertainty are enhanced
                    'fitter',fitter,...
                    'minuitOpt','min ; minos',...         % technical fitting options (minuit)
                    'pullFlag',pullFlag,...
                    'FSDFlag','KNM2',...          % final state distribution
                    'ELossFlag','KatrinT2A20',...            % energy loss function
                    'SysBudget',SystBudget,...                    % defines syst. uncertainties -> in GetSysErr.m;
                    'DopplerEffectFlag','FSD',...
                    'Twin_SameCDFlag','OFF',...
                    'Twin_SameIsotopFlag','OFF',...
                    'SynchrotronFlag','ON',...
                    'AngularTFFlag','ON',...
                    'TwinBias_Q',18573.7,...
                    'TwinBias_mnuSq',TwinBias_mnuSq,...
                    'FSD_Sigma',sqrt(0.0124+0.0025),...
                    'TwinBias_FSDSigma',sqrt(0.0124+0.0025),...
                    'BKG_PtSlope',3*1e-06,...
                    'TwinBias_BKG_PtSlope',3*1e-06);
            end

            U.exclDataStart = U.GetexclDataStart(range); % set region of interest
            if ~contains(fitPar,'mNu')
                U.ModelObj.mnuSq_i=mnuSqfix;
            end
            if strcmp(U.chi2,'chi2Stat')
                U.InitModelObj_Norm_BKG('Recompute','ON');
            end
            
            if ~isempty(TBDISBias)
                U.RunData.TBDIS = U.RunData.TBDIS + TBDISBias;
            end
            obj.M = U;
            
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
                savename=[matFilePath,sprintf('RelicChi2Scan_Twin_BiasmnuSq%g_Syst%s_range%g_%s_[0 %g]_%s',TwinBias_mnuSq,Syst,range,obj.Params,etafactor*10^etarange,fitPar)];
                if any(ismember(pullFlag,1))
                    savename = [savename,'_mnuSqPull'];
                end
                if DeltaChi2~=2.71
                    savename = [savename,sprintf('DeltaChi2_%g',DeltaChi2)];
                end
                if ~isempty(mnuSqfix)
                    savename = [savename,sprintf('_relicPeakPos_%g',mnuSqfix)];
                end
                savename=[savename,'.mat'];
                
                if exist(savename,'file') && strcmp(Recompute,'OFF')
                    if strcmp(Plot,'ON')
                        plotchi2scan(savename);
                    end
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
                       if strcmp(CheckErrors,'ON')
                           mnuSq_err(i) = obj.CorrectErr('Parameter','mNu','value',mnuSq(i),'eta',(i-1)*((etafactor*10^(etarange))/(Netabins-1)),'minchi2',Chi2(i),'factor',1+mnuSq_err(1)/mnuSq(1));
                           E0_err(i)    = obj.CorrectErr('Parameter','E0','value',E0(i),'eta',(i-1)*((etafactor*10^(etarange))/(Netabins-1)),'minchi2',Chi2(i),'factor',1+E0_err(1)/E0(1));
                           Bkg_err(i)   = obj.CorrectErr('Parameter','Bkg','value',Bkg(i),'eta',(i-1)*((etafactor*10^(etarange))/(Netabins-1)),'minchi2',Chi2(i),'factor',1+Bkg_err(1)/Bkg(1));
                       end
                    end

                    save(savename,'Chi2','Netabins','etafactor','etarange','mnuSq','mnuSq_err','E0','E0_err','Bkg','Bkg_err','Norm','Norm_err');
                    if strcmp(Plot,'ON')
                        plotchi2scan(savename);
                    end
                end
            end

            if strcmp(mode,'SEARCH')
                
                matFilePath = [getenv('SamakPath'),sprintf('RelicNuBkg/UpperLimits/')];
                if isempty(TBDISBias)
                    savename=[matFilePath,sprintf('RelicLimit_Twin_BiasmnuSq%g_Syst%s_range%g_%s_%s',TwinBias_mnuSq,Syst,range,obj.Params,fitPar)];
                else
                    savename=[matFilePath,sprintf('RelicLimit_Twin_BiasmnuSq%g_Syst%s_TBDISBias_range%g_%s_%s',TwinBias_mnuSq,Syst,range,obj.Params,fitPar)];
                end
                if any(ismember(pullFlag,1))
                    savename = [savename,'_mnuSqPull'];
                end
                if DeltaChi2~=2.71
                    savename = [savename,sprintf('DeltaChi2_%g',DeltaChi2)];
                end
                if ~isempty(mnuSqfix)
                    savename = [savename,sprintf('_relicPeakPos_%g',mnuSqfix)];
                end
                savename = [savename,'.mat'];
                
                if exist(savename,'file') && strcmp(Recompute,'OFF')
                    load(savename,'eta');
                    if strcmp(Plot,'ON')
                        relic_global_twin('eta',eta,'RunList',obj.Params,'fitPar',fitPar,'E0',U.TwinBias_Q,'mnuSq',U.TwinBias_mnuSq,'Syst',Syst);
                    end
                    sprintf('Final Result: eta = %g',eta)
                    obj.etaSensitivity = eta;
                else
                    if strcmp(RunList,'KNM1')
                        F = MultiRunAnalysis('RunList',RunList,... % runlist defines which runs are analysed -> set MultiRunAnalysis.m -> function: GetRunList()
                             'chi2',Chi2opt,...                 % uncertainties: statistical or stat + systematic uncertainties
                             'DataType','Twin',...                 % can be 'Real' or 'Twin' -> Monte Carlo
                             'fixPar',fitPar,...                   % free Parameter!!
                             'RadiativeFlag','ON',...              % theoretical radiative corrections applied in model
                             'NonPoissonScaleFactor',NonPoissonScaleFactor,...     % background uncertainty are enhanced
                             'fitter',fitter,...
                             'minuitOpt','min ; minos',...         % technical fitting options (minuit)
                             'pullFlag',pullFlag,...
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
                    elseif strcmp(RunList,'KNM2_Prompt')
                        F = MultiRunAnalysis('RunList',RunList,... % runlist defines which runs are analysed -> set MultiRunAnalysis.m -> function: GetRunList()
                            'chi2',Chi2opt,...                 % uncertainties: statistical or stat + systematic uncertainties
                            'DataType','Twin',...                 % can be 'Real' or 'Twin' -> Monte Carlo
                            'fixPar',fitPar,...                   % free Parameter!!
                            'RadiativeFlag','ON',...              % theoretical radiative corrections applied in model
                            'NonPoissonScaleFactor',NonPoissonScaleFactor,...     % background uncertainty are enhanced
                            'fitter',fitter,...
                            'minuitOpt','min ; minos',...         % technical fitting options (minuit)
                            'pullFlag',pullFlag,...
                            'FSDFlag','KNM2',...          % final state distribution
                            'ELossFlag','KatrinT2A20',...            % energy loss function
                            'SysBudget',SystBudget,...                    % defines syst. uncertainties -> in GetSysErr.m;
                            'DopplerEffectFlag','FSD',...
                            'Twin_SameCDFlag','OFF',...
                            'Twin_SameIsotopFlag','OFF',...
                            'SynchrotronFlag','ON',...
                            'AngularTFFlag','ON',...
                            'TwinBias_Q',18573.7,...
                            'TwinBias_mnuSq',TwinBias_mnuSq,...
                            'FSD_Sigma',sqrt(0.0124+0.0025),...
                            'TwinBias_FSDSigma',sqrt(0.0124+0.0025),...
                            'BKG_PtSlope',3*1e-06,...
                            'TwinBias_BKG_PtSlope',3*1e-06);
                    end

                    F.exclDataStart = F.GetexclDataStart(range); % set region of interest
                    if ~contains(fitPar,'mNu')
                        F.ModelObj.mnuSq_i=mnuSqfix;
                    end
                    if strcmp(F.chi2,'chi2Stat')
                        F.InitModelObj_Norm_BKG('Recompute','ON');
                    end
                    
                    if ~isempty(TBDISBias)
                        F.RunData.TBDIS = F.RunData.TBDIS + sqrt(diag(obj.M.FitCM));
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
                    chi2upper = U.FitResult.chi2min;
                    chi2lower = F.FitResult.chi2min;
                    if etalower==0
                        minchi2   = chi2lower;
                    end
                    chi2      = chi2upper;
                    mnuSq     = U.ModelObj.mnuSq_i+U.FitResult.par(1);
                    mnuSq_err = U.FitResult.err(1);
                    E0        = U.ModelObj.Q_i+U.FitResult.par(2);
                    E0_err    = U.FitResult.err(2);
                    Bkg       = U.ModelObj.BKG_RateSec_i+U.FitResult.par(3);
                    Bkg_err   = U.FitResult.err(3);
                    Norm      = U.FitResult.par(3+U.ModelObj.nPixels:3+2*U.ModelObj.nPixels-1) + 1;
                    Norm_err  = U.FitResult.err(3+U.ModelObj.nPixels:3+2*U.ModelObj.nPixels-1);
                    if chi2<DeltaChi2
                        sprintf('chi^2 too small! Set larger initial eta.')
                    elseif chi2>50
                        sprintf('Minos failed! Vary initial eta.')
                    else
                        if chi2<(DeltaChi2+minchi2)
                            eta = etalower;
                            chi2=chi2lower;
                        end
                        while abs(chi2-minchi2-DeltaChi2)>0.01
                            if chi2>(DeltaChi2+minchi2)
                                etaupper=eta;
                                chi2upper=chi2;
                                eta=((etaupper-etalower)./(chi2upper-chi2lower)).*(DeltaChi2+minchi2)-chi2upper.*((etaupper-etalower)./(chi2upper-chi2lower))+etaupper;
                                U.ModelObj.eta_i=eta;
                                U.ModelObj.eta=eta;
                                U.ModelObj.ComputeNormFactorTBDDS;
                                U.ModelObj.ComputeTBDDS;
                                U.ModelObj.ComputeTBDIS;
                                U.Fit;
                                mnuSq     = U.ModelObj.mnuSq_i+U.FitResult.par(1);
                                mnuSq_err = U.FitResult.err(1);
                                E0        = U.ModelObj.Q_i+U.FitResult.par(2);
                                E0_err    = U.FitResult.err(2);
                                Bkg       = U.ModelObj.BKG_RateSec_i+U.FitResult.par(3);
                                Bkg_err   = U.FitResult.err(3);
                                Norm      = U.FitResult.par(3+U.ModelObj.nPixels:3+2*U.ModelObj.nPixels-1) + 1;
                                Norm_err  = U.FitResult.err(3+U.ModelObj.nPixels:3+2*U.ModelObj.nPixels-1);
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
                                eta=((etaupper-etalower)./(chi2upper-chi2lower)).*(DeltaChi2+minchi2)-chi2upper.*((etaupper-etalower)./(chi2upper-chi2lower))+etaupper;
                                F.ModelObj.eta_i=eta;
                                F.ModelObj.eta=eta;
                                F.ModelObj.ComputeNormFactorTBDDS;
                                F.ModelObj.ComputeTBDDS;
                                F.ModelObj.ComputeTBDIS;
                                F.Fit;
                                mnuSq     = F.ModelObj.mnuSq_i+F.FitResult.par(1);
                                mnuSq_err = F.FitResult.err(1);
                                E0        = F.ModelObj.Q_i+F.FitResult.par(2);
                                E0_err    = F.FitResult.err(2);
                                Bkg       = F.ModelObj.BKG_RateSec_i+F.FitResult.par(3);
                                Bkg_err   = F.FitResult.err(3);
                                Norm      = F.FitResult.par(3+F.ModelObj.nPixels:3+2*F.ModelObj.nPixels-1) + 1;
                                Norm_err  = F.FitResult.err(3+F.ModelObj.nPixels:3+2*F.ModelObj.nPixels-1);
                                if F.FitResult.chi2min>50
                                    sprintf('Fit failed! Retrying...')
                                    etalower = etalower - delta;
                                    eta = etalower;
                                else
                                    chi2=F.FitResult.chi2min;
                                end
                            end
                        end
                        save(savename,'eta','mnuSq','mnuSq_err','E0','E0_err','Bkg','Bkg_err','Norm','Norm_err');
                        if strcmp(Plot,'ON')
                            relic_global_twin('eta',eta,'RunList',obj.Params,'fitPar',fitPar,'E0',U.TwinBias_Q,'mnuSq',U.TwinBias_mnuSq,'Syst',Syst);
                        end
                        sprintf('Final Result: eta = %g',eta)
                        obj.etaSensitivity = eta;
                    end
                end
            end
        end
        function Chi2Scan_2D(obj,varargin)
            p=inputParser;
            p.addParameter('Recompute','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('RunList','KNM1',@(x)ischar(x));                          % KNM1 or KNM2
            p.addParameter('Syst','ON',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('SystBudget',24,@(x)isfloat(x));
            p.addParameter('DataType','Twin',@(x)ismember(x,{'Twin','Real'}));
            p.addParameter('TwinBias_mnuSq',0,@(x)isfloat(x));
            p.addParameter('range',40,@(x)isfloat(x));
            p.addParameter('pullFlag',3);
            p.addParameter('Netabins',20,@(x)isfloat(x));
            p.addParameter('DeltaEta',3e11,@(x)isfloat(x) && x>0);
            p.addParameter('Nmnubins',20,@(x)isfloat(x));
            p.addParameter('Deltamnu',1,@(x)isfloat(x) && x>0);
            p.addParameter('Plot','ON',@(x)ismember(x,{'ON','OFF'}));
            p.parse(varargin{:});
            Recompute      = p.Results.Recompute;
            RunList        = p.Results.RunList;
            Syst           = p.Results.Syst;
            SystBudget     = p.Results.SystBudget;
            DataType       = p.Results.DataType;
            TwinBias_mnuSq = p.Results.TwinBias_mnuSq;
            range          = p.Results.range;
            pullFlag       = p.Results.pullFlag;
            Netabins       = p.Results.Netabins;
            DeltaEta       = p.Results.DeltaEta;
            Plot           = p.Results.Plot;
            Nmnubins       = p.Results.Nmnubins;
            Deltamnu       = p.Results.Deltamnu;
            
            matFilePath = [getenv('SamakPath'),sprintf('RelicNuBkg/Chi2Scans/')];
            savename=[matFilePath,sprintf('RelicChi2Scan_Twin_BiasmnuSq%g_Syst%s_range%g_%s_2D',TwinBias_mnuSq,Syst,range,obj.Params)];
            if strcmp(DataType,'Real')
                savename=[savename,'_Real'];
            end
            savename=[savename,'.mat'];
                
            if exist(savename,'file') && strcmp(Recompute,'OFF') && strcmp(Plot,'ON')
                plotchi2scan(savename);
                load(savename);
                DeltaChi2 = abs(Chi2 - 4.61);
                GlobalChi2Min = min(min(DeltaChi2));
                [a,b] = find(DeltaChi2==GlobalChi2Min);
                Uppermnu = mnuScanPoints(a);
                Uppereta = etaScanPoints(b);
                sprintf('m_{\\nu} = %.2g \n \\eta = %.2g \n \\Delta \\chi^{2} = %.3g',Uppermnu,Uppereta,Chi2(a,b))
                if (a==Nmnubins || b==Netabins) && Chi2(a,b)<4.61 && Chi2(a,b)==max(max(Chi2))
                    sprintf('Upper limit out of range! Scan range should be extended!')
                end
            else
            
                if strcmp(Syst,'ON')
                    Chi2opt='chi2CMShape';
                    switch RunList
                        case 'KNM1'
                            NonPoissonScaleFactor=1.064;
                        case 'KNM2_Prompt'
                            NonPoissonScaleFactor=1.112;
                    end
                else
                    Chi2opt='chi2Stat';
                    NonPoissonScaleFactor=1;
                end

                fitter = 'minuit';

                if strcmp(RunList,'KNM1')
                    U = MultiRunAnalysis('RunList',RunList,... % runlist defines which runs are analysed -> set MultiRunAnalysis.m -> function: GetRunList()
                        'chi2',Chi2opt,...                 % uncertainties: statistical or stat + systematic uncertainties
                        'DataType',DataType,...                 % can be 'Real' or 'Twin' -> Monte Carlo
                        'fixPar','E0 Norm Bkg',...                   % free Parameter!!
                        'RadiativeFlag','ON',...              % theoretical radiative corrections applied in model
                        'NonPoissonScaleFactor',NonPoissonScaleFactor,...     % background uncertainty are enhanced
                        'fitter',fitter,...
                        'minuitOpt','min ; minos',...         % technical fitting options (minuit)
                        'pullFlag',pullFlag,...
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
                    W = MultiRunAnalysis('RunList',RunList,... % runlist defines which runs are analysed -> set MultiRunAnalysis.m -> function: GetRunList()
                        'chi2',Chi2opt,...                 % uncertainties: statistical or stat + systematic uncertainties
                        'DataType',DataType,...                 % can be 'Real' or 'Twin' -> Monte Carlo
                        'fixPar','mNu E0 Norm Bkg eta',...                   % free Parameter!!
                        'RadiativeFlag','ON',...              % theoretical radiative corrections applied in model
                        'NonPoissonScaleFactor',NonPoissonScaleFactor,...     % background uncertainty are enhanced
                        'fitter',fitter,...
                        'minuitOpt','min ; minos',...         % technical fitting options (minuit)
                        'pullFlag',pullFlag,...
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
                elseif strcmp(RunList,'KNM2_Prompt')
                    U = MultiRunAnalysis('RunList',RunList,... % runlist defines which runs are analysed -> set MultiRunAnalysis.m -> function: GetRunList()
                        'chi2',Chi2opt,...                 % uncertainties: statistical or stat + systematic uncertainties
                        'DataType',DataType,...                 % can be 'Real' or 'Twin' -> Monte Carlo
                        'fixPar','E0 Norm Bkg',...                   % free Parameter!!
                        'RadiativeFlag','ON',...              % theoretical radiative corrections applied in model
                        'NonPoissonScaleFactor',NonPoissonScaleFactor,...     % background uncertainty are enhanced
                        'fitter',fitter,...
                        'minuitOpt','min ; minos',...         % technical fitting options (minuit)
                        'pullFlag',pullFlag,...
                        'FSDFlag','KNM2',...          % final state distribution
                        'ELossFlag','KatrinT2A20',...            % energy loss function
                        'SysBudget',SystBudget,...                    % defines syst. uncertainties -> in GetSysErr.m;
                        'DopplerEffectFlag','FSD',...
                        'Twin_SameCDFlag','OFF',...
                        'Twin_SameIsotopFlag','OFF',...
                        'SynchrotronFlag','ON',...
                        'AngularTFFlag','ON',...
                        'TwinBias_Q',18573.7,...
                        'TwinBias_mnuSq',TwinBias_mnuSq,...
                        'FSD_Sigma',sqrt(0.0124+0.0025),...
                        'TwinBias_FSDSigma',sqrt(0.0124+0.0025),...
                        'BKG_PtSlope',3*1e-06,...
                        'TwinBias_BKG_PtSlope',3*1e-06);
                    W = MultiRunAnalysis('RunList',RunList,... % runlist defines which runs are analysed -> set MultiRunAnalysis.m -> function: GetRunList()
                        'chi2',Chi2opt,...                 % uncertainties: statistical or stat + systematic uncertainties
                        'DataType',DataType,...                 % can be 'Real' or 'Twin' -> Monte Carlo
                        'fixPar','mNu E0 Norm Bkg eta',...                   % free Parameter!!
                        'RadiativeFlag','ON',...              % theoretical radiative corrections applied in model
                        'NonPoissonScaleFactor',NonPoissonScaleFactor,...     % background uncertainty are enhanced
                        'fitter',fitter,...
                        'minuitOpt','min ; minos',...         % technical fitting options (minuit)
                        'pullFlag',pullFlag,...
                        'FSDFlag','KNM2',...          % final state distribution
                        'ELossFlag','KatrinT2A20',...            % energy loss function
                        'SysBudget',SystBudget,...                    % defines syst. uncertainties -> in GetSysErr.m;
                        'DopplerEffectFlag','FSD',...
                        'Twin_SameCDFlag','OFF',...
                        'Twin_SameIsotopFlag','OFF',...
                        'SynchrotronFlag','ON',...
                        'AngularTFFlag','ON',...
                        'TwinBias_Q',18573.7,...
                        'TwinBias_mnuSq',TwinBias_mnuSq,...
                        'FSD_Sigma',sqrt(0.0124+0.0025),...
                        'TwinBias_FSDSigma',sqrt(0.0124+0.0025),...
                        'BKG_PtSlope',3*1e-06,...
                        'TwinBias_BKG_PtSlope',3*1e-06);
                end

                U.exclDataStart = U.GetexclDataStart(range); % set region of interest
                if strcmp(U.chi2,'chi2Stat')
                    U.InitModelObj_Norm_BKG('Recompute','ON');
                end
                
                W.exclDataStart = U.GetexclDataStart(range);
                W.Fit;
                mnuSq_best    = W.FitResult.par(1);
                eta_best      = W.FitResult.par(18)*1e10;
                GlobalChi2Min = W.FitResult.chi2min;

                mnuScanPoints = linspace(mnuSq_best-Deltamnu,mnuSq_best+Deltamnu,Nmnubins);
                etaScanPoints = linspace(eta_best-DeltaEta,eta_best+DeltaEta,Netabins);
                Chi2 = ones(Nmnubins,Netabins);
                E0   = ones(Nmnubins,Netabins);
                Norm = ones(Nmnubins,Netabins);
                Bkg  = ones(Nmnubins,Netabins);

                for i=1:Nmnubins
                    U.ModelObj.mnuSq_i = mnuScanPoints(i);
                    for j=1:Netabins
                        U.ModelObj.eta_i = etaScanPoints(j);
                        U.ModelObj.eta = etaScanPoints(j);
                        U.ModelObj.ComputeNormFactorTBDDS;
                        U.ModelObj.ComputeTBDDS;
                        U.ModelObj.ComputeTBDIS;
                        U.Fit;
                        Chi2(i,j) = U.FitResult.chi2min;
                        E0(i,j)   = U.ModelObj.Q_i+U.FitResult.par(2);
                        Bkg(i,j)  = U.ModelObj.BKG_RateSec_i+U.FitResult.par(3);
                        Norm(i,j) = U.FitResult.par(3+U.ModelObj.nPixels:3+2*U.ModelObj.nPixels-1) + 1;
                    end
                end
                save(savename,'mnuScanPoints','etaScanPoints','Chi2','GlobalChi2Min','E0','Bkg','Norm');
                DeltaChi2 = abs(Chi2 - GlobalChi2Min - 4.61);
                GlobalChi2Min = min(min(DeltaChi2));
                [a,b] = find(DeltaChi2==GlobalChi2Min);
                Uppermnu = mnuScanPoints(a);
                Uppereta = etaScanPoints(b);
                sprintf('m_{\\nu} = %.2g \n \\eta = %.2g \n \\Delta \\chi^{2} = %.3g',Uppermnu,Uppereta,Chi2(a,b))
                if (a==Nmnubins || b==Netabins) && Chi2(a,b)<4.61 && Chi2(a,b)==max(max(Chi2))
                    sprintf('Upper limit out of range! Scan range should be extended!')
                end
                plotchi2scan(savename);
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
            p.addParameter('pullFlag',3);
            p.addParameter('TwinBias_mnuSq',1,@(x)isfloat(x));
            p.addParameter('NetaBins',10,@(x)isfloat(x));
            p.addParameter('etarange',11,@(x)isfloat(x));
            p.addParameter('etafactor',5,@(x)isfloat(x));
            p.addParameter('Recompute','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('Plot','ON',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('DeltaChi2',2.71,@(x)isfloat(x));
            p.addParameter('CheckErrors','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('mnuSqfix','',@(x)(isfloat(x) && x>=0) || isempty(x));
            p.parse(varargin{:});

            range          = p.Results.range;
            fitPar         = p.Results.fitPar;
            Syst           = p.Results.Syst;
            SystBudget     = p.Results.SystBudget;
            pullFlag       = p.Results.pullFlag;
            TwinBias_mnuSq = p.Results.TwinBias_mnuSq;
            NetaBins       = p.Results.NetaBins;
            etarange       = p.Results.etarange;
            etafactor      = p.Results.etafactor;
            Recompute      = p.Results.Recompute;
            Plot           = p.Results.Plot;
            DeltaChi2      = p.Results.DeltaChi2;
            CheckErrors    = p.Results.CheckErrors;
            mnuSqfix       = p.Results.mnuSqfix;

            
            obj.Chi2Scan_Twin('Recompute',Recompute,...
                'range',range,...
                'RunList',obj.Params,...
                'fitPar',fitPar,...
                'Syst',Syst,...
                'SystBudget',SystBudget,...
                'pullFlag',pullFlag,...
                'TwinBias_mnuSq',TwinBias_mnuSq,...
                'Netabins',NetaBins,...
                'etarange',etarange,...
                'etafactor',etafactor,...
                'mode','SCAN',...
                'Plot',Plot,...
                'DeltaChi2',DeltaChi2,...
                'CheckErrors',CheckErrors,...
                'mnuSqfix',mnuSqfix);

            matFilePath = [getenv('SamakPath'),sprintf('RelicNuBkg/Chi2Scans/')];
            savename=[matFilePath,sprintf('RelicChi2Scan_Twin_BiasmnuSq%g_Syst%s_range%g_%s_[0 %g]_%s',TwinBias_mnuSq,Syst,range,obj.Params,etafactor*10^etarange,fitPar)];
            if any(ismember(pullFlag,1))
                savename = [savename,'_mnuSqPull'];
            end
            if DeltaChi2~=2.71
                savename = [savename,sprintf('DeltaChi2_%g',DeltaChi2)];
            end
            if ~isempty(mnuSqfix)
                savename = [savename,sprintf('_relicPeakPos_%g',mnuSqfix)];
            end
            savename=[savename,'.mat'];

            load(savename);

            if (Chi2(end)-Chi2(1))<DeltaChi2
                sprintf('Increase etafactor or etarange')
            else
                etavalues=(0:(Netabins-1))*((etafactor*10^(etarange))/(Netabins-1));
                m=find(Chi2<(DeltaChi2+Chi2(1)));
                n=find(Chi2>(DeltaChi2+Chi2(1)));
                etalower=etavalues(m(end));
                etaupper=etavalues(n(1));

                obj.Chi2Scan_Twin('Recompute',Recompute,...
                    'range',range,...
                    'RunList',obj.Params,...
                    'fitPar',fitPar,...
                    'Syst',Syst,...
                    'SystBudget',SystBudget,...
                    'pullFlag',pullFlag,...
                    'TwinBias_mnuSq',TwinBias_mnuSq,...
                    'etalower',etalower,...
                    'etaupper',etaupper,...
                    'minchi2',Chi2(1),...
                    'mode','SEARCH',...
                    'Plot',Plot,...
                    'DeltaChi2',DeltaChi2,...
                    'mnuSqfix',mnuSqfix);
            end
       end
        
       function Chi2Fake(obj,varargin)
            p = inputParser;
            p.addParameter('range',50,@(x)isfloat(x));
            p.addParameter('RunNr',1,@(x)isfloat(x));
            p.addParameter('NetaBins',10,@(x)isfloat(x));
            p.addParameter('etarange',10,@(x)isfloat(x));
            p.addParameter('etafactor',3,@(x)isfloat(x));
            p.addParameter('fitPar','mNu E0 Norm Bkg',@(x)ischar(x));
            p.addParameter('Init_Opt','',@(x)iscell(x) || isempty(x));        % use these options by switching to RunNr 10
            p.addParameter('Syst','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('pullFlag',3);
            p.addParameter('Recompute','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('Plot','ON',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('DeltaChi2',2.71,@(x)isfloat(x));
            
            p.parse(varargin{:});
            range     = p.Results.range;
            RunNr     = p.Results.RunNr;
            NetaBins  = p.Results.NetaBins;
            etarange  = p.Results.etarange;
            etafactor = p.Results.etafactor;
            fitPar    = p.Results.fitPar;
            Init_Opt  = p.Results.Init_Opt;
            Syst      = p.Results.Syst;
            pullFlag  = p.Results.pullFlag;
            Recompute = p.Results.Recompute;
            Plot      = p.Results.Plot;
            DeltaChi2 = p.Results.DeltaChi2;

            obj.Chi2Scan_Fake('Recompute',Recompute,...
                'range',range,...
                'RunNr',RunNr,...
                'fitPar',fitPar,...
                'Syst',Syst,...
                'pullFlag',pullFlag,...
                'Init_Opt',Init_Opt,...
                'Netabins',NetaBins,...
                'etarange',etarange,...
                'etafactor',etafactor,...
                'mode','SCAN',...
                'Plot',Plot,...
                'DeltaChi2',DeltaChi2);

            matFilePath = [getenv('SamakPath'),sprintf('RelicNuBkg/Chi2Scans/')];
            if RunNr==1
                savename=[matFilePath,sprintf('RelicChi2Scan_Fake_Syst%s_range%g_%s_[0 %g]_%s.mat',Syst,range,obj.Params,etafactor*10^etarange,fitPar)];
            elseif RunNr==10
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

            load(savename);

            if Chi2(end)<DeltaChi2
                sprintf('Increase etafactor or etarange')
            else
                etavalues=(0:(Netabins-1))*((etafactor*10^(etarange))/(Netabins-1));
                m=find(Chi2<DeltaChi2);
                n=find(Chi2>DeltaChi2);
                etalower=etavalues(m(end));
                etaupper=etavalues(n(1));

                obj.Chi2Scan_Fake('Recompute',Recompute,...
                    'range',range,...
                    'RunNr',RunNr,...
                    'fitPar',fitPar,...
                    'Syst',Syst,...
                    'pullFlag',pullFlag,...
                    'Init_Opt',Init_Opt,...
                    'etalower',etalower,...
                    'etaupper',etaupper,...
                    'mode','SEARCH',...
                    'Plot',Plot,...
                    'DeltaChi2',DeltaChi2);
            end
       end
       function [mNu,eta] = EtaFit(obj,varargin)
           p = inputParser();
           p.addParameter('NmNuBins',21,@(x)isfloat(x));
           p.addParameter('MinmNu',0,@(x)isfloat(x));
           p.addParameter('MaxmNu',2,@(x)isfloat(x));
           p.addParameter('Syst','ON',@(x)ismember(x,{'ON','OFF'}));
           p.addParameter('SysBudget',24,@(x)isfloat(x));
           p.addParameter('DataType','Twin',@(x)ismember(x,{'Twin','Real'}));
           p.addParameter('Recompute','OFF',@(x)ismember(x,{'ON','OFF'}));
           p.addParameter('Mode','MultimnuSq',@(x)ismember(x,{'MultimnuSq','SinglemnuSq'}));
           p.addParameter('fitPar','mNu E0 Norm Bkg eta',@(x)ischar(x));
           p.addParameter('TwinmnuSq',0,@(x)isfloat(x));
           p.addParameter('DeltaChi2',1,@(x)isfloat(x));
           p.addParameter('CorrectErr','ON',@(x)ismember(x,{'ON','OFF'}));
           p.parse(varargin{:});
           NmNuBins   = p.Results.NmNuBins;
           MinmNu     = p.Results.MinmNu;
           MaxmNu     = p.Results.MaxmNu;
           Syst       = p.Results.Syst;
           SysBudget  = p.Results.SysBudget;
           DataType   = p.Results.DataType;
           Recompute  = p.Results.Recompute;
           Mode       = p.Results.Mode;
           fitPar     = p.Results.fitPar;
           TwinmnuSq  = p.Results.TwinmnuSq;
           DeltaChi2  = p.Results.DeltaChi2;
           CorrectErr = p.Results.CorrectErr;
           
           matFilePath = [getenv('SamakPath'),sprintf('RelicNuBkg/UpperLimits/')];
           if MinmNu~=0
               savename    = [matFilePath,sprintf('EtaFit_mNu_%g_%g_Syst%s_DataType%s',MinmNu,MaxmNu,Syst,DataType)];
           else
               savename    = [matFilePath,sprintf('EtaFit_mNu0_%g_Syst%s_DataType%s',MaxmNu,Syst,DataType)];
           end
           if strcmp(Mode,'SinglemnuSq')
               savename    = [savename,sprintf('SinglemnuSq_TwinBias_mnuSq%g',TwinmnuSq)];
           end
           if ~strcmp(fitPar,'mNu E0 Norm Bkg eta')
               savename    = [savename,sprintf('_%s',fitPar)];
           end
           if DeltaChi2~=2.71
               savename    = [savename,sprintf('DeltaChi2_%i',DeltaChi2)];
           end
           savename        = [savename,'.mat'];
           
           if ~exist(savename,'file') || strcmp(Recompute,'ON')
               eta       = zeros(1,NmNuBins);
               Chi2      = zeros(1,NmNuBins);
               mnuSq     = zeros(1,NmNuBins);
               mnuSq_err = zeros(1,NmNuBins);
               Norm      = zeros(1,NmNuBins);
               Norm_err  = zeros(1,NmNuBins);
               E0        = zeros(1,NmNuBins);
               E0_err    = zeros(1,NmNuBins);
               Bkg       = zeros(1,NmNuBins);
               Bkg_err   = zeros(1,NmNuBins);
               etaFit    = eta;
               if strcmp(Syst,'ON')
                   chi2opt = 'chi2CMShape';
                   if strcmp(obj.Params,'KNM1')
                       SysBudget=24;
                       NPfac    = 1.064;
                   elseif strcmp(obj.Params,'KNM2_Prompt')
                       SysBudget=40;
                       NPfac    =1.1102;
                   elseif strcmp(obj.Params,'TDR')
                       SysBudget=67;
                       NPfac    =1;
                   end
               else
                   chi2opt = 'chi2Stat';
                   NPfac   = 1;
               end

               if strcmp(Mode,'MultimnuSq')
                   mNu  = linspace(MinmNu,MaxmNu,NmNuBins);
                   for i=1:NmNuBins
                       if strcmp(obj.Params,'KNM1')
                            U = MultiRunAnalysis('RunList',obj.Params,... % runlist defines which runs are analysed -> set MultiRunAnalysis.m -> function: GetRunList()
                                'chi2',chi2opt,...                 % uncertainties: statistical or stat + systematic uncertainties
                                'DataType',DataType,...                 % can be 'Real' or 'Twin' -> Monte Carlo
                                'fixPar',fitPar,...        % free Parameter!!
                                'RadiativeFlag','ON',...              % theoretical radiative corrections applied in model
                                'NonPoissonScaleFactor',NPfac,...     % background uncertainty are enhanced
                                'fitter','minuit',...                 % minuit standard, matlab to be tried
                                'minuitOpt','min ; minos',...         % technical fitting options (minuit)
                                'FSDFlag','SibilleFull',...           % final state distribution
                                'ELossFlag','KatrinT2',...            % energy loss function
                                'SysBudget',SysBudget,...                    % defines syst. uncertainties -> in GetSysErr.m;
                                'DopplerEffectFlag','FSD',...
                                'Twin_SameCDFlag','OFF',...
                                'Twin_SameIsotopFlag','OFF',...
                                'SynchrotronFlag','ON',...
                                'AngularTFFlag','OFF',...
                                'TwinBias_Q',18573.73,...
                                'TwinBias_mnuSq',mNu(i));
                       elseif strcmp(obj.Params,'KNM2_Prompt')
                           U = MultiRunAnalysis('RunList',obj.Params,... % runlist defines which runs are analysed -> set MultiRunAnalysis.m -> function: GetRunList()
                                'chi2',chi2opt,...                 % uncertainties: statistical or stat + systematic uncertainties
                                'DataType',DataType,...                 % can be 'Real' or 'Twin' -> Monte Carlo
                                'fixPar',fitPar,...                   % free Parameter!!
                                'RadiativeFlag','ON',...              % theoretical radiative corrections applied in model
                                'NonPoissonScaleFactor',NPfac,...     % background uncertainty are enhanced
                                'fitter','minuit',...
                                'minuitOpt','min ; minos',...         % technical fitting options (minuit)
                                'FSDFlag','KNM2',...          % final state distribution
                                'ELossFlag','KatrinT2A20',...            % energy loss function
                                'SysBudget',SysBudget,...                    % defines syst. uncertainties -> in GetSysErr.m;
                                'DopplerEffectFlag','FSD',...
                                'Twin_SameCDFlag','OFF',...
                                'Twin_SameIsotopFlag','OFF',...
                                'SynchrotronFlag','ON',...
                                'AngularTFFlag','ON',...
                                'TwinBias_Q',18573.7,...
                                'TwinBias_mnuSq',mNu(i),...
                                'FSD_Sigma',sqrt(0.0124+0.0025),...
                                'TwinBias_FSDSigma',sqrt(0.0124+0.0025),...
                                'BKG_PtSlope',3*1e-06,...
                                'TwinBias_BKG_PtSlope',3*1e-06);
                       end

                        U.exclDataStart=U.GetexclDataStart(40);
                        U.ModelObj.mnuSq_i = U.TwinBias_mnuSq;
                        if strcmp(U.chi2,'chi2Stat')
                            U.InitModelObj_Norm_BKG('Recompute','ON');
                        end
                        obj.M = U;
                        U.Fit;
                        etaFit(i)=U.FitResult.par(18)*1e10;
                        if strcmp(CorrectErr,'ON') || DeltaChi2~=1
                            eta(i) = obj.CorrectErr('Parameter','eta',...
                                'value',U.FitResult.par(18)*1e10,...
                                'eta',U.FitResult.par(18)*1e10,...
                                'minchi2',U.FitResult.chi2min,...
                                'factor',(1+U.FitResult.err(18)/U.FitResult.par(18)),...
                                'fitPar',fitPar,...
                                'DeltaChi2',DeltaChi2);
                        else
                            eta(i) = U.FitResult.err(18)*1e10;
                        end
                        Chi2(i)      = U.FitResult.chi2min;
                        mnuSq(i)     = U.ModelObj.mnuSq_i+U.FitResult.par(1);
                        mnuSq_err(i) = U.FitResult.err(1);
                        E0(i)        = U.ModelObj.Q_i+U.FitResult.par(2);
                        E0_err(i)    = U.FitResult.err(2);
                        Bkg(i)       = U.ModelObj.BKG_RateSec_i+U.FitResult.par(3);
                        Bkg_err(i)   = U.FitResult.err(3);
                        Norm(i)      = U.FitResult.par(3+U.ModelObj.nPixels:3+2*U.ModelObj.nPixels-1) + 1;
                        Norm_err(i)  = U.FitResult.err(3+U.ModelObj.nPixels:3+2*U.ModelObj.nPixels-1);
                        sprintf('Final Result: eta = %g',eta(i))
                   end
               else
                   mNu  = linspace(0,MaxmNu,NmNuBins).^2;
                   for i=1:NmNuBins
                       if strcmp(obj.Params,'KNM1')
                            U = MultiRunAnalysis('RunList',obj.Params,... % runlist defines which runs are analysed -> set MultiRunAnalysis.m -> function: GetRunList()
                                'chi2',chi2opt,...                 % uncertainties: statistical or stat + systematic uncertainties
                                'DataType',DataType,...                 % can be 'Real' or 'Twin' -> Monte Carlo
                                'fixPar',fitPar,...        % free Parameter!!
                                'RadiativeFlag','ON',...              % theoretical radiative corrections applied in model
                                'NonPoissonScaleFactor',NPfac,...     % background uncertainty are enhanced
                                'fitter','minuit',...                 % minuit standard, matlab to be tried
                                'minuitOpt','min ; minos',...         % technical fitting options (minuit)
                                'FSDFlag','SibilleFull',...           % final state distribution
                                'ELossFlag','KatrinT2',...            % energy loss function
                                'SysBudget',SysBudget,...                    % defines syst. uncertainties -> in GetSysErr.m;
                                'DopplerEffectFlag','FSD',...
                                'Twin_SameCDFlag','OFF',...
                                'Twin_SameIsotopFlag','OFF',...
                                'SynchrotronFlag','ON',...
                                'AngularTFFlag','OFF',...
                                'TwinBias_Q',18573.73,...
                                'TwinBias_mnuSq',TwinmnuSq);
                       elseif strcmp(obj.Params,'KNM2_Prompt')
                           U = MultiRunAnalysis('RunList',obj.Params,... % runlist defines which runs are analysed -> set MultiRunAnalysis.m -> function: GetRunList()
                                'chi2',chi2opt,...                 % uncertainties: statistical or stat + systematic uncertainties
                                'DataType',DataType,...                 % can be 'Real' or 'Twin' -> Monte Carlo
                                'fixPar',fitPar,...                   % free Parameter!!
                                'RadiativeFlag','ON',...              % theoretical radiative corrections applied in model
                                'NonPoissonScaleFactor',NPfac,...     % background uncertainty are enhanced
                                'fitter','minuit',...
                                'minuitOpt','min ; minos',...         % technical fitting options (minuit)
                                'FSDFlag','KNM2',...          % final state distribution
                                'ELossFlag','KatrinT2A20',...            % energy loss function
                                'SysBudget',SysBudget,...                    % defines syst. uncertainties -> in GetSysErr.m;
                                'DopplerEffectFlag','FSD',...
                                'Twin_SameCDFlag','OFF',...
                                'Twin_SameIsotopFlag','OFF',...
                                'SynchrotronFlag','ON',...
                                'AngularTFFlag','ON',...
                                'TwinBias_Q',18573.7,...
                                'TwinBias_mnuSq',TwinmnuSq,...
                                'FSD_Sigma',sqrt(0.0124+0.0025),...
                                'TwinBias_FSDSigma',sqrt(0.0124+0.0025),...
                                'BKG_PtSlope',3*1e-06,...
                                'TwinBias_BKG_PtSlope',3*1e-06);
                       end

                        U.exclDataStart=U.GetexclDataStart(40);
                        U.ModelObj.mnuSq_i = mNu(i);
                        if strcmp(U.chi2,'chi2Stat')
                            U.InitModelObj_Norm_BKG('Recompute','ON');
                        end
                        obj.M = U;
                        U.Fit;
                        etaFit(i)=U.FitResult.par(18)*1e10;
                        if strcmp(CorrectErr,'ON') || DeltaChi2~=1
                            eta(i) = obj.CorrectErr('Parameter','eta',...
                                'value',U.FitResult.par(18)*1e10,...
                                'eta',U.FitResult.par(18)*1e10,...
                                'minchi2',U.FitResult.chi2min,...
                                'factor',(1+U.FitResult.err(18)/U.FitResult.par(18)),...
                                'fitPar',fitPar,...
                                'DeltaChi2',DeltaChi2);
                        else
                            eta(i) = U.FitResult.err(18)*1e10;
                        end
                        Chi2(i)      = U.FitResult.chi2min;
                        mnuSq(i)     = U.ModelObj.mnuSq_i+U.FitResult.par(1);
                        mnuSq_err(i) = U.FitResult.err(1);
                        E0(i)        = U.ModelObj.Q_i+U.FitResult.par(2);
                        E0_err(i)    = U.FitResult.err(2);
                        Bkg(i)       = U.ModelObj.BKG_RateSec_i+U.FitResult.par(3);
                        Bkg_err(i)   = U.FitResult.err(3);
                        Norm(i)      = U.FitResult.par(3+U.ModelObj.nPixels:3+2*U.ModelObj.nPixels-1) + 1;
                        Norm_err(i)  = U.FitResult.err(3+U.ModelObj.nPixels:3+2*U.ModelObj.nPixels-1);
                        sprintf('Final Result: eta = %g',eta(i))
                   end
               end
               save(savename,'mNu','eta','etaFit','Chi2','mnuSq','mnuSq_err','E0','E0_err','Bkg','Bkg_err','Norm','Norm_err');
           else
               load(savename,'eta');
               obj.etaSensitivity = eta;
               plotchi2scan(savename);
           end
       end
   end
   methods
       function FakeSystBreakdown(obj,varargin)
           p=inputParser;
           p.addParameter('TwinBias_mnuSq',0,@(x)isfloat(x));
           p.addParameter('range',40,@(x)isfloat(x));
           p.addParameter('Init_Opt','',@(x)iscell(x) || isempty(x));
           p.addParameter('Plot','ON',@(x)ismember(x,{'ON','OFF'}));
           p.addParameter('Recompute','OFF',@(x)ismember(x,{'ON','OFF'}));
           p.parse(varargin{:});
           TwinBias_mnuSq = p.Results.TwinBias_mnuSq;
           range          = p.Results.range;
           Init_Opt       = p.Results.Init_Opt;
           Plot           = p.Results.Plot;
           Recompute      = p.Results.Recompute;
           
           matFilePath = [getenv('SamakPath'),sprintf('RelicNuBkg/Misc/')];
           savename = [matFilePath,sprintf('FakeSensitivityBreakdown_%s_mnuSq%g_range%g',obj.Params,TwinBias_mnuSq,range)];
           savename = [savename,'.mat'];
           
           if strcmp(obj.Params,'TDR')
                initfile=@ref_RelicNuBkg_DesignReport;
                TwinBias_Q=18575;
                RingList=1:14;
                SysBudget=67;
                NonPoissonScaleFactor=1;
            elseif strcmp(obj.Params,'Formaggio')
                initfile=@ref_RelicNuBkg_Formaggio;
                TwinBias_Q=18575;
                RingList=1:14;
            elseif strcmp(obj.Params,'KNM1')
                initfile=@ref_RelicNuBkg_KNM1;
                TwinBias_Q=18573.73;
                RingList=1:12;
                SysBudget=24;
                NonPoissonScaleFactor=1.064;
            elseif strcmp(obj.Params,'KNM2_Prompt')
                initfile=@ref_RelicNuBkg_KNM2;
                TwinBias_Q=18573.7;
                RingList=1:12;
                NonPoissonScaleFactor=1.1120;
            end

           if exist(savename,'file') && strcmp(Recompute,'OFF')
               load(savename,'X','Y');
           else
           
               A = RunAnalysis('RunNr',10,...             
                    'FakeInitFile',initfile,...
                    'Init_Opt',Init_Opt,...
                    'RecomputeFakeRun','ON',...
                    'chi2','chi2CMShape',...                 % uncertainties: statistical or stat + systematic uncertainties
                    'NonPoissonScaleFactor',NonPoissonScaleFactor,...
                    'DataType','Fake',...                 % can be 'Real' or 'Twin' -> Monte Carlo
                    'RingList',RingList,...
                    'TwinBias_Q',TwinBias_Q,...
                    'fixPar','mNu E0 Norm Bkg eta',...                   % free Parameter!!
                    'minuitOpt','min ; minos',...         % technical fitting options (minuit)
                    'FSDFlag','SibilleFull',...          % final state distribution                        !!check ob initfile hier überschrieben wird
                    'ELossFlag','KatrinT2',...            % energy loss function
                    'SysBudget',SysBudget,...                    % defines syst. uncertainties -> in GetSysErr.m;
                    'DopplerEffectFlag','FSD',...
                    'SynchrotronFlag','OFF',...
                    'AngularTFFlag','OFF');
                    
               A.exclDataStart = A.GetexclDataStart(range);
               A.ModelObj.mnuSq_i = A.TwinBias_mnuSq;
               obj.M = A;

               A.Fit;
               ErrTotal = A.FitResult.err(18)*1e10;
               A.ComputeCM('SysEffects',struct(),'BkgCM','OFF','BkgPtCM','OFF');
               A.Fit;
               ErrNP = A.FitResult.err(18)*1e10;
               A.NonPoissonScaleFactor = 1;
               A.ComputeCM('SysEffects',struct(),'BkgCM','OFF','BkgPtCM','OFF');
               A.Fit;
               ErrStat = A.FitResult.err(18)*1e10;
               A.ComputeCM('SysEffects',struct('FSD','ON'),'BkgCM','OFF','BkgPtCM','OFF');
               A.Fit;
               ErrFSD = A.FitResult.err(18)*1e10;
               A.ComputeCM('SysEffects',struct('RF','ON'),'BkgCM','OFF','BkgPtCM','OFF');
               A.Fit;
               ErrRF = A.FitResult.err(18)*1e10;
               A.ComputeCM('SysEffects',struct('TASR','ON'),'BkgCM','OFF','BkgPtCM','OFF');
               A.Fit;
               ErrTASR = A.FitResult.err(18)*1e10;
               A.ComputeCM('SysEffects',struct('Stack','ON'),'BkgCM','OFF','BkgPtCM','OFF');
               A.Fit;
               ErrStack = A.FitResult.err(18)*1e10;
               A.ComputeCM('SysEffects',struct('FPDeff','ON'),'BkgCM','OFF','BkgPtCM','OFF');
               A.Fit;
               ErrFPD = A.FitResult.err(18)*1e10;
               A.ComputeCM('SysEffects',struct('TC','ON'),'BkgCM','OFF','BkgPtCM','OFF');
               A.Fit;
               ErrTC = A.FitResult.err(18)*1e10;
               A.ComputeCM('SysEffects',struct(),'BkgCM','ON','BkgPtCM','OFF');
               A.Fit;
               ErrBkg = A.FitResult.err(18)*1e10;
               ErrFSD = sqrt(ErrFSD.^2-ErrStat.^2);
               ErrRF = sqrt(ErrRF.^2-ErrStat.^2);
               ErrTASR = sqrt(ErrTASR.^2-ErrStat.^2);
               ErrStack = sqrt(ErrStack.^2-ErrStat.^2);
               ErrFPD = sqrt(ErrFPD.^2-ErrStat.^2);
               ErrTC = sqrt(ErrTC.^2-ErrStat.^2);
               ErrBkg = sqrt(ErrBkg.^2-ErrStat.^2);
               ErrNP = sqrt(ErrNP.^2-ErrStat.^2);

               X = categorical({'Total','Stat','FSD','RF','TASR','Stack','FPD','TC','Bkg','NP'});
               X = reordercats(X,{'Total','Stat','FSD','RF','TASR','Stack','FPD','TC','Bkg','NP'});
               Y = [ErrTotal ErrStat ErrFSD ErrRF ErrTASR ErrStack ErrFPD ErrTC ErrBkg ErrNP];
               save(savename,'X','Y');
           end
           if strcmp(Plot,'ON')
               bar(X,Y);
           end
           PrettyFigureFormat;
       end
       function SystBreakdown(obj,varargin)
           p=inputParser;
           p.addParameter('TwinBias_mnuSq',0,@(x)isfloat(x));
           p.addParameter('range',40,@(x)isfloat(x));
           p.addParameter('DataType','Twin',@(x)ismember(x,{'Twin','Real'}));
           p.addParameter('Plot','ON',@(x)ismember(x,{'ON','OFF'}));
           p.addParameter('Recompute','OFF',@(x)ismember(x,{'ON','OFF'}));
           p.addParameter('CorrectErr','ON',@(x)ismember(x,{'ON','OFF'}));
           p.parse(varargin{:});
           TwinBias_mnuSq = p.Results.TwinBias_mnuSq;
           range          = p.Results.range;
           DataType       = p.Results.DataType;
           Plot           = p.Results.Plot;
           Recompute      = p.Results.Recompute;
           CorrectErr     = p.Results.CorrectErr;
           
           matFilePath = [getenv('SamakPath'),sprintf('RelicNuBkg/Misc/')];
           savename = [matFilePath,sprintf('SensitivityBreakdown_%s_mnuSq%g_range%g',obj.Params,TwinBias_mnuSq,range)];
           if strcmp(DataType,'Real')
               savename = [savename,'_RealData'];
           end
           savename = [savename,'.mat'];
           
           fitter = 'minuit';
           if strcmp(obj.Params,'KNM1')
               SysBudget=24;
               NP       = 1.064;
           elseif strcmp(obj.Params,'KNM2_Prompt')
               SysBudget=40;
               NP       =1.112;
           elseif strcmp(obj.Params,'TDR')
               SysBudget=67;
               NP       =1;
           end

           if exist(savename,'file') && strcmp(Recompute,'OFF')
               load(savename,'X','Y');
               if strcmp(Plot,'ON')
                   bar(X,Y);
               end
           else
           
               if strcmp(obj.Params,'KNM1')
                   A = MultiRunAnalysis('RunList',obj.Params,... % runlist defines which runs are analysed -> set MultiRunAnalysis.m -> function: GetRunList()
                            'chi2','chi2CMShape',...                 % uncertainties: statistical or stat + systematic uncertainties
                            'DataType',DataType,...                 % can be 'Real' or 'Twin' -> Monte Carlo
                            'fixPar','mNu E0 Norm Bkg eta',...                   % free Parameter!!
                            'RadiativeFlag','ON',...              % theoretical radiative corrections applied in model
                            'NonPoissonScaleFactor',NP,...     % background uncertainty are enhanced
                            'fitter',fitter,...
                            'minuitOpt','min ; minos',...         % technical fitting options (minuit)
                            'FSDFlag','SibilleFull',...          % final state distribution
                            'ELossFlag','KatrinT2',...            % energy loss function
                            'SysBudget',SysBudget,...                    % defines syst. uncertainties -> in GetSysErr.m;
                            'DopplerEffectFlag','FSD',...
                            'Twin_SameCDFlag','OFF',...
                            'Twin_SameIsotopFlag','OFF',...
                            'SynchrotronFlag','ON',...
                            'AngularTFFlag','OFF',...
                            'TwinBias_Q',18573.73,...
                            'TwinBias_mnuSq',TwinBias_mnuSq);
               elseif strcmp(obj.Params,'KNM2_Prompt')
                   A = MultiRunAnalysis('RunList',obj.Params,... % runlist defines which runs are analysed -> set MultiRunAnalysis.m -> function: GetRunList()
                            'chi2','chi2CMShape',...                 % uncertainties: statistical or stat + systematic uncertainties
                            'DataType',DataType,...                 % can be 'Real' or 'Twin' -> Monte Carlo
                            'fixPar','mNu E0 Norm Bkg eta',...                   % free Parameter!!
                            'RadiativeFlag','ON',...              % theoretical radiative corrections applied in model
                            'NonPoissonScaleFactor',NP,...     % background uncertainty are enhanced
                            'fitter',fitter,...
                            'minuitOpt','min ; minos',...         % technical fitting options (minuit)
                            'FSDFlag','KNM2',...          % final state distribution
                            'ELossFlag','KatrinT2A20',...            % energy loss function
                            'SysBudget',SysBudget,...                    % defines syst. uncertainties -> in GetSysErr.m;
                            'DopplerEffectFlag','FSD',...
                            'Twin_SameCDFlag','OFF',...
                            'Twin_SameIsotopFlag','OFF',...
                            'SynchrotronFlag','ON',...
                            'AngularTFFlag','ON',...
                            'TwinBias_Q',18573.7,...
                            'TwinBias_mnuSq',TwinBias_mnuSq,...
                            'FSD_Sigma',sqrt(0.0124+0.0025),...
                            'TwinBias_FSDSigma',sqrt(0.0124+0.0025),...
                            'BKG_PtSlope',3*1e-06,...
                            'TwinBias_BKG_PtSlope',3*1e-06);
               end
                    
               A.exclDataStart = A.GetexclDataStart(range);
               A.ModelObj.mnuSq_i = A.TwinBias_mnuSq;
               obj.M = A;

               A.Fit;
               if strcmp(CorrectErr,'ON')
                   ErrTotal = obj.CorrectErr('Parameter','eta','value',A.FitResult.par(18)*1e10,'eta',A.FitResult.par(18)*1e10,'minchi2',A.FitResult.chi2min,'factor',(1+A.FitResult.err(18)/A.FitResult.par(18)));
               else
                   ErrTotal = 1e10*(abs(A.FitResult.errPos(5))+abs(A.FitResult.errNeg(5)))/2;
               end
               A.ComputeCM('SysEffects',struct(),'BkgCM','OFF','BkgPtCM','OFF');
               A.Fit;
               if strcmp(CorrectErr,'ON')
                    ErrNP = obj.CorrectErr('Parameter','eta','value',A.FitResult.par(18)*1e10,'eta',A.FitResult.par(18)*1e10,'minchi2',A.FitResult.chi2min,'factor',1.5,'SystSelect','NP');
               else
                   ErrNP = 1e10*(abs(A.FitResult.errPos(5))+abs(A.FitResult.errNeg(5)))/2;
               end
               A.NonPoissonScaleFactor = 1;
               A.ComputeCM('SysEffects',struct(),'BkgCM','OFF','BkgPtCM','OFF');
               A.Fit;
               if strcmp(CorrectErr,'ON')
                    ErrStat = obj.CorrectErr('Parameter','eta','value',A.FitResult.par(18)*1e10,'eta',A.FitResult.par(18)*1e10,'minchi2',A.FitResult.chi2min,'factor',1.5,'SystSelect','None');
               else
                   ErrStat = 1e10*(abs(A.FitResult.errPos(5))+abs(A.FitResult.errNeg(5)))/2;
               end
               A.ComputeCM('SysEffects',struct('FSD','ON'),'BkgCM','OFF','BkgPtCM','OFF');
               A.Fit;
               if strcmp(CorrectErr,'ON')
                    ErrFSD = obj.CorrectErr('Parameter','eta','value',A.FitResult.par(18)*1e10,'eta',A.FitResult.par(18)*1e10,'minchi2',A.FitResult.chi2min,'factor',(1+A.FitResult.err(18)/A.FitResult.par(18)),'SystSelect','FSD');
               else 
                   ErrFSD = 1e10*(abs(A.FitResult.errPos(5))+abs(A.FitResult.errNeg(5)))/2;
               end
               A.ComputeCM('SysEffects',struct('RF','ON'),'BkgCM','OFF','BkgPtCM','OFF');
               A.Fit;
               if strcmp(CorrectErr,'ON')
                    ErrRF = obj.CorrectErr('Parameter','eta','value',A.FitResult.par(18)*1e10,'eta',A.FitResult.par(18)*1e10,'minchi2',A.FitResult.chi2min,'factor',(1+A.FitResult.err(18)/A.FitResult.par(18)),'SystSelect','RF');
               else 
                   ErrRF = 1e10*(abs(A.FitResult.errPos(5))+abs(A.FitResult.errNeg(5)))/2;
               end
               A.ComputeCM('SysEffects',struct('TASR','ON'),'BkgCM','OFF','BkgPtCM','OFF');
               A.Fit;
               if strcmp(CorrectErr,'ON')
                    ErrTASR = obj.CorrectErr('Parameter','eta','value',A.FitResult.par(18)*1e10,'eta',A.FitResult.par(18)*1e10,'minchi2',A.FitResult.chi2min,'factor',(1+A.FitResult.err(18)/A.FitResult.par(18)),'SystSelect','TASR');
               else
                   ErrTASR = 1e10*(abs(A.FitResult.errPos(5))+abs(A.FitResult.errNeg(5)))/2;
               end
               A.ComputeCM('SysEffects',struct('Stack','ON'),'BkgCM','OFF','BkgPtCM','OFF');
               A.Fit;
               if strcmp(CorrectErr,'ON')
                    ErrStack = obj.CorrectErr('Parameter','eta','value',A.FitResult.par(18)*1e10,'eta',A.FitResult.par(18)*1e10,'minchi2',A.FitResult.chi2min,'factor',(1+A.FitResult.err(18)/A.FitResult.par(18)),'SystSelect','Stack');
               else
                   ErrStack = 1e10*(abs(A.FitResult.errPos(5))+abs(A.FitResult.errNeg(5)))/2;
               end
               A.ComputeCM('SysEffects',struct('FPDeff','ON'),'BkgCM','OFF','BkgPtCM','OFF');
               A.Fit;
               if strcmp(CorrectErr,'ON')
                    ErrFPD = obj.CorrectErr('Parameter','eta','value',A.FitResult.par(18)*1e10,'eta',A.FitResult.par(18)*1e10,'minchi2',A.FitResult.chi2min,'factor',1.5,'SystSelect','FPDeff');
               else 
                   ErrFPD = 1e10*(abs(A.FitResult.errPos(5))+abs(A.FitResult.errNeg(5)))/2;
               end
               A.ComputeCM('SysEffects',struct('TC','ON'),'BkgCM','OFF','BkgPtCM','OFF');
               A.Fit;
               if strcmp(CorrectErr,'ON')
                    ErrTC = obj.CorrectErr('Parameter','eta','value',A.FitResult.par(18)*1e10,'eta',A.FitResult.par(18)*1e10,'minchi2',A.FitResult.chi2min,'factor',1.5,'SystSelect','TC');
               else 
                   ErrTC = 1e10*(abs(A.FitResult.errPos(5))+abs(A.FitResult.errNeg(5)))/2;
               end
               if strcmp(obj.Params,'KNM2_Prompt')
                   A.ComputeCM('SysEffects',struct('LongPlasma','ON','BkgCM','OFF','BkgPtCM','OFF'));
                   A.Fit;
                   if strcmp(CorrectErr,'ON')
                        ErrPlasma = obj.CorrectErr('Parameter','eta','value',A.FitResult.par(18)*1e10,'eta',A.FitResult.par(18)*1e10,'minchi2',A.FitResult.chi2min,'factor',1.5,'SystSelect','LongPlasma');
                   else 
                        ErrPlasma = 1e10*(abs(A.FitResult.errPos(5))+abs(A.FitResult.errNeg(5)))/2;
                   end 
                   A.ComputeCM('SysEffects',struct(),'BkgCM','OFF','BkgPtCM','ON');
                   A.Fit;
                   if strcmp(CorrectErr,'ON')
                        ErrPT = obj.CorrectErr('Parameter','eta','value',A.FitResult.par(18)*1e10,'eta',A.FitResult.par(18)*1e10,'minchi2',A.FitResult.chi2min,'factor',1.5,'SystSelect','BkgPtCM');
                   else 
                        ErrPT = 1e10*(abs(A.FitResult.errPos(5))+abs(A.FitResult.errNeg(5)))/2;
                   end 
               else
                   ErrPlasma=0;ErrPT=0;
               end
               A.ComputeCM('SysEffects',struct(),'BkgCM','ON','BkgPtCM','ON');
               A.Fit;
               if strcmp(CorrectErr,'ON')
                    ErrBkg = obj.CorrectErr('Parameter','eta','value',A.FitResult.par(18)*1e10,'eta',A.FitResult.par(18)*1e10,'minchi2',A.FitResult.chi2min,'factor',1.5,'SystSelect','BkgCM');
               else
                   ErrBkg = 1e10*(abs(A.FitResult.errPos(5))+abs(A.FitResult.errNeg(5)))/2;
               end
               ErrFSD = sqrt(ErrFSD.^2-ErrStat.^2);
               ErrRF = sqrt(ErrRF.^2-ErrStat.^2);
               ErrTASR = sqrt(ErrTASR.^2-ErrStat.^2);
               ErrStack = sqrt(ErrStack.^2-ErrStat.^2);
               ErrFPD = sqrt(ErrFPD.^2-ErrStat.^2);
               ErrTC = sqrt(ErrTC.^2-ErrStat.^2);
               ErrPlasma = sqrt(ErrPlasma.^2-ErrStat.^2);
               ErrPT = sqrt(ErrPT.^2-ErrStat.^2);
               ErrBkg = sqrt(ErrBkg.^2-ErrStat.^2);
               ErrNP = sqrt(ErrNP.^2-ErrStat.^2);

               X = categorical({'Total','Stat','FSD','RF','TASR','Stack','FPD','TC','Bkg','NP','Plasma','PT'});
               X = reordercats(X,{'Total','Stat','FSD','RF','TASR','Stack','FPD','TC','Bkg','NP','Plasma','PT'});
               Y = [ErrTotal ErrStat ErrFSD ErrRF ErrTASR ErrStack ErrFPD ErrTC ErrBkg ErrNP ErrPlasma ErrPT];
               save(savename,'X','Y');
           end
           if strcmp(Plot,'ON')
               bar(X,Y);
           end
           PrettyFigureFormat;
       end
       function Error = CorrectErr(obj,varargin)
           p = inputParser;
           p.addParameter('Parameter','',@(x)ismember(x,{'mNu','E0','Norm','Bkg','eta'}));
           p.addParameter('value',0,@(x)isfloat(x));
           p.addParameter('eta',0,@(x)isfloat(x));
           p.addParameter('minchi2',0,@(x)isfloat(x));
           p.addParameter('factor',0,@(x)isfloat(x));
           p.addParameter('SystSelect','OFF',@(x)ismember(x,{'RF','TASR','Stack','FSD','TC','FPDeff','LongPlasma','BkgPtCM','BkgCM','NP','None','OFF'}));
           p.addParameter('DeltaChi2',1,@(x)isfloat(x));
           p.addParameter('fitPar','mNu E0 Norm Bkg',@(x)ischar(x));
           p.parse(varargin{:});
           Parameter  = p.Results.Parameter;
           value      = p.Results.value;
           eta        = p.Results.eta;
           minchi2    = p.Results.minchi2;
           factor     = p.Results.factor;
           SystSelect = p.Results.SystSelect;
           DeltaChi2  = p.Results.DeltaChi2;
           fitPar     = p.Results.fitPar;
           
           Chi2upper  = minchi2;
           Chi2lower  = minchi2;
           valueupper = value;
           valuelower = value;
           value0     = value;
           obj.SetParam('Parameter',Parameter,'INIT',1,'SystSelect',SystSelect,'fitPar',fitPar)
           obj.M.ModelObj.eta_i = eta;
           obj.M.ModelObj.eta   = eta;
           if strcmp(Parameter,'mNu')
               factor = 2;
           elseif strcmp(Parameter,'E0')
               factor = 1+5e-6;
           elseif strcmp(Parameter,'Norm')
               factor = 1.005;
           elseif strcmp(Parameter,'Bkg')
               factor = 1.004;
           end
           
           if strcmp(Parameter,'eta')
               while Chi2upper-minchi2 < DeltaChi2 && abs(Chi2upper-minchi2-DeltaChi2) > 0.01
                   Chi2lower  = Chi2upper;
                   valuelower = valueupper;
                   obj.SetParam('Parameter',Parameter,'value',abs(factor*value));
                   obj.M.ModelObj.ComputeNormFactorTBDDS;
                   obj.M.ModelObj.ComputeTBDDS;
                   obj.M.ModelObj.ComputeTBDIS;
                   obj.M.Fit;
                   valueupper = abs(factor*value);
                   if (obj.M.FitResult.chi2min-Chi2upper)<0.5*DeltaChi2
                       factor=2*factor;
                   end
                   Chi2upper  = obj.M.FitResult.chi2min;
                   value = valueupper;
               end
           else
               while Chi2upper-minchi2 < DeltaChi2
                   Chi2lower  = Chi2upper;
                   valuelower = valueupper;
                   obj.SetParam('Parameter',Parameter,'value',abs(factor*value));
                   obj.M.ModelObj.ComputeNormFactorTBDDS;
                   obj.M.ModelObj.ComputeTBDDS;
                   obj.M.ModelObj.ComputeTBDIS;
                   obj.M.Fit;
                   Chi2upper  = obj.M.FitResult.chi2min;
                   valueupper = abs(factor*value);
                   value = valueupper;
               end
           end
           while abs(Chi2upper-minchi2-DeltaChi2)>0.01 && abs(Chi2lower-minchi2-DeltaChi2)>0.01
               obj.SetParam('Parameter',Parameter,'value',((valueupper-valuelower)./(Chi2upper-Chi2lower)).*(DeltaChi2+minchi2)-Chi2upper.*((valueupper-valuelower)./(Chi2upper-Chi2lower))+valueupper);
               value = ((valueupper-valuelower)./(Chi2upper-Chi2lower)).*(DeltaChi2+minchi2)-Chi2upper.*((valueupper-valuelower)./(Chi2upper-Chi2lower))+valueupper;
               obj.M.ModelObj.ComputeNormFactorTBDDS;
               obj.M.ModelObj.ComputeTBDDS;
               obj.M.ModelObj.ComputeTBDIS;
               obj.M.Fit;
               if obj.M.FitResult.chi2min-minchi2>DeltaChi2
                   valueupper = ((valueupper-valuelower)./(Chi2upper-Chi2lower)).*(DeltaChi2+minchi2)-Chi2upper.*((valueupper-valuelower)./(Chi2upper-Chi2lower))+valueupper;
                   Chi2upper  = obj.M.FitResult.chi2min;
               else
                   valuelower = ((valueupper-valuelower)./(Chi2upper-Chi2lower)).*(DeltaChi2+minchi2)-Chi2upper.*((valueupper-valuelower)./(Chi2upper-Chi2lower))+valueupper;
                   Chi2lower  = obj.M.FitResult.chi2min;
               end
           end
           
           Error = value-value0;
           
       end
       function SystBias(obj,varargin)
            p=inputParser;
            p.addParameter('range',40,@(x)isfloat(x));
            p.addParameter('fitPar','mNu E0 Norm Bkg',@(x)ischar(x));
            p.addParameter('Syst','ON',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('SystBudget',24,@(x)isfloat(x));
            p.addParameter('TwinBias_mnuSq',1,@(x)isfloat(x));
            p.addParameter('pullFlag',3);
            p.addParameter('Recompute','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('Plot','ON',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('DeltaChi2',2.71,@(x)isfloat(x));
            p.parse(varargin{:});

            range          = p.Results.range;
            fitPar         = p.Results.fitPar;
            Syst           = p.Results.Syst;
            SystBudget     = p.Results.SystBudget;
            TwinBias_mnuSq = p.Results.TwinBias_mnuSq;
            pullFlag       = p.Results.pullFlag;
            Recompute      = p.Results.Recompute;
            DeltaChi2      = p.Results.DeltaChi2;
            
            SystEffects = ["RF","TASR","Stack","FSD","TC","FPDeff","BkgCM"];
            Bias = zeros(1,numel(SystEffects));
            matFilePath = [getenv('SamakPath'),sprintf('RelicNuBkg/Misc/')];
            savename=[matFilePath,sprintf('SystEffectOnSensitivity_%s_mnuSq%g.mat',obj.Params,TwinBias_mnuSq)];
            
            if exist(savename,'file') && strcmp(Recompute,'OFF')
                load(savename);
            else
            
               obj.M = MultiRunAnalysis('RunList',obj.Params,... % runlist defines which runs are analysed -> set MultiRunAnalysis.m -> function: GetRunList()
                'chi2','chi2CMShape',...                 % uncertainties: statistical or stat + systematic uncertainties
                'DataType','Twin',...                 % can be 'Real' or 'Twin' -> Monte Carlo
                'fixPar',fitPar,...                   % free Parameter!!
                'RadiativeFlag','ON',...              % theoretical radiative corrections applied in model
                'NonPoissonScaleFactor',1.064,...     % background uncertainty are enhanced
                'fitter','minuit',...
                'minuitOpt','min ; minos',...         % technical fitting options (minuit)
                'pullFlag',pullFlag,...
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
            
            TBDISBias = zeros(numel(obj.M.RunData.qU),numel(SystEffects));
            
            for i=1:numel(SystEffects)
                if ~strcmp(SystEffects(i),'Bkg')
                    obj.M.ComputeCM('SysEffects',struct(SystEffects(i),'ON'),'BkgCM','OFF');
                    TBDISBias(:,i) = sqrt(diag(obj.M.FitCM));
                else
                    obj.M.ComputeCM('BkgCm','ON');
                    TBDISBias(:,i) = sqrt(diag(obj.M.FitCM));
                end
            end
            
            obj.Chi2Scan_Twin('Recompute','OFF',...
                    'range',range,...
                    'RunList',obj.Params,...
                    'fitPar',fitPar,...
                    'Syst','OFF',...
                    'TwinBias_mnuSq',TwinBias_mnuSq,...
                    'pullFlag',pullFlag,...
                    'etalower',0,...
                    'etaupper',2e11,...
                    'mode','SEARCH',...
                    'Plot','OFF',...
                    'DeltaChi2',DeltaChi2);
                
                   etaStat = obj.etaSensitivity;

                for i=1:numel(SystEffects)

                   obj.Chi2Scan_Twin('Recompute',Recompute,...
                       'TBDISBias',TBDISBias(:,i),...
                       'range',range,...
                       'RunList',obj.Params,...
                       'fitPar',fitPar,...
                       'Syst','OFF',...
                       'TwinBias_mnuSq',TwinBias_mnuSq,...
                       'pullFlag',pullFlag,...
                       'etalower',0,...
                       'etaupper',2.001e11,...
                       'mode','SEARCH',...
                       'Plot','OFF',...
                       'DeltaChi2',DeltaChi2);

                   Bias(i) = obj.etaSensitivity - etaStat;
                end
                save(savename,'SystEffects','Bias');
            end
            X = categorical({'RF','TASR','Stack','FSD','TC','FPDeff','BkgCM'});
            X = reordercats(X,{'RF','TASR','Stack','FSD','TC','FPDeff','BkgCM'});
            bar(X,Bias);
            PrettyFigureFormat;
       end
   end
    
end