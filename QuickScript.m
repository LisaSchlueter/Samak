%CheckmnuSqDist('Nfit',1000);

% Nfit=1000;
% statfluct   = zeros(numel(M.RunData.qU),Nfit);
% statfluct_1=mvnrnd(M.RunData.TBDIS',M.FitCM,Nfit)';
% for j=1:Nfit
%     for i=1:numel(M.RunData.qU)
%         gm=gmdistribution(M.RunData.TBDIS(i),(M.NonPoissonScaleFactor)*M.RunData.TBDIS(i));
%         statfluct(i,j) = random(gm)-M.RunData.TBDIS(i);
%     end
%     statfluct_1(:,j)=statfluct_1(:,j)-M.RunData.TBDIS;
% end

% A=ref_RelicNuBkg_KNM1('TTFSD','OFF','HTFSD','OFF','DTFSD','OFF','mnuSq_i',2);
% B=ref_RelicNuBkg_KNM1('mnuSq_i',2);
% a=semilogy(A.Te-A.Q,A.TBDDS,'LineWidth',2);
% hold on;
% b=semilogy(A.Te-A.Q,A.TBDDS_R,'LineWidth',2);
% c=semilogy(B.Te-B.Q,B.TBDDS,'LineWidth',2);
% d=semilogy(B.Te-B.Q,B.TBDDS_R,'LineWidth',2);
% ylabel('Rate per Energy Bin');
% xlabel('Energy - E_{0} (eV)');
% legend([a b c d],'\beta decay, FSD off','C\nuB, FSD off','\beta decay, FSD on','C\nuB, FSD on');
% legend boxoff;
% xlim([-10 4]);
% ylim([1e-10 1e14]);
% PrettyFigureFormat;

%  load('./RelicNuBkg/Results_mNuFree.mat');
%  %eta_corrected = eta;
%  eta_Chi2_corrected = eta_Chi2;
%  load('./Results_Local_OnlyChi2.mat');
%  eta_local = eta_Chi2;
%  load('./RelicNuBkg/Results_mNuFree_uncorrected.mat');
%  %plot(linspace(0,2,21),eta_fit,'LineWidth',2);
%  hold on;
%  plot(linspace(0,2,21),eta_Chi2,'LineWidth',2);
%  %plot(linspace(0,2,21),eta_corrected,'LineWidth',2);
%  a=plot(linspace(0,2,21),eta_Chi2_corrected,'LineWidth',2);
%  b=plot(linspace(0,2,21),eta_local,'LineWidth',2);
%  ylim([1.39e11 1.95e11]);
%  xlabel('m_{\nu}^{2} (eV^2)');
%  ylabel('1\sigma upper Limit of \eta');
%  %lgd=legend([a b],{'Server','Local'},'Location','northeast','box','off');
%  %lgd.FontSize = 12;
%  PrettyFigureFormat;
%  hold off;

%% INFO:
    % - rL=R1 line (v=0): plot(R2,m1.*(r2-R1).^2./R1.^2)
    % - correct order of surf plot for rL, v: surf(R2,m2,v)
    % - rotating potential contour plot: contour(x,y,U',50)
close all;
%clear all;
%startangle=linspace(140,190,10);
%launchangle=linspace(-25,25,10);
%v=TransitVel('rho2',4000,'R2',8660000,'startangle',startangle,'launchangle',launchangle);
vel=TransitVel('rho2',4000,'R2',8660000,'startangle',157,'launchangle',-14,'Startheight',0);
%TransitVel('rho2',5000,'R2',30927000);
%TransitVel('startangle',168,'launchangle',-14);
% N=10;
% rho2=linspace(3000,11342,N);
% R2  =linspace(0.5*3092700,10*3092700,N);
% v=zeros(N,N);
% for i=1:N
%     for j=1:N
%         v(i,j)=TransitVel('rho2',rho2(i),'R2',R2(j),'plot','OFF');
%     end
% end

r_outer = 0.8; %meter
r_inner = 0.795;
mass = 4/3*pi*(r_outer^3-r_inner^3)*8920
energy = 0.5*mass*vel^2

function vel=TransitVel(varargin)
    p=inputParser;
    p.addParameter('rho2',11342,@(x)isfloat(x));
    p.addParameter('R2',3092700,@(x)isfloat(x));
    p.addParameter('launchangle',0,@(x)isfloat(x));
    p.addParameter('startangle',180,@(x)isfloat(x));
    p.addParameter('Plot','ON',@(x)ismember(x,{'ON','OFF'}));
    p.addParameter('Startheight',0,@(x)isfloat(x));
    p.parse(varargin{:});
    %% Settings
    rho2        = p.Results.rho2;
    R2          = p.Results.R2; %eV²
    launchangle = p.Results.launchangle;
    startangle  = p.Results.startangle;
    Plot        = p.Results.Plot;
    Startheight = p.Results.Startheight;
    
    G=6.67428e-11;
    rho1=11342; %kg/m³
    R1=3092700; %m
    r2=R1+R2+1612e3;
    m1=4*pi./3.*R1^3.*rho1; %kg
    m2= 4.*pi./3.*R2.^3.*rho2;%linspace(0,100*7.117e24,N);
    T=2.*pi.*sqrt(r2.^3./(G.*(m1+m2)));
    r_B=r2.*m2./(m1+m2); %barycentre
    rL=r2./(1+sqrt(m2/m1));
    v=sqrt(-2*G*(m1./rL+m2./(r2-rL+1e-100)-m1./R1-m2./(r2-R1)));
    if m2>m1.*(r2-R1).^2/R1.^2
        v=0;
    end
    vnull=sqrt(-2*G*(m1./rL-m1./R1));
    
    Ngrid=1000;
    x=linspace(-2*R2,2*R2,Ngrid);
    y=linspace(-r2./2-r_B,1.5.*r2-r_B,Ngrid);
    v_r=zeros(Ngrid,Ngrid);
    U=zeros(Ngrid,Ngrid);
    U_L=zeros(Ngrid,Ngrid);
    F=zeros(Ngrid,Ngrid,2);
    vx=[0 0 0];
    for i=1:numel(x)
        for j=1:numel(y)
            v_r(i,j)=(2.*pi.*sqrt(x(i)^2+(y(j))^2))./T;
            U(i,j)=-G.*m1./(sqrt(x(i).^2+(y(j)+r_B).^2))-G.*m2./(sqrt(x(i).^2+(y(j)-(r2-r_B)).^2))-...
                v_r(i,j).^2./2;
            U_L(i,j)=U(i,j);
            F(i,j,:)=G.*m1./(x(i).^2+(y(j)+r_B).^2).*[-x(i) -y(j)-r_B]./sqrt(x(i)^2+(-y(j)-r_B)^2)+...
                G.*m2./(x(i).^2+(y(j)-(r2-r_B)).^2).*[-x(i) -y(j)+(r2-r_B)]./sqrt(x(i)^2+(y(j)-(r2-r_B))^2)+...
                v_r(i,j).^2./sqrt(x(i)^2+(y(j))^2).*[x(i) y(j)]./sqrt(x(i)^2+(y(j))^2)+...
                2.*[2*pi./T.*vx(2) -vx(1).*2*pi./T];                                  %F(x,y)=[F(x,y,1),F(x,y,2)]
            if x(i)^2+(y(j)+r_B)^2<R1^2 || (y(j)-(r2-r_B))^2+x(i)^2<R2^2
                U(i,j)=0;
                F(i,j,:)=0;
            end
        end
    end
    if strcmp(Plot,'ON')
        figure(1);
        hold on;
        Nplot=20;
        for i=1:Nplot
            for j=1:Nplot
                plot([x(Ngrid/Nplot*i-(Ngrid/Nplot-1)) x(Ngrid/Nplot*i-(Ngrid/Nplot-1))+1e5*F(Ngrid/Nplot*i-Ngrid/Nplot+1,Ngrid/Nplot*j-Ngrid/Nplot+1,1)],...
                    [y(Ngrid/Nplot*j-Ngrid/Nplot+1) y(Ngrid/Nplot*j-Ngrid/Nplot+1)+1e5*F(Ngrid/Nplot*i-Ngrid/Nplot+1,Ngrid/Nplot*j-Ngrid/Nplot+1,2)]);
            end
        end
        contour(x,y,U');
    end
    sprintf('Transit velocity (static): %g m/s',v)
    
    Force=@(x,y,vx) G.*m1./(x.^2+(y+r_B).^2).*[-x -y-r_B]./sqrt(x^2+(-y-r_B)^2)+...
                G.*m2./(x.^2+(y-(r2-r_B)).^2).*[-x -y+(r2-r_B)]./sqrt(x^2+(y-(r2-r_B))^2)+...
                ((2.*pi.*sqrt(x^2+(y)^2))./T).^2./sqrt(x^2+(y)^2).*[x y]./sqrt(x^2+(y)^2)+...
                2.*[2*pi./T.*vx(2) -vx(1).*2*pi./T];
    Potential=@(x,y) -G.*m1./(sqrt(x.^2+(y+r_B).^2))-G.*m2./(sqrt(x.^2+(y-(r2-r_B)).^2))-...
                ((2.*pi.*sqrt(x^2+(y)^2))./T).^2./2;
    PotentialDiffL1=zeros(1,Ngrid);
    PotentialDiffL2=zeros(1,Ngrid);
    L1range=linspace(-r_B+R1,-r_B+R1+1612e3,Ngrid+1);
    L2range=linspace(y(1),-r_B-R1,Ngrid+1);
    for i=2:Ngrid+1
       PotentialDiffL1(i-1)=abs(Potential(0,L1range(i))-Potential(0,L1range(i-1)));
       PotentialDiffL2(i-1)=abs(Potential(0,L2range(i))-Potential(0,L2range(i-1)));
    end
    L1Pos=L1range((PotentialDiffL1)==min((PotentialDiffL1)));
    L2Pos=L2range((PotentialDiffL2)==min((PotentialDiffL2)));
    L1=Potential(0,L1Pos(end));
    L2=Potential(0,L2Pos(end));
    if strcmp(Plot,'ON')
        contour(x,y,U_L',[L1 L2]);
    end
    phi=linspace(0,2*pi,3600);
    for i=1:numel(phi)
        surface(i,:)=[R1.*sin(phi(i)) -r_B-R1.*cos(phi(i))];
        SurfF(i,:)=Force(R1.*sin(phi(i)),-r_B-R1.*cos(phi(i)),[0 0]);
        g(i)=sqrt(SurfF(i,1).^2+SurfF(i,2).^2);
        slope(i)=acos((SurfF(i,1)*(-surface(i,1))+SurfF(i,2)*(-surface(i,2)-r_B))./(sqrt(SurfF(i,1).^2+SurfF(i,2).^2)*sqrt(surface(i,1).^2+(-surface(i,2)-r_B).^2)));
    end
    TimeStep=0.1;
    if v==0
        v=100;
    end
    v_0=v;
    for i=1:numel(startangle)
        for j=1:numel(launchangle)
            Position=[(R1+Startheight).*sin(phi(10*round(startangle(i)))) -r_B-(R1+Startheight).*cos(phi(10*round(startangle(i))))];
            velocity=[v_0*sin(pi/180*launchangle(j)) v_0*cos(pi/180*launchangle(j))];
            cntr2=1;
            while (Position(2)-(r2-r_B))^2+Position(1)^2>R2^2 && abs(Position(1))<x(end)
                v_0=cntr2*v;
                Position=[(R1+Startheight).*sin(phi(10*round(startangle(i)))) -r_B-(R1+Startheight).*cos(phi(10*round(startangle(i))))];
                velocity=[v_0*sin(pi/180*launchangle(j)) v_0*cos(pi/180*launchangle(j))];
                cntr=1;
                clear trajectory;
                clear velcurve;
                while (Position(1)^2+(Position(2)+r_B)^2>R1^2 && (Position(2)-(r2-r_B))^2+Position(1)^2>R2^2 && abs(Position(1))<x(end)) || cntr==1
                    acc=Force(Position(1),Position(2),velocity);
                    Position=Position+TimeStep.*velocity;
                    velocity=velocity+TimeStep.*acc;
                    trajectory(cntr,:)=Position;
                    velcurve(cntr,:)=velocity;
                    cntr=cntr+1;
                end
                cntr2=cntr2+1;
            end
            cntr3=1;
            Position=[(R1+Startheight).*sin(phi(10*round(startangle(i)))) -r_B-(R1+Startheight).*cos(phi(10*round(startangle(i))))];
            velocity=[v_0*sin(pi/180*launchangle(j)) v_0*cos(pi/180*launchangle(j))];
            while (Position(2)-(r2-r_B))^2+Position(1)^2>R2^2 && cntr3<11
                v_0=(cntr2-2+0.1*cntr3)*v;
                Position=[(R1+Startheight).*sin(phi(10*round(startangle(i)))) -r_B-(R1+Startheight).*cos(phi(10*round(startangle(i))))];
                velocity=[v_0*sin(pi/180*launchangle(j)) v_0*cos(pi/180*launchangle(j))];
                cntr=1;
                clear trajectory;
                clear velcurve;
                while (Position(1)^2+(Position(2)+r_B)^2>R1^2 && (Position(2)-(r2-r_B))^2+Position(1)^2>R2^2 && abs(Position(1))<x(end)) || cntr==1
                    acc=Force(Position(1),Position(2),velocity);
                    Position=Position+TimeStep.*velocity;
                    velocity=velocity+TimeStep.*acc;
                    trajectory(cntr,:)=Position;
                    velcurve(cntr,:)=velocity;
                    cntr=cntr+1;
                end
                cntr3=cntr3+1;
            end
            cntr4=1;
            Position=[(R1+Startheight).*sin(phi(10*round(startangle(i)))) -r_B-(R1+Startheight).*cos(phi(10*round(startangle(i))))];
            velocity=[v_0*sin(pi/180*launchangle(j)) v_0*cos(pi/180*launchangle(j))];
            while (Position(2)-(r2-r_B))^2+Position(1)^2>R2^2 && cntr4<11
                v_0=(cntr2-2+0.1*(cntr3-2)+0.01*cntr4)*v;
                Position=[(R1+Startheight).*sin(phi(10*round(startangle(i)))) -r_B-(R1+Startheight).*cos(phi(10*round(startangle(i))))];
                velocity=[v_0*sin(pi/180*launchangle(j)) v_0*cos(pi/180*launchangle(j))];
                cntr=1;
                clear trajectory;
                clear velcurve;
                while (Position(1)^2+(Position(2)+r_B)^2>R1^2 && (Position(2)-(r2-r_B))^2+Position(1)^2>R2^2 && abs(Position(1))<x(end)) || cntr==1
                    acc=Force(Position(1),Position(2),velocity);
                    Position=Position+TimeStep.*velocity;
                    velocity=velocity+TimeStep.*acc;
                    trajectory(cntr,:)=Position;
                    velcurve(cntr,:)=velocity;
                    cntr=cntr+1;
                end
                cntr4=cntr4+1;
            end
            vel(i,j)=(cntr2-2+0.1*(cntr3-2)+0.01*(cntr4-1))*v;
            traveltime(i,j)=(cntr*TimeStep)/60;
        end
    end
    if strcmp(Plot,'ON')
        plot(trajectory(:,1),trajectory(:,2),'LineWidth',2);
        hold off;
        figure(2);
        plot(phi.*(180/pi),g,'LineWidth',2);
        xlabel('Latitude (degrees)');
        ylabel('Local g');
        figure(3);
        plot(phi.*(180/pi),slope.*(180/pi),'LineWidth',2);
        xlabel('Latitude (degrees)');
        ylabel('Local apparent slope (degrees)');
    end
    sprintf('Transit velocity (rotating): %g m/s',min(min(vel)))
    sprintf('Travel time: %g min',min(min(traveltime)))
end

% D = MultiRunAnalysis('RunList','KNM1',...         % runlist defines which runs are analysed -> set MultiRunAnalysis.m -> function: GetRunList()
%             'chi2','chi2CMShape',...              % uncertainties: statistical or stat + systematic uncertainties
%             'DataType','Twin',...                 % can be 'Real' or 'Twin' -> Monte Carlo
%             'fixPar','mNu E0 Norm Bkg',...        % free Parameter!!
%             'RadiativeFlag','ON',...              % theoretical radiative corrections applied in model
%             'NonPoissonScaleFactor',1.064,...     % background uncertainty are enhanced
%             'minuitOpt','min ; minos',...         % technical fitting options (minuit)
%             'FSDFlag','SibilleFull',...           % final state distribution
%             'ELossFlag','KatrinT2',...            % energy loss function
%             'SysBudget',24,...                    % defines syst. uncertainties -> in GetSysErr.m;
%             'DopplerEffectFlag','FSD',...
%             'Twin_SameCDFlag','OFF',...
%             'Twin_SameIsotopFlag','OFF',...
%             'SynchrotronFlag','ON',...
%             'AngularTFFlag','OFF',...
%             'TwinBias_Q',18573.73,...
%             'TwinBias_mnuSq',0);
% 
%     M = MultiRunAnalysis('RunList','KNM1',...         % runlist defines which runs are analysed -> set MultiRunAnalysis.m -> function: GetRunList()
%                 'chi2','chi2CMShape',...              % uncertainties: statistical or stat + systematic uncertainties
%                 'DataType','Twin',...                 % can be 'Real' or 'Twin' -> Monte Carlo
%                 'fixPar','E0 Norm Bkg',...    % free Parameter!!
%                 'RadiativeFlag','ON',...              % theoretical radiative corrections applied in model
%                 'NonPoissonScaleFactor',1.064,...     % background uncertainty are enhanced
%                 'minuitOpt','min ; minos',...         % technical fitting options (minuit)
%                 'FSDFlag','SibilleFull',...           % final state distribution
%                 'ELossFlag','KatrinT2',...            % energy loss function
%                 'SysBudget',24,...                    % defines syst. uncertainties -> in GetSysErr.m;
%                 'DopplerEffectFlag','FSD',...
%                 'Twin_SameCDFlag','OFF',...
%                 'Twin_SameIsotopFlag','OFF',...
%                 'SynchrotronFlag','ON',...
%                 'AngularTFFlag','OFF',...
%                 'TwinBias_Q',18573.73,...
%                 'TwinBias_mnuSq',-1);
%             
% R = MultiRunAnalysis('RunList','KNM1',...         % runlist defines which runs are analysed -> set MultiRunAnalysis.m -> function: GetRunList()
%                 'chi2','chi2CMShape',...              % uncertainties: statistical or stat + systematic uncertainties
%                 'DataType','Twin',...                 % can be 'Real' or 'Twin' -> Monte Carlo
%                 'fixPar','E0 Norm Bkg eta',...    % free Parameter!!
%                 'RadiativeFlag','ON',...              % theoretical radiative corrections applied in model
%                 'NonPoissonScaleFactor',1.064,...     % background uncertainty are enhanced
%                 'minuitOpt','min ; minos',...         % technical fitting options (minuit)
%                 'FSDFlag','SibilleFull',...           % final state distribution
%                 'ELossFlag','KatrinT2',...            % energy loss function
%                 'SysBudget',24,...                    % defines syst. uncertainties -> in GetSysErr.m;
%                 'DopplerEffectFlag','FSD',...
%                 'Twin_SameCDFlag','OFF',...
%                 'Twin_SameIsotopFlag','OFF',...
%                 'SynchrotronFlag','ON',...
%                 'AngularTFFlag','OFF',...
%                 'TwinBias_Q',18573.73,...
%                 'TwinBias_mnuSq',-1);
%             
% M.exclDataStart = M.GetexclDataStart(40);
% D.exclDataStart = D.GetexclDataStart(40);
% M.ModelObj.mnuSq_i=4;
% R.ModelObj.mnuSq_i=4;
% R.Fit;
% M.Fit;
% D.Fit;
% 
% %% Global variables
% times = M.ModelObj.qUfrac*M.ModelObj.TimeSec;
% qU    = M.ModelObj.qU; qU    = qU-M.TwinBias_Q; % Energy axis
% NormF = (D.FitResult.par(3+D.ModelObj.nPixels:3+2*D.ModelObj.nPixels-1) + 1)./(M.FitResult.par(3+M.ModelObj.nPixels:3+2*M.ModelObj.nPixels-1) + 1);
% 
% % Spectrum mnu free
% IS  = D.ModelObj.TBDIS; 
% YID = IS./times;
% DIS = D.RunData.TBDIS;
% DIS = DIS./times;
% 
% % Spectrum mnu fixed
% YIM = M.ModelObj.TBDIS; 
% YIM = YIM./times;
% MIS = M.RunData.TBDIS;
% MIS = MIS./times;
% 
% % Spectrum relics
% YIR = R.ModelObj.TBDIS; 
% YIR = YIR./times;
% RIS = R.RunData.TBDIS;
% RIS = RIS./times;
% 
% % Error bar
% err  = (diag(sqrt(D.FitCMShape)));
% err  = err./times;
% err  = err./YIM;
% 
% %% Constraining everything to qULimiteV
% qULimit = -40;
% YID=YID(qU>qULimit);
% DIS=DIS(qU>qULimit);
% MIS=MIS(qU>qULimit);
% YIM=YIM(qU>qULimit);
% RIS=RIS(qU>qULimit);
% YIR=YIR(qU>qULimit);
% err=err(qU>qULimit);
% times=times(qU>qULimit);
% qU=qU(qU>qULimit);
% 
% LocalFontSize = 20;
% 
% fig = figure('Renderer','painters');
% set(fig, 'Units', 'normalized', 'Position', [0.001, 0.001,0.45, 0.8]);
% 
% plot_title = sprintf('Relic neutrino model: \\eta = %.2g',D.ModelObj.eta);
% prlG = [81 126 102]/255;
% prlB = [50 148 216]/255;
% FitStyleArg = {'o','Color','k','LineWidth',1.0,'MarkerFaceColor',rgb('Black'),'MarkerSize',4,'MarkerEdgeColor',rgb('Black')};
% 
% % Ratio
% RSP  = (YIM./DIS);
% RSPd = (YIR./DIS);
% 
% % Plot
% hr1 = plot(qU,RSP,'color',rgb('Salmon'),'LineWidth',3,'LineStyle','-');
% hold on;
% hr2 = plot(qU,ones(1,numel(qU)),'color',prlB,'LineWidth',3,'LineStyle',':');
% hr3 = errorbar(qU,(MIS./DIS),err,FitStyleArg{:},'CapSize',0);
% hr4 = plot(qU,RSPd,'color',rgb('Orange'),'LineWidth',3,'LineStyle','-');
% yl2=ylabel('Ratio');
% katrinsim   = sprintf('\\eta=0, m_{\\nu}^2=0 eV^2');
% sterilemod  = sprintf('Fit: \\eta=%.2g, m_{\\nu}^2=4 eV^2 (fix)',R.ModelObj.eta);
% TwinSim     = sprintf('Twin data: m_{\\nu}^2=-1 eV^2');
% norelics    = sprintf('Fit: \\eta=0 (fix), m_{\\nu}^2=4 eV^2 (fix)');
% hl=legend([hr2 hr3 hr1 hr4],{katrinsim,TwinSim,norelics,sterilemod},'Location','southeast','box','off');
% hl.NumColumns=2;
% hl.FontSize = LocalFontSize-2;
% 
% xlim([min(qU-5) 10]);
% %ylim([min((DIS./YI-1)./err) max((DIS./YI-1)./err)]);
% 
% PRLFormat;
% set(gca,'FontSize',LocalFontSize);
% set(get(gca,'XLabel'),'FontSize',LocalFontSize+4);
% set(get(gca,'YLabel'),'FontSize',LocalFontSize+4); 
% hl.Position(2) = 0.333;
% hold off;