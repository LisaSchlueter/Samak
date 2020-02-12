function sim_hvdata_allpixels(varargin)
%
%     Simulation INTEGRAL SPECTRUM 83mKr Lines
%          All Pixels Simultaneously
%
%  For each pixels
%  - Background fixed with pixel-wise fit : 148 par. fixed
%  - Amplitude fixed with pixel-wise fit : 148 par. fixed
%  A common fit for: E0 and W
%
%          Th. Lasserre - CEA Saclay
%                February 2018
%

% Parser
p = inputParser;
p.addParameter('display','ON',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('nsim',100);
p.parse(varargin{:});
display    =    p.Results.display;
nsim       =    p.Results.nsim;

% Initalization
A = init_hvdata_pixel('FPD_Pixel',0,'nPixels',40);
A.TimeSec = 1;
A.qUfrac  = ones(numel(A.qUfrac),1);
A.SetKrDSBinning();
A.DisplayKrInfo;
A.L3_32_Phi0allPixels    = [44.6232      47.1541      46.3954      46.6741       45.761      46.6343      47.0666      45.9541      46.2233      45.9166      45.9511      45.5604      47.8608      45.5325      45.3468      46.2513      47.6844      48.7293      47.6418      46.6697      46.9427      50.1574      46.6698       44.659      46.4393       46.455      47.6318      47.1173       46.598      45.5871      46.1051       47.319      46.3033      46.8243      46.2841      47.1165      45.7393       46.505      44.4374      47.6845 ];
A.BKG_RateSecallPixels   = [21.646      20.1148      20.1998      20.3904      21.3103      20.2357      19.2365      20.2849      20.4684      20.9596      20.6376       21.125      19.9332      20.5295      21.5831      20.7355      20.4928      19.1473      19.2369      19.6851       19.459      18.6784       20.078      20.5546      20.7442      20.4938       19.193      20.4008      19.9845      21.0239      19.9944      19.3057      20.1331      19.8468      20.8075      19.8936      20.3485      20.5098      21.5158      19.6338  ];
A.L3_32_E_i = 30472.569  ;
A.L3_32_W_i = 1.149      ;
A.ComputeKrDSallPixels();

    switch display
        case 'ON'
            figure(1)
    end

%progressbar('...Samak Simulate Krypton IS All Pixels...');
for j=1:1:nsim    
    A.ComputeKrISallPixels();

    switch display
        case 'ON'
            for i=1:1:A.nPixels
                %fprintf('Simulation: %.0f   \t Pixel: %.0f \n',j,i);
                tmpis = (A.KrISallPixels(:,i));
                %disp([A.qU',tmpis]);
                plot(A.qU',tmpis,'Color','Black','LineWidth',1,'LineStyle','-');
                hold on
            end
            hold off
            grid on
            xlabel('qU (eV)','FontSize',10);
            ylabel('Counts','FontSize',10);
    end
    %progressbar(j/nsim);

end

end