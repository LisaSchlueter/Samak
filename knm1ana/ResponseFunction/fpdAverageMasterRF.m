function [MasterTe MasterRFfpd] = fpdAverageMasterRF()
% Read Reference Knm1 Response Function
% Average over pixel list
%
% fpdRF:
%  - column = qu
%  - line   = te
%
% te: kinetic energy
% fpdRF: matrix providing rf value
%           column = pixel
%           line   = qu
% 
PixList = [1    2    3    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19   20   21   22   23   24   25   26   27   28   29   30   31   32   33   34   35   36   37   38   39   40   41   42   43   44   45   46   47   48   49   50   51   52   53   54   55   56   57   58   59   60   61   62   63   64   65   66   67   68   69   70   71   72   73   74   75   76   77   78   79   80   81   82   83   84   85   86   87   88   89   90   91   92   93   94   95   96   97  100  102  103  104  105  106  107  108  109  110  114  115  116  117  118  119  120  121  132  133  134];

[MasterTe MasterRF] = readMasterRF();

MasterRFfpd         = squeeze(mean(MasterRF(:,:,PixList),3))';

MasterTe            = mean(MasterTe(:,PixList),2);

fign = figure('Name','Samak','NumberTitle','off','rend','painters','pos',[10 10 1200 600]);
maintitle=sprintf('Master Response Function - \\rho.d = %.4g mol/cm^2 - Bs = %.2f T - Ba = %.3g G',...
    1.1e17,2.52,6.3);
a=annotation('textbox', [0 0.9 1 0.1], ...
    'String', maintitle, ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center');
a.FontSize=20;a.FontWeight='bold';
h = plot(MasterTe,MasterRFfpd(:,:),'LineWidth',2);
set(gca, 'YScale', 'lin');
xlabel('T_e (eV)','FontSize',18);
ylabel('transmission','FontSize',18);
grid on
PrettyFigureFormat
set(gca,'FontSize',18);
end