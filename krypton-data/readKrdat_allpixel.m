% 
% Read Krypton Data
% Th. Lasserre, August 2017
% 
% Input   : Line 
% Output for 1 Line ALL pixels: 
%     * retarding energy qU in eV
%     * count rate in cps
%     * uncertainty of count rate in cps
% 
% 
% - spectra of L3-32 line --> 33051_L3-32 folder
% - spectra of K-32 line  --> 33054_K-32 folder
% 
% - Pixel-wise (148): 
%     * L3-32: 33051_channelXXX.txt
%     * K32  : 33054_channelXXX.txt
% 


%Line  = 'L3-32';
Line  = 'K32';
pub   = 0;
fign  = 1;

switch Line
    case 'L3-32'
folder = '33051_L3-32';
prefix = '33051';
    case 'K32'
folder = '33054_K-32';
prefix = '33054';
end

krmap = zeros(148,2,31);

for(i=0:1:147)

Pixel = num2str(i);
krfile = sprintf('./%s/%s_channel%s.txt',folder,prefix,Pixel);
krplot = sprintf('./figs/%s_channel%s.eps',folder,prefix,Pixel);
krtitle = sprintf('%s Line - Pixel %s',Line,Pixel);
data = importdata(krfile);

% read qU
if (i==0)    qU = data(:,1);
end

% read pixel
krmap(i+1,1,:)    = data(:,2);
krmap(i+1,2,:)    = data(:,3);
end


figure(fign)
for(j=1:1:148)
    
    krtitle = sprintf('%s Line - Pixel %s',Line,num2str(j));

    Count    = flip(krmap(j,1,:));
    CountErr = flip(krmap(j,2,:));
    
    h = errorbar(qU(:),Count(:),CountErr(:),...
        'ks','MarkerSize',3,'MarkerFaceColor',.8*[1 1 1],'LineWidth',1);
    title(krtitle,'FontSize',14);
    a = legend(h,'Krypton Data');% legend(a,'boxoff');
 
    pause(1)
end
grid on
set(gca, 'YScale', 'lin');
xlabel('qU (eV)','FontSize',14);
ylabel('Counts per Energy Bin','FontSize',14);
set(gca,'FontSize',12);
PrettyFigureFormat;
