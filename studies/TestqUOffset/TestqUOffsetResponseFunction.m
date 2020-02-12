%% Test how qU-Offset can be implemented best
M =  MultiRunAnalysis('RunList','KNM1_m149mvRW','AnaFlag','Ring','RingList',7:9,'PullFlag',5);
RF_i = M.ModelObj.RF;
TBDDSandRF_i = reshape(repmat(M.ModelObj.TBDDS,1,M.ModelObj.nqU,1),[M.ModelObj.nTe,M.ModelObj.nqU,M.ModelObj.nPixels]).*RF_i;
TBDIS_i = squeeze(simpsons(M.ModelObj.Te,TBDDSandRF_i));

qUOffset = [0.01 0.1 5];
%% Method 1: Shift Response Function by Interpolation of Response Function, keep Te as it is
thisqUOffset = 3; % for plot

%interpolate RF
RF_Offset = zeros(M.ModelObj.nTe,M.ModelObj.nqU,numel(qUOffset));
for i=1:numel(qUOffset)
    RF_Offset(:,:,i) = interp1(M.ModelObj.qU(:,i),squeeze(M.ModelObj.RF(:,:,i))',M.ModelObj.qU(:,i)+qUOffset(i),'spline')';
end

% for i=1:numel(qUOffset)
%     RF_Offset(:,:,i) = interp1(M.ModelObj.Te,squeeze(M.ModelObj.RF(:,:,i)),M.ModelObj.Te+qUOffset(i),'spline');
% end

for q=1:M.ModelObj.nqU
RF_Offset(M.ModelObj.Te<=M.ModelObj.qU(q),q,:) =0;
end

RF_Offset(RF_Offset<0)=0;

%calculate TBDIS with RF
TBDDSandRF_method1 = reshape(repmat(M.ModelObj.TBDDS,1,M.ModelObj.nqU,1),[M.ModelObj.nTe,M.ModelObj.nqU,M.ModelObj.nPixels]).*RF_Offset;
TBDIS_method1 = squeeze(simpsons(M.ModelObj.Te,TBDDSandRF_method1));

% plot response function: shift versus interpolation: Result: Interpolation of RF not good, finer binning would be required
plot(M.ModelObj.Te+qUOffset(thisqUOffset),squeeze(RF_i(:,5,thisqUOffset)));
hold on;
plot(M.ModelObj.Te,squeeze(RF_Offset(:,5,thisqUOffset)));
hold off;
legend('shifted','interpolated');
%% Method 2: Sfift Response Function by Shifting Te-vector and interpolating TBDDS at qU=qU+Offset

%interpolate TBDDS
TBDDS = zeros(M.ModelObj.nTe,numel(qUOffset));
for i=1:1:numel(qUOffset)
    TBDDS(:,i) = interp1(M.ModelObj.Te,M.ModelObj.TBDDS(:,i),M.ModelObj.Te+qUOffset(i));
end
TBDDS(TBDDS<0)=0;
TBDDS(isnan(TBDDS))=0;

% compute TBDIS with shifted Te vector and initial RF
TBDDSandRF_method2 = reshape(repmat(TBDDS,1,M.ModelObj.nqU,1),[M.ModelObj.nTe,M.ModelObj.nqU,M.ModelObj.nPixels]).*RF_i;
TBDIS_method2 = zeros(M.ModelObj.nqU,M.ModelObj.nPixels);
for i=1:M.ModelObj.nPixels
TBDIS_method2(:,i) = squeeze(simpsons(M.ModelObj.Te+qUOffset(:,i),TBDDSandRF_method2(:,:,i)));
end

% % plot
plot(M.ModelObj.Te-qUOffset(thisqUOffset),M.ModelObj.TBDDS(:,thisqUOffset));
hold on;
plot(M.ModelObj.Te,TBDDS(:,thisqUOffset));
hold off;
legend('shifted','interpolated');  
%% Compare TBDIS from two methods
plot(M.ModelObj.qU,(TBDIS_method1(:,thisqUOffset)-TBDIS_method2(:,thisqUOffset))./sqrt(TBDIS_method1(:,thisqUOffset)));
ylabel('norm. residuals (\sigma)'); xlabel('qU (eV)');

% plot(M.ModelObj.qU,TBDIS_method1(:,thisqUOffset));%
% hold on;
% plot(M.ModelObj.qU,TBDIS_method2(:,thisqUOffset));%
% legend('inerpolate RF','shift Te, interpolate TBDDS');
% ylabel('counts'); xlabel('qU (eV)');
% hold off;
%% Compare Product of TBDDS*RF
plot(M.ModelObj.Te,TBDDSandRF_method1(:,5,thisqUOffset));
hold on;
plot(M.ModelObj.Te+qUOffset(thisqUOffset),TBDDSandRF_method2(:,5,thisqUOffset));
plot(M.ModelObj.Te,TBDDSandRF_i(:,5,thisqUOffset));
hold off;                   
legend('inerpolate RF','shift Te, interpolate TBDDS','initial');                    
ylabel('TBDDS * RF');                    
                    
%% Compare Response Functions
plot(M.ModelObj.Te,RF_Offset(:,5,thisqUOffset)); hold on;
plot(M.ModelObj.Te+qUOffset(thisqUOffset),RF_i(:,5,thisqUOffset)); 
plot(M.ModelObj.Te,RF_i(:,5,thisqUOffset));
legend('inerpolate RF','shift Te','initial');  
hold off;

%% COmpare TBDIS
plot(M.ModelObj.qU(:,thisqUOffset),M.ModelObj.TBDIS(:,thisqUOffset));
hold on;
plot(M.ModelObj.qU(:,thisqUOffset),TBDIS_method2(:,thisqUOffset));
hold off;
                    
                    