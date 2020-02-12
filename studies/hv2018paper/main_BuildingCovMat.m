addpath(genpath('../../../Samak2.0'));
clear all;

A = init_hvdata_pixel('nPixels', 20);
A.TimeSec = 1;
A.qUfrac  = ones(numel(A.qUfrac),1);
A.SetKrDSBinning(A.qU(1),A.qU(end),1e-2);

initFunc = @fit_hvdata_pixel;
A.HVRipplesP2PV = 0.416;

%MultiCM = zeros(A.nqU*A.nPixels, A.nqU*A.nPixels); %Multipixel Covariance Matrix init

for i= 13:24
fprintf('CovMat for Pixel %u \n',i)
A.FPD_Pixel = i;
MyCovMat = BuildCovMat_HVRipple(A, 1e4, initFunc);
mystr = sprintf('./CovMat/CovMatrix_hv2018paper_Pixel%u.mat', A.FPD_Pixel);
save(mystr,'MyCovMat','-mat');  

%MyCovMat_tmp = load(sprintf('./CovMat/CovMatrix_hv2018paper_Pixel%u.mat', A.FPD_Pixel));
%MyCovMat = struct2cell(MyCovMat_tmp);
%MultiCM((i-1)*A.nqU+1:(i-1)*A.nqU+A.nqU,(i-1)*A.nqU+1:(i-1)*A.nqU+A.nqU) = MyCovMat{:,:};

%mystr = sprintf('./CovMat/MultiCovMatrix_hv2018paper_all%uPixels.mat', A.nPixels);
%save(mystr,'MultiCM','-mat');  
end





