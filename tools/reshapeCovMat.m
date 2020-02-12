NqU      = 40;
Nring    = 4;
qUstart  = 14;
qUstop   = NqU;
nqU_used = NqU-qUstart+1;

% Cut matrix in blocs
M=a.obj.CovMat;
M4=mat2tiles(M,[NqU NqU]); % 40qU values

% exclude data
for i=1:1:numel(M4)
MTMP   = M4{i};    
M4C{i} = MTMP(qUstart:qUstop,qUstart:qUstop);
end

% rebuild matrix with excluded data
M4CN = [M4C{1} M4C{5} M4C{9} M4C{13};
        M4C{2},M4C{6},M4C{10},M4C{14};
        M4C{3},M4C{7},M4C{11},M4C{15};
        M4C{4},M4C{8},M4C{12},M4C{16}];
    

% check lisa
exclIndex = sort(reshape(repmat(qUstart:qUstop,[Nring,1])+[0:Nring-1]'.*NqU,[nqU_used*Nring,1]));

figure(1)
subplot(1,2,1)
imagesc(M4CN);colorbar
subplot(1,2,2)
imagesc(M4CN-M(exclIndex,exclIndex));colorbar