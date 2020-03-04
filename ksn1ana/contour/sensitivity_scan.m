%% KATRIN Sterile neutrino sensitivity

%% Settings
% constants
range = 40; % eV below the endpoint
d=3;

start_decade_m = -1; start_decade_s = -2;
stop_decade_m = 4; stop_decade_s = 0;

% variables
m4=[];
sith4=[];

for n=start_decade_m:stop_decade_m-2
    m4=[m4,logspace(n,n+1,d)];
    m4(length(m4))=[];
end
m4=[m4,logspace(stop_decade_m-1,stop_decade_m,d)];

for n=start_decade_s:stop_decade_s-2
    sith4=[sith4,logspace(n,n+1,d)];
    sith4(length(sith4))=[];
end
sith4=[sith4,logspace(stop_decade_s-1,stop_decade_s,d)];

chiX  = zeros(length(sith4),length(m4));

%% For loop
km=0;
for m=m4
    m
    km=km+1;
    ks=0;
    for s=sith4
        s
        ks=ks+1;
        X2=fit_chi(m,s);
        
        chiX(ks,km)=X2;
    end
end

%% sin(2th)
sith42=[];

for elt=sith4
    sith42=[sith42,1-(1-2*elt*elt)^2];
end

%% Plot

contour(sith4,m4,chiX',[4.61 4.61])
%contour(sort(sith42),m4,chiX',[4.61 4.61])
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');