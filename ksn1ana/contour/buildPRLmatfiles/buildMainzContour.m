%% 95CL
SinSquare2Theta_X  = [0.051704    0.056937    0.073569     0.07705     0.07963    0.085703    0.094552     0.12715     0.12799     0.14072     0.14606     0.14965     0.17967     0.21786     0.25137     0.26213      0.3912     0.42176     0.52752     0.71515     0.74621     0.87553      0.9947     0.99574];
DmSquare41_Y       = [20199           29956           31595          3033.1            9894          1980.8           32970          1001.8          5111.6           900.6           775.8           34041          200.56          396.56          102.52           35903           37069          69.124           37867           38682          30.439          19.878           10.38           39938];
[i j]=sort(DmSquare41_Y);
a=sort(DmSquare41_Y);
b=(SinSquare2Theta_X(j));
DmSquare41_Y=a;
SinSquare2Theta_X=b;
CL = '95';
Reference = 'European Physical Journal C volume 73, Article number: 2323 (2013)';
save('coord_Mainz_95CL.mat','SinSquare2Theta_X','DmSquare41_Y','CL','Reference')

%% 90CL
SinSquare2Theta_X = [0.039382    0.042451    0.066352    0.067811    0.071662    0.080131    0.081241    0.091953     0.12579     0.14282     0.14885     0.19043     0.19898      0.3247     0.43074      0.5109       0.629     0.98938      1.0011];
DmSquare41_Y      = [20415           30277          3033.1          5003.9          55.269           32970          1780.6          7111.6           872.2           202.7           35903            97.2          465.26          71.369          47.108           39096          28.554          10.053           40366];
[i j]=sort(DmSquare41_Y);
a=sort(DmSquare41_Y);
b=(SinSquare2Theta_X(j));
DmSquare41_Y=a;
SinSquare2Theta_X=b;
CL = '90';
Reference = 'European Physical Journal C volume 73, Article number: 2323 (2013)';
save('coord_Mainz_90CL.mat','SinSquare2Theta_X','DmSquare41_Y','CL','Reference')