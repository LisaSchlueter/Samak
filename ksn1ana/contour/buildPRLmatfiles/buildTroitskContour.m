%% 95CL
%SinSquareTheta_X  = [0.95586     0.67306     0.30156     0.10686    0.027868   0.0083969   0.0072729   0.0078932   0.0085653   0.0071548   0.0045174];
SinSquareTheta_X  = [0.5 0.5     0.30156     0.10686    0.027868   0.0083969   0.0072729   0.0078932   0.0085653   0.0071548   0.0045174];
m4Square_Y        = [1.00441       2.0228       3.0128       5.0125      10.0436      19.9489       29.594       39.828       49.933       69.924       98.775].^2;
[i j]=sort(m4Square_Y);
a=sort(m4Square_Y);
b=(SinSquareTheta_X(j));
m4Square_Y=a;
SinSquareTheta_X=b;

SinSquare2Theta_X  = 4*SinSquareTheta_X.*(1-SinSquareTheta_X);
DmSquare41_Y       = m4Square_Y;

CL = '95';
Reference = 'https://link.springer.com/content/pdf/10.1134/S0021364013020033.pdf';
save('coord_Troitsk_95CL.mat','SinSquare2Theta_X','DmSquare41_Y','SinSquareTheta_X','m4Square_Y','CL','Reference');
