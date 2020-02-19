function [parKatrinT2, errKatrinT2] = GetElossPara
d = importdata([getenv('SamakPath'),'inputs/ELossFunction/ELoss_KatrinT2_parametrization.dat']);
amp1 	= d.data(2,1); %3.13912034e-02;
pos1 	= d.data(3,1);%1.19359493e+01;
sig1 	= d.data(4,1);%1.79712352e-01;
amp2 	= d.data(5,1);%2.98192064e-01;
pos2 	= d.data(6,1);%1.28267238e+01;
sig2 	= d.data(7,1);%4.70786732e-01;
amp3 	= d.data(8,1);%7.64887721e-02;
pos3 	= d.data(9,1);%1.49725935e+01;
sig3 	= d.data(10,1);%8.69999541e-01;

amp1err 	= d.data(2,2); %1.33586229e-03;
pos1err 	= d.data(3,2); %9.06581516e-03;
sig1err 	= d.data(4,2); %0.04777842e-03;
amp2err 	= d.data(5,2); %8.74058824e-04;
pos2err 	= d.data(6,2); %2.29909787e-03;
sig2err 	= d.data(7,2); %2.43500480e-03;
amp3err 	= d.data(8,2); %4.38750490e-04;
pos3err 	= d.data(9,2); %4.75464867e-03;
sig3err 	= d.data(10,2); %1.33838628e-02;

mutof 	= 0;  % not used
norm 	= 0;  % not used
mu1 	= 0;  % not used
mu2 	= 0;  % not used
mu3 	= 0;  % not used
%tail 	= 0;  % not used
n1 	    = 0;  %not used
mutoferr 	=  0;   % not used
normerr 	=  0;   % not used
mu1err   	=  0;   % not used
mu2err  	=  0;   % not used
mu3err  	=  0;   % not used
n1err 	    =  0;   % not used

parKatrinT2 = [amp1, pos1, sig1, amp2, pos2, sig2, amp3, pos3, sig3, mutof, norm, mu1, mu2, mu3, n1];
errKatrinT2 = [amp1err, pos1err, sig1err, amp2err, pos2err, sig2err, amp3err, pos3err, sig3err, mutoferr, normerr, mu1err, mu2err, mu3err, n1err];
end