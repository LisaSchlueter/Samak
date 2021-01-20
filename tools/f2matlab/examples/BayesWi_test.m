function bayes(varargin)


clear all; %clear functions;

global unit2fid;  if ~isempty(unit2fid), unit2fid=[]; end
persistent maxit reject tepoch xref ; 

format_10=['abc','%4x', '\n ' ,repmat([repmat(['%20.13e'] ,1,6), '\n ' , '\n ' ] ,1,2), '\n ' , '\n ' ];
format_11=['%20.13e','tepoch=','%20.13e','tepoch=','%20.13e'];
format_6=[ '\n ' , '\n ' ,'%2x','BAYES FILTER', '\n ' ,'%2x','epoch time :','%20.13e', '\n ' ,'%2x','Previous Estimated State Vector:', '\n ' ,'%2x',repmat('%20.13e',1,3), '\n ' ,'%2x',repmat('%20.13e',1,3), '\n ' ,'%2x','maximum iterations:','%3d', '\n ' ,'%2x','reject if gt ','%12.5e',' sigma ', '\n ' , '\n ' , '\n ' ];

 if isempty(tepoch), tepoch=0; end;
 if isempty(xref), xref=zeros(1,6); end;
 if isempty(maxit), maxit=0; end;
 if isempty(reject), reject=0; end;

tepoch=sqrt(2.0);
xref(:)=[sqrt(3.0),2.,3.,4.,5.,6.];
maxit=sqrt(5.0);
reject=sqrt(6.0);

thismlfid=fopen(strtrim('baydebug.out'),'r+');
  unit2fid=[unit2fid;9,thismlfid];


[writeErrFlag]=writeFmt(1,[format_11],'tepoch','tepoch','tepoch');


[writeErrFlag]=writeFmt(9,[format_10],'xref','xref');


[writeErrFlag]=writeFmt(9,['abc','%4x', '\n ' ,repmat([repmat(['%20.13e'] ,1,6), '\n ' , '\n ' ] ,1,2), '\n ' , '\n ' ],'xref','xref');



[writeErrFlag]=writeFmt(9,[format_6],'tepoch','xref','maxit','reject');


[writeErrFlag]=writeFmt(1,['%g','%c'],'tepoch','''ggggggg''');

[writeErrFlag]=writeFmt(9,['%g','%c','%g','%c','%g'],'tepoch','''asd''','xref','''tytyty''','maxit');
[writeErrFlag]=writeFmt(9,['%g'],'tepoch');

[writeErrFlag]=writeFmt(1,['%c'],'''wrote it''');

try;
fclose(unit2fid(find(unit2fid(:,1)==9,1,'last'),2));
unit2fid=unit2fid(find(unit2fid(:,1)~=9),:);
end;


[writeErrFlag]=writeFmt(1,['%g','%c','%g','%c','%g'],'tepoch','''asd''','xref','''tytyty''','maxit');

[writeErrFlag]=writeFmt(1,[format_6],'tepoch','xref','maxit','reject');




end %program bayes

