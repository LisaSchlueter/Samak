function Write2Txt(varargin)
p=inputParser;
p.addParameter('filename','',@(x)ischar(x));
p.addParameter('variable','',@(x)isfloat(x));
p.addParameter('variableName','',@(x)ischar(x));
p.addParameter('nCol',1,@(x)isfloat(x) & x>0);
p.parse(varargin{:});
filename      = p.Results.filename;
variable      = p.Results.variable;
variableName  = p.Results.variableName;
nCol          = p.Results.nCol;

system(['touch ',filename,'.txt'])

fileID = fopen([filename,'.txt'],'w');
fprintf(fileID,'%s\n',variableName);
if nCol==1
    fprintf(fileID,'%.5g\n',variable);
elseif nCol==2
    fprintf(fileID,'%.20e %12.20e\n',variable);
elseif nCol==3
    fprintf(fileID,'%.3f %.3f %.3f\n',variable);
end
fclose(fileID);
fprintf('variables writen to txt file: %s.txt \n',filename)
end
