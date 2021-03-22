function Write2Txt(varargin)
p=inputParser;
p.addParameter('filename','',@(x)ischar(x));
p.addParameter('variable','',@(x)isfloat(x));
p.addParameter('variableName','',@(x)ischar(x));
p.addParameter('nCol',1,@(x)isfloat(x) & x>0);
p.addParameter('Format','txt',@(x)ismember(x,{'txt','dat'}));
p.parse(varargin{:});
filename      = p.Results.filename;
variable      = p.Results.variable;
variableName  = p.Results.variableName;
nCol          = p.Results.nCol;
Format        = p.Results.Format;
system(['touch ',filename,'.',Format])

fileID = fopen([filename,'.',Format],'w');
fprintf(fileID,'%s\n',variableName);
if nCol==1
    fprintf(fileID,'%.5g\n',variable);
elseif nCol==2
    fprintf(fileID,'%.10e %12.10e\n',variable);
elseif nCol==3
    fprintf(fileID,'%.10e %.10e %.10e\n',variable);
elseif nCol==8
    fprintf(fileID,'%.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f\n',variable);
end
fclose(fileID);
fprintf('variables writen to txt file: %s.%s \n',filename,Format)
end
