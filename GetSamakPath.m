function GetSamakPath
% can the relative path to Samak3.0 directory, according to your current position
s = what('Samak3.0');
SamakAbsPath = s.path; %absolute path of Samak folder
SamakRelPath = relativepath(SamakAbsPath,pwd);
setenv('SamakPath',SamakRelPath);
end
