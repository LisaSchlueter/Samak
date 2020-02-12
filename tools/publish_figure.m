function publish_figure(fh,filename,fixLSeps)
if nargin<3
    fixLSeps = false;
end

if ~fixLSeps
    print(fh,'-depsc2','-loose','-r300',filename);
    
else
    print(fh,'-depsc2','-loose','-r300','tmp.eps');
    fixPSlinestyle('tmp.eps',filename);
    system('rm tmp.eps');
end
end