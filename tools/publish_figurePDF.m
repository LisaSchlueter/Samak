function publish_figurePDF(fh,filename)

set(fh,'Units','Inches');
pos = get(fh,'Position');
set(fh,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

print(fh,'-dpdf','-r300','-bestfit',filename);
  % system(['pdfcrop -margins 10 ',filename]);
end