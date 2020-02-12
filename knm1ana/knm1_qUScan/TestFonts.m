figure(1)
f = listfonts;
x = linspace(0,10,numel(f));
y = linspace(0,10000,numel(f));
plot(x,y,'Color',rgb('White'));
for i=2:numel(f)
    hold on;
    try
text(0,y(i),'Samak','FontName',f{i})
    catch
    end
end