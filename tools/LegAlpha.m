function LegAlpha(leg,alpha)

legend boxon
leg.EdgeColor = 'none';
set(leg.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;alpha]));
         
end