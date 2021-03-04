function PrettyLegendFormat(leg,varargin)
p = inputParser;
p.addParameter('alpha',0.7,@(x)isfloat(x));
p.parse(varargin{:});
alpha = p.Results.alpha;
leg.EdgeColor = rgb('White');
set(leg.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;alpha])); % make legend semi-transparent
end