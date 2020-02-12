function [qU , qUfrac] = FlatMTD_S_creator(varargin)
% FlatMTD_S_creator('Range',30,'Save','ON')

% --------------------------- PARSER START ------------------------------%
p = inputParser;
p.addParameter('E0',18575);
p.addParameter('Range',90);
p.addParameter('Save','OFF',@(x)ismember(x,{'ON','OFF'}));
p.parse(varargin{:});
E0     = p.Results.E0;        % eV
Range  = p.Results.Range;     % eV
Save   = p.Results.Save;      %
% --------------------------- PARSER ENDS  ------------------------------%

% MTD Name
TD  = ['FlatSignal' num2str(Range)];

% Build qU vector
if Range <= 30
    qU     = (E0-Range:2:E0-4)';
else
    qU     = [E0-Range:5:E0-35 E0-30:2:E0-4]';
end
%qU     = ([E0-120:5:E0-35 E0-30:1:E0-20 E0-19:0.5:E0+2 E0+5 E0+10 E0+20])';

% Allocate Time / qU
qUfrac = ones(numel(qU),1)./numel(qU);

% Save
switch Save
    case 'ON'
        save([TD '.mat'],'TD','qU','qUfrac','E0');
end

end
