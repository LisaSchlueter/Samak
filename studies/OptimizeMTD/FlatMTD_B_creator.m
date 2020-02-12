function [qU , qUfrac] = FlatMTD_B_creator(varargin)
% FlatMTD_B_creator('Range',30,'Save','ON')

% --------------------------- PARSER START ------------------------------%
p = inputParser;
p.addParameter('E0',18575,@(x)isfloat(x));
p.addParameter('BqU',[5 10 20]);
p.addParameter('Save','OFF',@(x)ismember(x,{'ON','OFF'}));
p.parse(varargin{:});
E0     = p.Results.E0;        % eV
BqU    = p.Results.BqU;       % 
Save   = p.Results.Save;      %
% --------------------------- PARSER ENDS  ------------------------------%

% MTD Name
TD  = ['FlatBackground'];

%qU  = E0+BqU;
qU   = [E0+BqU];

% Allocate Time / qU
qUfrac = ones(numel(qU),1)./numel(qU);

% Save
switch Save
    case 'ON'
        save([TD '.mat'],'TD','qU','qUfrac','E0');
end

end
