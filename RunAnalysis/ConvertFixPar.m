% little auxiliary script to convert
% input: string, which contains parameters name which shall be free(!) e.g. 'mNu, E0, B, N'
% output: string of parameter which are fixed, which can be used in RunAnalysis class, Fit class,...
function fixPar = ConvertFixPar(varargin)
p=inputParser;

p.addParameter('freePar','',@(x)ischar(x));
p.addParameter('nPixels',1,@(x)isfloat(x)); % number of pixels or rings
p.addParameter('nPar','',@(x)isfloat(x)); % number of maximal available parameter is fit
p.addParameter('Mode','Normal',@(x)ismember(x,{'Normal','Reverse'}));
p.parse(varargin{:});

freePar = p.Results.freePar;
nPixels = p.Results.nPixels;
nPar    = p.Results.nPar;
Mode    = p.Results.Mode;
% ------------parser end -------------------

if strcmp(Mode,'Normal')
    if isempty(nPar)
        nPar = 3*nPixels+10;
    end
    
    fixPar_v = 1:nPar; % array of all parameter
    
    % if neutrino mass is free parameter, erase from list of fixed parameter
    if contains(freePar,'mNu') % same for whole FPD (==uniform)
        findIndex = ismember(fixPar_v,1);
        fixPar_v(findIndex)=[];
    end
    
    % endpoint
    if contains(freePar,'E0') % same for whole FPD (==uniform)
        findIndex = ismember(fixPar_v,2);
        fixPar_v(findIndex)=[];
    end
    
    % background
    if contains(freePar,'Bkg') %ringwise
        findIndex = ismember(fixPar_v,3:2+nPixels);
        fixPar_v(findIndex)=[];
    end
    
    % signal normalization
    if contains(freePar,'Norm') %ringwise
        findIndex = ismember(fixPar_v,3+nPixels:3+2*nPixels-1);
        fixPar_v(findIndex)=[];
    end
    
    % FSD: final state distribution ground state-excited state normalization factors
    if contains(freePar,'FSD')  % same for whole FPD (==uniform)
        findIndex = ismember(fixPar_v,3+2*nPixels:2*nPixels+8);
        fixPar_v(findIndex)=[];
    end
    
    % qU Offset
    if contains(freePar,'qU') %ringwise
        % fix qU of the first ring to zero -> anchor point
        findIndex = ismember(fixPar_v,2*nPixels+10:3*nPixels+8); %2*nPixels+9
        fixPar_v(findIndex)=[];
    end
    
    % background slope
    if contains(freePar,'BkgSlope')  % same for whole FPD (==uniform)
        findIndex = ismember(fixPar_v,3*nPixels+9);
        fixPar_v(findIndex)=[];
    end
    
    % ringwise effective neutrino mass --> sigma^2 energy smearing of DS
    if contains(freePar,'mTSq')  % ringwise
        % first neutrino mass offset (3*nPixels+10) fixed to zero -> anchor point
        findIndex = ismember(fixPar_v,3*nPixels+11:4*nPixels+9);
        fixPar_v(findIndex)=[];
    end
    
    % Fraction T- ions
    if contains(freePar,'FracTm')  % ringwise
        findIndex = ismember(fixPar_v,4*nPixels+10);
        fixPar_v(findIndex)=[];
    end
    fixPar = sprintf('fix %.0f ;',fixPar_v);
    
elseif strcmp(Mode,'Reverse')
  %convert fixPar back to readable string
  %input: (freePar)  given to fit class e.g. 'fix 1 ;fix 5 ;' --> actual parameters which are fixed!
  %output: fixPar (actual parameters which are free );
  fixPar = '';
  
  if ~contains(freePar,'fix 1 ;')
      fixPar = [fixPar,'mNu'];
  end
  
  if ~contains(freePar,'fix 2 ;')
      fixPar = [fixPar,'E0'];
  end
  
  if ~contains(freePar,arrayfun(@(x) sprintf('fix %.0f ;',x),3:2+nPixels,'UniformOutput',0))
      fixPar = [fixPar,'Bkg'];
  end
  
  if ~contains(freePar,arrayfun(@(x) sprintf('fix %.0f ;',x),3+nPixels:3+2*nPixels-1,'UniformOutput',0))
      fixPar = [fixPar,'Norm'];
  end
  
  if ~contains(freePar,arrayfun(@(x) sprintf('fix %.0f ;',x),3+2*nPixels:2*nPixels+8,'UniformOutput',0))
      fixPar = [fixPar,'FSD'];
  end
  
  if ~contains(freePar,arrayfun(@(x) sprintf('fix %.0f ;',x),2*nPixels+10:3*nPixels+8,'UniformOutput',0)) && nPixels>1 % multiring
      fixPar = [fixPar,'qU'];
  end
  
  if ~contains(freePar,sprintf('fix %.0f ;',3*nPixels+9))
      fixPar = [fixPar,'BkgSlope'];
  end
  
  if ~contains(freePar,arrayfun(@(x) sprintf('fix %.0f ;',x),3*nPixels+11:4*nPixels+9,'UniformOutput',0)) && nPixels>1 % multiring
      fixPar = [fixPar,'mTSq'];
  end
  
  if ~contains(freePar,sprintf('fix %.0f ;',4*nPixels+10))
      fixPar = [fixPar,'FracTm'];
  end
  
end
end