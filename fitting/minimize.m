%% Class with static functions for minimization and model interface.
%
% 2017 Nov - Marc Korzeczek

classdef minimize
    
    methods(Static)
        
        function chi2 = chi2(parValues,fitArg)
            % parValues = vector of nuisance parameters
            % fitArg    = list containing:
            %              - array of experimental data,
            %              - array of uncertainties,
            %              - function handle to update theory, (according to)
            %              - list of parameters (passed to it)
            %              - (optional) function for auxiliary chi2 term
            
           %directly pass mandatory input (1:4)
           %define data and uncertainties arrays
           [data, sigma, theFun, parNames] = fitArg{1:4};
            
           %calculate the array of theoretical values from handle to theory function by
           %parsing fit parameters with the syntax:  'parName',parValue
           parValueCell=num2cell(parValues); 
          
           if numel(parValues)==numel(parNames) %single pixel fit
           parArgIn={parNames{:}; parValueCell{:}}; 
           else   %multipixel fit;
           nPixel = (numel(parValueCell)-2)/2; %nPixel for Krypon (2 common, (nPixel*2) shared parameter)       
           parArgIn={parNames{1:2}; parValueCell{1:2}};
           parArgIn={parArgIn{:} parNames{3} [parValueCell{3:2+nPixel}] parNames{4} [parValueCell{3+nPixel:2+2*nPixel}]};                                
           end
           
           theory= theFun(parArgIn{:});
           % Calculate auxiliary chi2 value if and auxiliary function handle is provided
           auxChi2=0;
           if numel(fitArg)==5
               auxFun=fitArg{5};
               auxChi2=auxFun(parValueCell{:});
           end            
            
            % Calculate chi2 value w/ provided uncertainties array or covariance matrix   
            chi2 = sum((data - theory)' * ( sigma \ (data - theory) )) ;         
            
            %return final chi2 value
            chi2=chi2 +auxChi2;
        end
        
    end
    
end
