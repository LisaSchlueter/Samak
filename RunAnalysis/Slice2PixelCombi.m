function [PixList,SlicePixList,SliceAngPos] = Slice2PixelCombi(myPixList,RingMerge)
% function to synchronize list of pixel (myPixList) with list of azimuthal patches  
% input: - myPixList: list of  pixel, which you want to analyze
%        - myRingList:list of rings, which you want to analyze
% output: - PixList: list of selected pixels, which are within the selected rings
%         - RingPixList: cell of rings, which shows which pixel belongs to which rings

filen = [getenv('SamakPath'),'RunAnalysis/PixelAngPosition.mat'];
if exist(filen,'file')
   load(filen,'PixPos_deg');
else
   [~,~,PixPos_deg] = FPDViewer(zeros(148,1),'ReDrawSkeleton','ON');
   PixNum = 1:148;
  
   save(filen,'PixPos_deg','PixNum');
end

% find pixels that have the same angular position
PixPos_deg_u = unique(PixPos_deg);

nSlice_i = numel(PixPos_deg_u); % original number of slices
PixList_all = 1:148;
Pix_plot = NaN.*ones(148,1);

if strcmp(RingMerge,'Slice') 
    nSlice = nSlice_i;
    SlicePixList = cell(nSlice,1);
    SliceAngPos  = zeros(nSlice,1);
    for i=1:nSlice
        SlicePixList{i} = PixList_all(PixPos_deg==PixPos_deg_u(i));
        Pix_plot(SlicePixList{i}) = i; % for sanity plot
        SliceAngPos(i,:) = PixPos_deg_u(i);
    end
elseif strcmp(RingMerge,'Slice2')
    nSlice = nSlice_i/2;
    SlicePixList = cell(nSlice,1);
    SliceAngPos  = zeros(nSlice,2);
    for i=1:nSlice
        SlicePixList{i} = PixList_all(PixPos_deg==PixPos_deg_u(2*i-1) | PixPos_deg==PixPos_deg_u(2*i));
        Pix_plot(SlicePixList{i}) = i; % for sanity plot
        SliceAngPos(i,1) = PixPos_deg_u(2*i-1);
        SliceAngPos(i,2) = PixPos_deg_u(2*i);
    end
elseif strcmp(RingMerge,'Slice3')
    % stack 3 into 1
    nSlice = nSlice_i/3;
    SlicePixList = cell(nSlice,1);
    SliceAngPos  = zeros(nSlice,3);
    for i=1:nSlice
        SlicePixList{i} = PixList_all(PixPos_deg==PixPos_deg_u(3*i-2) | PixPos_deg==PixPos_deg_u(3*i-1)...
                                    | PixPos_deg==PixPos_deg_u(3*i));
        Pix_plot(SlicePixList{i}) = i; % for sanity plot
        SliceAngPos(i,1) = PixPos_deg_u(3*i-2);
        SliceAngPos(i,2) = PixPos_deg_u(3*i-1);
        SliceAngPos(i,3) = PixPos_deg_u(3*i);
    end    
elseif strcmp(RingMerge,'Slice4')
    % stack 4 slides into 1
    nSlice = nSlice_i/4;
    SlicePixList = cell(nSlice,1);
    SliceAngPos  = zeros(nSlice,4);
    for i=1:nSlice
        SlicePixList{i} = PixList_all(PixPos_deg==PixPos_deg_u(4*i-3) | PixPos_deg==PixPos_deg_u(4*i-2)...
                                    | PixPos_deg==PixPos_deg_u(4*i-1) | PixPos_deg==PixPos_deg_u(4*i));
        Pix_plot(SlicePixList{i}) = i; % for sanity plot
        SliceAngPos(i,1) = PixPos_deg_u(4*i-3);
        SliceAngPos(i,2) = PixPos_deg_u(4*i-2);
        SliceAngPos(i,3) = PixPos_deg_u(4*i-1);
        SliceAngPos(i,4) = PixPos_deg_u(4*i);
    end
end

%  FPDViewer(Pix_plot) % sanity plot
%  colormap(hsv)
%%
 SlicePixList = cellfun(@(x) intersect(myPixList,x),SlicePixList,'UniformOutput',0); %cell, which shows which pixels belong to which ring
 PixList = intersect(myPixList,PixList_all);
end