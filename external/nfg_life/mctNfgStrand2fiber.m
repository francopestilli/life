function fg  = mctNfgStrand2fiber(strands)
% Transforms NFG strands to a MicroTrack fiber group 
%
%  fg = mctNfgStrand2fiber(strands)
%
%  strands is a collection of NFG strands (fibers) returned by
%  mctNfgLoadStrands.m
%
% by Franco Pestilli
%
% see also:
%     mctNfgLoadStrands.m 
%
% (C) 2012 Stanford VISTA team. 

% Create output fiber group
fg = dtiNewFiberGroup;
    
% Strip off first and last points from strand
for ll = 1:size(strands,1)
    strand = strands{ll,1}';
    fg.fibers{ll} = strand(:,2:end-1);
end

return