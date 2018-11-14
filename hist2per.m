function hp = hist2per(hn)
% hist2per	- convert histogram count to percentages
%--------------------------------------------------------------------------------
% Input(s): 	hn - histogram values in counts
% Output(s):	hp - histogram values in % (assuming that sum(hn)=100%)
% Usage:	hp = hist2per(hn)
%
% Last modified 18.11.03
% Copyright (c) 2003 Igor Kagan					 
% kigor@tx.technion.ac.il
% http://igoresha.virtualave.net
%--------------------------------------------------------------------------------

hp = hn/nansum(hn)*100;
