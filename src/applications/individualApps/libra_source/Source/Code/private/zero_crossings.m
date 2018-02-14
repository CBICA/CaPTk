function [num_xings]=zero_crossings(t,bw,gw,sig)
%  compute zero-crossings code
%
%  Version info:
%  $Rev: 493 $:     Revision of last commit
%  $Author: hsiehm@UPHS.PENNHEALTH.PRV $:  Author of last commit
%  $Date: 2015-12-23 10:42:16 -0500 (Wed, 23 Dec 2015) $:    Date of last commit
%
%  Contact:
%     CBIG Group <software at cbica.upenn.edu>

[h, ~]=hist(t,-4:bw:4);

%make smoothing window - heavy filtering, we dont care about very small
%variations for this
w=gw;
g=gausswin(w+1,sig);

%smooth hist and take 1st order derivative: zero-crossing=local min/max
hc=conv(h(1:end),g,'same');
d=diff(hc);

%Find (+) -> (-) zero crossings
e=d>0;
z=(e(1:end-1)-e(2:end))==1; 

%How many are there that are not just noise crossings
%...1)What are the potential crossings
%...2)Require significant second derivative to not be a noise crossing
d2=diff(d);%figure(3);plot(d2);
amp=max(abs(d2(z)));
min_d2=amp*0.001;%*0.001;

num_xings=sum(abs(d2(z))>=min_d2);
