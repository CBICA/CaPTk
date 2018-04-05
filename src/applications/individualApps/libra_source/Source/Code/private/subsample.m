function tsub = subsample(tvalues, dispflag)
%  Subsample intensity values for faster FCM clustering
%  SINTAX:
%      tsub = subsample(tvalues)
%
%  NOTE: This function assumes that the values in tvalues are z-scored
%
%  Developed by Said Pertuz on Apr 18 2014
%
%  Version info:
%  $Rev: 450 $:     Revision of last commit
%  $Author: hsiehm@UPHS.PENNHEALTH.PRV $:  Author of last commit
%  $Date: 2015-11-20 11:34:20 -0500 (Fri, 20 Nov 2015) $:    Date of last commit
%
%  Contact:
%     CBIG Group <software at cbica.upenn.edu>

B = 5000;   %Bins in histogram
N = 20000;   %Size of sub-sample

if nargin<2, dispflag = false; end

y = linspace(min(tvalues),max(tvalues), B);
h = hist(tvalues, y);
h = h/sum(h);
x = cumsum(h);
xi = linspace(0, 1, N);
[xu, idx] = unique(x);
tsub = interp1(xu, y(idx), xi, 'pchip');

if dispflag
    figure
    y = linspace(min(tvalues),max(tvalues), 100);
    h1 = hist(tvalues, y);
    h2 = hist(tsub, y);
    plot(y, h1/sum(h1), y, h2/sum(h2))
    legend({'Original','Sub-sample'})
end
