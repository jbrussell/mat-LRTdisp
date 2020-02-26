% Program to plot MINEOS model cards
% NJA, 2014
%
% 10/27/2014 Modified to handle either longer period S models (S0_50) or
% shorter period S models (S50_150).
% 
% 10/24/2016 JBR - Reads from .q file instead of .asc file, so that it
% includes the q-corrected period (tq) and phase velocity (phvq)

function [mode] = readMINEOS_qfile_allper(QIN,MODE)

cardname = QIN;

% Read .q file
dat = {};
imode = MODE+1;
com = ['awk ''{ if ($1 ==',num2str(imode-1),' && $10 != "") print $0}'' ',QIN];
[log3, dat{imode}] = system(com);
dat{imode} = str2num(dat{imode});
nn =  dat{imode}(:,1);
ll =  dat{imode}(:,2);
w =   dat{imode}(:,3)/(2*pi)*1000; %convert rad/s ---> mhz
qq =  dat{imode}(:,4);
phi = dat{imode}(:,5);
cv =  dat{imode}(:,6);
gv =  dat{imode}(:,7);
cvq = dat{imode}(:,8);
Tq =  dat{imode}(:,9);
T =   dat{imode}(:,10);

mode.fname=cardname;
mode.n = nn;
mode.l = ll;
mode.wrad = w;
mode.w = w/(2*pi)*1000;
mode.T = T;
mode.Tq = Tq;
mode.grv = gv;
mode.q = qq;
mode.phv = cv;
mode.phvq = cvq;

end