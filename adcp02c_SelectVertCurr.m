% % Fall 2013 - Canaveral Swale East (AWAC - range 30 m)
% good2=find(z_p0<=10.5);
% % Correct time from EST to UTC (add 4 hours during Daylight Saving)
% t=t+4/24; % in days
% tJ=tJ+4/24; % in days
% 
% % Fall 2013 - Canaveral Swale West (Aquadopp - range 10 m)
% good2=find(z_p0<=10);
% % Correct time from EST to UTC (add 4 hours during Daylight Saving)
% t=t+4/24; % in days
% tJ=tJ+4/24; % in days
% 
% % Spring 2014 - Canaveral Swale East (AWAC - range 30 m)
% good2=find(z_p0<=10);
% 
% % Spring 2014 - Canaveral Swale West (Aquadopp - range 10 m)
% good2=find(z_p0<=10);
% 
% % Fall 2014 - Canaveral Swale East (AWAC - range 30 m)
% good2=find(z_p0<=10.5);
% 
% % Fall 2014 - Chester Swale West (Aquadopp - range 10 m)
% good2=find(z_p0<=10);
% 
% % Winter 2014-2015 - Chester Swale West (Aquadopp - range 10 m)
% good2=find(z_p0<=10);
% 
% % Summer A 2015 - Canaveral Swale East (Aquadopp - range 10 m)
% good2=find(z_p0<=10);
% 
% % Summer A 2015 - Chester Swale West (AWAC - range 30 m)
% good2=find(z_p0<=8);
% 
% % Fall B 2015 - Canaveral Swale West (AWAC - range 30 m)
% good2=find(z_p0<=10.5);
% 
% % Fall B 2015 - Chester Swale West (Aquadopp - range 10 m)
% good2=find(z_p0<=8.5);
% 
% % Winter 2015-2016 - Canaveral Swale West (AWAC - range 30 m)
% good2=find(z_p0<=10);
% 
if strcmp(which_instr,'AWAC')
    good2=find(z_p0<=(min(P)*cosd(25))); % Feb-Apr 2009 - Matanzas (AWAC)
elseif strcmp(which_instr,'AqDp')
    good2=find(z_p0<=(min(P)*cosd(25))); % Feb-Apr 2009 - Matanzas (AqDp)
end
