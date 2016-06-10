% % Fall 2013 - Canaveral Swale East (AWAC)
% indexwaveclean=20:3158;
%     % Correct time from EST to UTC (add 4 hours during Daylight Saving)
%     waveclean.t=wave.t(:,indexwaveclean)+4/24;
%     waveclean.tJ=wave.tJ(:,indexwaveclean)+4/24;
%     
% % Fall 2013 - Canaveral Swale West (Aquadopp)
% indexwaveclean=13:1711;
%     % Correct time from EST to UTC (add 4 hours during Daylight Saving)
%     waveclean.t=wave.t(:,indexwaveclean)+4/24;
%     waveclean.tJ=wave.tJ(:,indexwaveclean)+4/24;
% 
% % Spring 2014 - Canaveral Swale East (AWAC)
% indexwaveclean=12:753;
%     waveclean.t=wave.t(:,indexwaveclean);
%     waveclean.tJ=wave.tJ(:,indexwaveclean);
% 
% % Spring 2014 - Canaveral Swale West (Aquadopp)
% indexwaveclean=12:753;
%     waveclean.t=wave.t(:,indexwaveclean);
%     waveclean.tJ=wave.tJ(:,indexwaveclean);
% 
% % Fall 2014 - Canaveral Swale East (AWAC)
% indexwaveclean=11:1843;
%     waveclean.t=wave.t(:,indexwaveclean);
%     waveclean.tJ=wave.tJ(:,indexwaveclean);
% 
% % Fall 2014 - Chester Swale West (Aquadopp)
% indexwaveclean=37:884;
%     waveclean.t=wave.t(:,indexwaveclean);
%     waveclean.tJ=wave.tJ(:,indexwaveclean);
% 
% % Winter 2014-2015 - Chester Swale West (Aquadopp)
% indexwaveclean=63:901;
%     waveclean.t=wave.t(:,indexwaveclean);
%     waveclean.tJ=wave.tJ(:,indexwaveclean);
% 
% % Summer A 2015 - Canaveral Swale East (Aquadopp)
% indexwaveclean=9:532;
%     waveclean.t=wave.t(:,indexwaveclean);
%     waveclean.tJ=wave.tJ(:,indexwaveclean);
% 
% % Summer A 2015 - Chester Swale West (AWAC)
% indexwaveclean=50:530;
%     waveclean.t=wave.t(:,indexwaveclean);
%     waveclean.tJ=wave.tJ(:,indexwaveclean);
% 
% % Fall  B 2015 - Canaveral Swale East (AWAC)
% indexwaveclean=5:676;
%     waveclean.t=wave.t(:,indexwaveclean);
%     waveclean.tJ=wave.tJ(:,indexwaveclean);
% 
% % Fall B 2015 - Chester Swale West (Aquadopp)
% indexwaveclean=2:467;
%     waveclean.t=wave.t(:,indexwaveclean);
%     waveclean.tJ=wave.tJ(:,indexwaveclean);
% 
% % Winter 2015-2016 - Canaveral Swale East (AWAC)
% indexwaveclean=18:1337;
%     waveclean.t=wave.t(:,indexwaveclean);
%     waveclean.tJ=wave.tJ(:,indexwaveclean);
% 
if strcmp(which_instr,'AWAC')
    indexwaveclean=7:1685; % Feb-Apr 2009 - Matanzas (AWAC)
    waveclean.t=wave.t(:,indexwaveclean);
    waveclean.tJ=wave.tJ(:,indexwaveclean);
elseif strcmp(which_instr,'AqDp')
    indexwaveclean=13:1684; % Feb-Apr 2009 - Matanzas (AqDp)
    waveclean.t=wave.t(:,indexwaveclean);
    waveclean.tJ=wave.tJ(:,indexwaveclean);
end