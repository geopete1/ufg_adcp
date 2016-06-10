% % First and last 'good' data (MUST BE THE SAME AS ABOVE IN 'wave')
% % Fall 2013 - Canaveral Swale East (AWAC)
% indexgoodwp=20:3158;
%     % Correct time from EST to UTC (add 4 hours during Daylight Saving)
%     h_corr=4;
%     
% % Fall 2013 - Canaveral Swale West (Aquadopp)
% indexgoodwp=13:1711;
%     % Correct time from EST to UTC (add 4 hours during Daylight Saving)
%     h_corr=4;
% 
% % Spring 2014 - Canaveral Swale East (AWAC)
% indexgoodwp=12:753;
%     h_corr=0;
% 
% % Spring 2014 - Canaveral Swale West (Aquadopp)
% indexgoodwp=12:753;
%     h_corr=0;
% 
% % Fall 2014 - Canaveral Swale East (AWAC)
% indexgoodwp=11:1843;
%     h_corr=0;
% 
% % Fall 2014 - Chester Swale West (Aquadopp)
% indexgoodwp=37:884;
%     h_corr=0;
% 
% % Winter 2014-2015 - Chester Swale West (Aquadopp)
% indexgoodwp=63:901;
%     h_corr=0;
% 
% % Summer A 2015 - Canaveral Swale East (Aquadopp)
% indexgoodwp=9:532;
%     h_corr=0;
% 
% % Summer A 2015 - Chester Swale West (AWAC)
% indexgoodwp=50:530;
%     h_corr=0;
% 
% % Fall B 2015 - Canaveral Swale East (AWAC)
% indexgoodwp=6:676;
%     h_corr=0;
% 
% % Fall B 2015 - Chester Swale West (Aquadopp)
% indexgoodwp=2:467;
%     h_corr=0;
% 
% % Winter 2015-2016 - Canaveral Swale East (AWAC)
% indexgoodwp=18:1337;
%     h_corr=0;
% 
if strcmp(which_instr,'AWAC')
    indexgoodwp=7:1685; % Feb-Apr 2009 - Matanzas (AWAC)
    h_corr=0;
elseif strcmp(which_instr,'AqDp')
    indexgoodwp=13:1684; % Feb-Apr 2009 - Matanzas (AqDp)
    h_corr=0;
end
