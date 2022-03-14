function angle = solaralt( t,date,hour,location )
%solaralt: solar altitude according to time and location.
%give the t=simulation time, date(month,date), start time (hour) and location(latitude,longitude,hourcircle), ans=[solar altitude,time angle, solar declination].

% Day_month=zeros(1,13);
% Day_number=[0;31;28;31;30;31;30;31;31;30;31;30;31];
% % for i=1:date(1)
% %     Day_month(i)=1;
% % end
% Day_month(1:date(1))=1;
% %Tsaving    = 3600; %daylight saving time.
% Day_start  = Day_month*Day_number+date(2);
Day_number=[31;28;31;30;31;30;31;31;30;31;30;31];
Day_start  = sum(Day_number(1:date(1)-1))+date(2);

dtime      = (location(2)-location(3))/15*60*60;                                         %time difference between local time and clock time                  [second]
day_num    = floor(Day_start+(hour+t/(60*60))/24);                                       %day number
gam        = (day_num-81)*2*pi/365;                                                      %the year angle which is zero at the vernal equinox about March 21  [rad]
dEqofT     = 7.13*cos(gam)+1.84*sin(gam)+0.69*cos(2*gam)-9.92*sin(2*gam);                %A fit curve for equation of time                                   [min]
solar_dec  = 23.45*sin(gam)*pi/180;                                                      %solar declination                                                  [rad]
Angle_time = (hour+(t+dtime)/3600-dEqofT/60-12)*2*pi/24;                                 %time angle                                                         [rad]
solar_alt1 = asin(sin(location(1)*pi/180)*sin(solar_dec)+cos(location(1)*pi/180)*cos(solar_dec).*cos(Angle_time));    %solar altitude                        [rad]
solar_alt  = solar_alt1.*switch01(solar_alt1,-1);                                                                     %switch function of solar altitude

angle      =[solar_alt',Angle_time',solar_dec'];

end

