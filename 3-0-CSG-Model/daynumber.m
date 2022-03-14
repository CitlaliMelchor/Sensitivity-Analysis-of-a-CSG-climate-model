function daynumber = daynumber( t,date,hour )
%daynumber: calculate the day number in one year (Feburary has 28 days).
%give the t=simulation time, date(month,date), and start time (hour), ans=day number.

% Day_month=zeros(1,13);
% Day_number=[0;31;28;31;30;31;30;31;31;30;31;30;31];
% for i=1:date(1)
%     Day_month(i)=1;
% end
% %Tsaving    = 3600; %daylight saving time.
% Day_start  = Day_month*Day_number+date(2);

Day_number=[31;28;31;30;31;30;31;31;30;31;30;31];
Day_start  = sum(Day_number(1:date(1)-1))+date(2);

daynumber    = floor(Day_start+(hour+t/(60*60))/24);                                       %day number

n = ceil((daynumber-365)/365);

daynumber = daynumber - n* 365;

end