//
//  DvTime.h
//  Identical to QmjdTime, renamed to avoid clash with QDOS
//
//  Dvos - time and time interval object classes
//
//  Tony Allen
//  November 2012
//
	/* MJD measures from the start of 17 Nov 1858 */
	
	/* These utilities use the Gregorian calendar after 4 Oct 1582 (Julian) i.e. from 15 Oct 1582 Gregoriam 
	 Note C libraries use Gregorian only from 14 Sept 1752 
	 More detailed discussion can be found at http://aa.usno.navy.mil/data/docs/JulianDate.php
	 These routines have been compared with the results of the US Naval Observatory online converter.
	 Modified Julian Date (MJD) = Julian Date (JD) - 2400000.5
	 
	 In all routines, specifying a day, hour, minute or second field greater than would be valid is
	 handled with modulo arithmetic and safe. 
	 Thus 2006-12-32 00:62:00.0 will safely, and correctly, be treated as 2007-01-01 01:02:00.0 
	 
	*/
	
#include "DvTime.h"
#include "DvString.h"

#include <stdio.h>  // for sprintf
#include <math.h>  // for pow()
#include <iomanip>  // for cout precision
#include <sstream>


static const double MJDtoJD = 2400000.5;
static const double SecInDay = 86400.; /* we ignore leap seconds */
static const int MonthStartDOY[] = {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334};
static const int MonthStartDOY_L[] = {0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335};

// static variables, common to all instances
static int MJDStartGregorian = -100840;
static int YearStartGregorian = 1582;
static int MonthStartGregorian = 10;
static int DayStartGregorian = 15;
static int DOYStartGregorian = 288;	



static string ISO_FORMAT_D = "%04d-%02d-%02dT%02d:%02d:%02d";
static string ISO_FORMAT_DN = "-%04d-%02d-%02dT%02d:%02d:%02d";
static string ISO_FORMAT_DZ = "%04d-%02d-%02dT%02d:%02d:%02dZ";
static string ISO_FORMAT_DZN = "-%04d-%02d-%02dT%02d:%02d:%02dZ";
static string ISO_FORMAT_S = "%04d-%02d-%02d %02d:%02d:%02d";
static string ISO_FORMAT_SN = "-%04d-%02d-%02d %02d:%02d:%02d";
static int secResolution = DEFAULT_DECIMAL_SEC;
static std::ostringstream secFrac;
	
	

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ class DvTime
// :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

// Constructors
DvTime::DvTime(double sec, int day){
	// default day is Jan 01, 2000 = MJD 51544
	base_day = day;
	time_sec = sec;
}

DvTime::DvTime(const DvTime& ttime){
	base_day = ttime.base_day;
	time_sec = ttime.time_sec;
}

DvTime::DvTime(const char* iso_time)
{
	// from string in ISO time format
	if( !setFromISOstring(iso_time))
	{
		base_day = MJD2000;
		time_sec = 0.0;
	}
    
}                
	
DvTime::DvTime(const DvString &iso_time )
{
    // from string in ISO time format
	if( !setFromISOstring(iso_time))
	{
		base_day = MJD2000;
		time_sec = 0.0;
	}
}                
	
DvTime::DvTime(int year, int month, int day, int hour, int min, double sec)
{
	setFromUT( year, month, day, hour, min, sec); 
}
	

int DvTime::setSecResolution(int resolution)
{
	int oldRes = secResolution;	
	secResolution = resolution;

    secFrac.str("");
    secFrac.clear();
	secFrac << std::fixed;
	secFrac << std::setprecision(secResolution);
	
/*	ISO_FORMAT_D  = "%04d-%02d-%02dT%02d:%02d:%02d";
	ISO_FORMAT_DN = "-%04d-%02d-%02dT%02d:%02d:%02d";
	ISO_FORMAT_DZ  = ISO_FORMAT_D+"Z";
	ISO_FORMAT_DZN = ISO_FORMAT_DN+"Z";
	ISO_FORMAT_S  = "%04d-%02d-%02d %02d:%02d:%02d";
	ISO_FORMAT_SN = "-%04d-%02d-%02d %02d:%02d:%02d";
*/

	return oldRes;
}

int DvTime::getSecResolution(){return secResolution;} // returns previous resolution

	
void DvTime::setFromUT(int year, int month, int day, int hour, int min, double sec)
{	
	int leaps, year4, year100, year400;
	// convert Gregorian or Julian date plus time to MJD 
	// MJD measures from the start of 17 Nov 1858 
	
	// The globals YearStartGregorian etc determine after which date the Gregorian calendar is used 
	// Default is to use Gregorian after 4 Oct 1582 (Julian) i.e. from 15 Oct 1582 Gregorian 
	// Note Unix libraries use Gregorian only from 14 Sept 1752 onwards 
	
	if(month < 1 || month > 12) return;

	if(year <=0) {
		year4 = year - 4;
		year100 = year - 100;
		year400 = year - 400;
	} 
	else {
		year4 = year - 1;
		year100 = year - 1;
		year400 = year - 1;
	}  

	if(year < YearStartGregorian || 
	       (year == YearStartGregorian && month < MonthStartGregorian) || 
		   (year == YearStartGregorian && month == MonthStartGregorian && day < DayStartGregorian) )
	{
		/* count leap years on Julian Calendar */
		/* MJD for Jan 1 0000 (correctly Jan 01, BCE 1) is  - 678943, count from there */
	
		leaps = year4 / 4;
		if(year%4 == 0)
			base_day = MJD_0000J + year * 365 + leaps + MonthStartDOY_L[month-1] + day;
		else
			base_day = MJD_0000J + year * 365 + leaps + MonthStartDOY[month-1] + day;
	}
	else
	{
		/* count leap years Gregorian Calendar - modern dates */
		/* Algorithm below for  17 Nov 1858 (0 MJD) gives */
		/* leaps = 450 and hence base_day of 678941, so subtract it to give MJD day  */
		
		leaps = year4 / 4 - year100 / 100 + year400 / 400;

		if( (year%4 == 0 && year%100 != 0) || (year%4 == 0 && year%400 == 0) )
			base_day = MJD_0000G + year * 365 + leaps + MonthStartDOY_L[month-1] + day;
		else
			base_day = MJD_0000G + year * 365 + leaps + MonthStartDOY[month-1] + day;
	
	}	
		
	time_sec = sec + ( (double) min  +  (double) hour * 60. ) * 60.;

	if(time_sec >= SecInDay)
	{
		int extraDays = (int) (time_sec / SecInDay);
		base_day += extraDays;
		time_sec -= extraDays * SecInDay;
	}
	
	return;
}

bool DvTime::setFromISOstring(const char* ISOstring)
{
	return this->setFromISOstring(DvString(ISOstring));
}

bool DvTime::setFromISOstring(DvString ISOstring)
{
	double seconds;
	int   y, m, d, h, min;
    bool valid = true;
    
	ISOstring.clean(); // remove quotes and spaces outside quotes
	ISOstring.trimmed(); // remove spaces inside quotes
	
	// ISO is "1995-01-23 02:33:17.235" or "1995-01-23T02:33:17.235" or "1995-01-23T02:33:17.235Z" 

	// parse off year 
	DvString token = ISOstring.before('-');
	ISOstring = ISOstring.after('-');
	if(token.empty()) return false;
	
	y = strtol(token.c_str(), NULL, 10);
	
	// parse off month 	
	token = ISOstring.before('-');
	ISOstring = ISOstring.after('-');
	if(token.empty() || ISOstring.length() < 4) return false;
	
	m = strtol(token.c_str(), NULL, 10);
	
	// parse off day 	
	token = ISOstring.substr(0, 2);
	ISOstring = ISOstring.substr(3);
	if(token.empty()) return false;
	
	d = strtol(token.c_str(), NULL, 10);

	// parse off hour 	
	token = ISOstring.front(":");
	if(token.length() != 2) valid = false; // return invalid, but set as best we can
	ISOstring = ISOstring.after(':');
	if(token[0] == 'T') token = token.after('T');
	
	if(token.empty()) h=0;
	else h = strtol(token.c_str(), NULL, 10);

	// parse off minute 	
	token = ISOstring.front(":");
	if(token.length() != 2) valid = false; // return invalid, but set as best we can
	ISOstring = ISOstring.after(':');
	
	if(token.empty()) min=0;
	else min = strtol(token.c_str(), NULL, 10);

	// seconds 	
	if(ISOstring.length() < 2) valid = false; // return invalid, but set as best we can
	if(ISOstring.empty()) seconds=0;
	else seconds = strtod(ISOstring.c_str(), NULL);

	setFromUT(y, m, d, h, min, seconds);
	
	return valid;
}


void DvTime::setFromDOY(int year, int dofy, int hour, int min, double sec)
{	
	/* Set from Day Of Year format */
	
	/* convert Gregorian date plus time to MJD */
	/* MJD measures from the start of 17 Nov 1858 */
	
	/* default is to use Gregorian after 4 Oct 1582 (Julian) i.e. from 15 Oct 1582 Gregorian */
	/* Note C libraries use Gregorian only from 14 Sept 1752 onwards*/

	int leaps, extraDays;
	int year4, year100, year400;
	
	if(year <=0) {
		year4 = year - 4;
		year100 = year - 100;
		year400 = year - 400;
	} 
	else {
		year4 = year - 1;
		year100 = year - 1;
		year400 = year - 1;
	}  

	if(year < YearStartGregorian || (year == YearStartGregorian && dofy < DOYStartGregorian)  )
	{
		/* count leap years on Julian Calendar */
		/* MJD for Jan 1 0000 (correctly Jan 01, BCE 1) is  - 678943, count from there */
	
		leaps = year4 / 4;
		base_day = MJD_0000J + year * 365 + leaps + dofy;
		
	}
	else
	{
		/* count leap years Gregorian Calendar - modern dates */
		/* Algorithm below for  17 Nov 1858 (0 MJD) gives */
		/* leaps = 450 and hence base_day of 678941, so subtract it to give MJD day  */
		
		leaps = year4 / 4 - year100 / 100 + year400 / 400;
		base_day = MJD_0000G + year * 365 + leaps + dofy;
	
	}	
		
	time_sec = sec + ( (double) min  +  (double) hour * 60. ) * 60.;

	if(time_sec >= SecInDay)
	{
		extraDays = (int) (time_sec / SecInDay);
		base_day += extraDays;
		time_sec -= extraDays * SecInDay;
	}
	
	return;

}


void DvTime::setFromBCE(int yearBCE, int month, int day, int hour, int min, double sec)
{
	/* utility to allow user to input dates BCE (BC) */
	
	int year = 1 - yearBCE;
	setFromUT(year, month, day, hour, min, sec);
	
}

void DvTime::setFromMJD(double ModifiedJulianDate)
{
	/* convert MJD double into MJD structure */
	base_day = (int) ModifiedJulianDate;
	time_sec = (ModifiedJulianDate - base_day) * SecInDay;
}

void DvTime::setFromEpoch2000(double sec)
{
    /* Same as default constructor from double */
    
    base_day = MJD2000;
    time_sec = sec;
}

void DvTime::setFromTAI2000(DvInt8 nsec)
{
    // Note DvTime has Msec accuracy, so there is a drift between cdf TT2000 conversions and DvTime
    double sec = nsec * 1.e-9;

    // TAI to UTC corrections
         if(sec >  487980868.184) sec -= 36;
    else if(sec >  393372867.184) sec -= 35;
    else if(sec >  283040066.184) sec -= 34;
    else if(sec >  188345665.184) sec -= 33;
    else if(sec > -32579135.816)  sec -= 32;
    else if(sec > -80012736.816)  sec -= 31;
    else if(sec > -127273537.816) sec -= 30;
    else if(sec > -174707138.816) sec -= 29;
    else if(sec > -206243139.816) sec -= 28;
    else if(sec > -237779140.816) sec -= 27;
    else if(sec > -285039941.816) sec -= 26;
    else if(sec > -316575942.816) sec -= 25;
    else if(sec > -379734343.816) sec -= 24;
    else if(sec > -458703944.816) sec -= 23;
    else if(sec > -521862345.816) sec -= 22;
    else if(sec > -553398346.816) sec -= 21;
    else if(sec > -584934347.816) sec -= 20;
    else if(sec > -632195148.816) sec -= 19;
    else if(sec > -663731149.816) sec -= 18;
    else if(sec > -695267150.816) sec -= 17;
    else if(sec > -726803151.816) sec -= 16;
    else if(sec > -758425552.816) sec -= 15;
    else if(sec > -789961553.816) sec -= 14;
    else if(sec > -821497554.816) sec -= 13;
    else if(sec > -853033555.816) sec -= 12;
    else if(sec > -868931156.816) sec -= 11;
    else if(sec > -884655957.816) sec -= 10;
    else if(sec > -1008207961.629022)  sec -= (sec + 1072958363.501534)  * 0.002592 * 0.0000115   + 4.2131700;
    else if(sec > -1073958363.501534)  sec -= (sec + 1072958363.501534)  * 0.002592 * 0.0000115   + 4.3131700;
    else if(sec > -1084499163.660294)  sec -= (sec + 1104494364.275222)  * 0.001296 * 0.0000115   + 3.8401300;
    else if(sec > -1089855963.840646)  sec -= (sec + 1104494364.275222)  * 0.001296 * 0.0000115   + 3.7401300;
    else if(sec > -1100396764.098758)  sec -= (sec + 1104494364.275222)  * 0.001296 * 0.0000115   + 3.6401300;
    else if(sec > -1105494364.275222)  sec -= (sec + 1104494364.275222)  * 0.001296 * 0.0000115   + 3.5401300;
    else if(sec > -1116035164.533334)  sec -= (sec + 1104494364.275222)  * 0.001296 * 0.0000115   + 3.4401300;
    else if(sec > -1129254364.831622)  sec -= (sec + 1104494364.275222)  * 0.001296 * 0.0000115   + 3.3401300;
    else if(sec > -1137116765.049558)  sec -= (sec + 1104494364.275222)  * 0.001296 * 0.0000115   + 3.2401300;
    else if(sec > -1142387165.1181596) sec -= (sec + 1199188765.9695804) * 0.0011232 * 0.0000115  + 1.9458580;
    else if(sec > -1200188765.9695804) sec -= (sec + 1199188765.9695804) * 0.0011232 * 0.0000115  + 1.8458580;
    else if(sec > -1213407966.167782)  sec -= (sec + 1230724766.392534)  * 0.001296 * 0.0000115   + 1.3728180;
    else if(sec > -1231724766.392534)  sec -= (sec + 1230724766.392534)  * 0.001296 * 0.0000115   + 1.4228180;
    else if(sec > -1263347166.871870)  sec -= (sec + 1230724766.392534)  * 0.001296 * 0.0000115   + 1.4178180;
    // sec (TAI2000) is negative before year 2000
    
    time_sec = sec + TAI2000offset;
    
    base_day = MJD2000;
    

}

void DvTime::setFromTT2000(DvInt8 nsec)
{
    setFromTAI2000(nsec);
    
    time_sec -= 32.184;
}

void DvTime::setFromEpoch1970(double epoch1970)
{
	/* convert MJD double into MJD structure */
	base_day = (int) (epoch1970 / SecInDay) + MJD_1970;
	time_sec = epoch1970 - (base_day - MJD_1970) * SecInDay;
}

void DvTime::setFromJD(double JulianDate)
{
	/* break JD double into MJD based structure
	   Note Julian Day starts Noon, so convert to MJD first */

	base_day = (int) (JulianDate - MJDtoJD) ;
	time_sec = (JulianDate - MJDtoJD - (double) base_day) * SecInDay; 
}

void DvTime::setFromCDFepoch(double cdfepoch){

	/* convert cdf epoch double into MJD structure 
	   Note that cdfepoch is msec from 0 AD on the Gregorian calendar */
	   
	double seconds = cdfepoch * 0.001;

	base_day = (int) (seconds / 86400.0);
	time_sec = seconds - base_day * SecInDay;
	base_day += MJD_0000G;

}

void DvTime::setFromCDFepoch16(double *cdfepoch){

	/* convert cdf epoch double into MJD structure 
	   Note that cdfepoch is msec from 0 AD on the Gregorian calendar */
	   
	double seconds = cdfepoch[0];

	base_day = (int) (seconds / 86400.0);
	time_sec = seconds - base_day * SecInDay;
	base_day += MJD_0000G;
	
	time_sec += cdfepoch[1] * 1.e-12;

}

double DvTime::getCDFepoch(){

	/* convert MJD structure into cdf epoch double 
	   Note that cdfepoch is msec from 0 AD on the Gregorian Calendar */
	   
	int days = base_day - MJD_0000G;
	double seconds = days * SecInDay + time_sec;
	return seconds * 1000.;
}

void DvTime::getCDFepoch16(double *epoch16){

	/* convert MJD structure into cdf epoch16 doubles
	   Note that cdfepoch16 is picosec from 0 AD on the Gregorian Calendar */
	   
	int days = base_day - MJD_0000G;
	
	int int_sec = (int) time_sec;
	
	epoch16[0] = days * SecInDay + (double) int_sec;
	
	epoch16[1] = (time_sec - int_sec) * 1.e12;
	
	return;
}

double DvTime::getEpoch2000(){

	// convert MJD structure into old qdos double double 
	   
	int days = base_day - MJD2000;
	double seconds = days * SecInDay + time_sec;
	return seconds;
}

double DvTime::getEpoch1970(){

	// convert MJD structure into old qdos double double 
	   
	int days = base_day - MJD_1970;
	double seconds = days * SecInDay + time_sec;
	return seconds;
}


const char * DvTime::getDayOfWeek()
{
	static const char *dow = {"Wed\0Thu\0Fri\0Sat\0Sun\0Mon\0Tue"};
	int d = base_day % 7;
	return &(dow[d*4]);
}

const char * DvTime::getLongDayOfWeek()
{
	static const char *dow = {"Wednesday\0Thursday\0\0Friday\0\0\0\0Saturday\0\0Sunday\0\0\0\0Monday\0\0\0\0Tuesday"};
	int d = base_day % 7;
	return &(dow[d*10]);
}

const char * DvTime::getMonth( int m)
{
	static const char *months = {"Jan\0Feb\0Mar\0Apr\0May\0Jun\0Jul\0Aug\0Sep\0Oct\0Nov\0Dec"};
	return &(months[(m-1)*4]);
}

const char * DvTime::getLongMonth( int m)
{
	static const char *months = {"January\0\0\0February\0\0March\0\0\0\0\0April\0\0\0\0\0May\0\0\0\0\0\0\0June\0\0\0\0\0\0July\0\0\0\0\0\0August\0\0\0\0September\0October\0\0\0November\0\0December"};
	return &(months[(m-1)*10]);
}

DvTime DvTime::normalized() const
{
  int extra_days;
  
  // Calculate normalized version (i.e., 0. <= MJDout->time_sec < 86400.) 
  if(time_sec >= 0) {
    extra_days = (int) (time_sec / SecInDay);
  } 
  else {
    /* allow for negative seconds push into previous day even if less than 1 day */
    extra_days = (int) (time_sec / SecInDay) - 1 ;
  }
  
  return DvTime(time_sec - extra_days * SecInDay, base_day + extra_days);
}

bool DvTime::getYandD(int *year, int *dofy)
{
	int j, year4, year100, year400;
	bool incorrect=false;
	bool isLeapyear;
	
	j = this->base_day;

	if( j < MJDStartGregorian)
	{
		/* Julian Dates and Julian proleptic calendar */
		/* Shift j epoch to 0000-01-01 for the Julian proleptic calendar.*/
		 j -= MJD_0000J;

		/* 365.25 is the exact period of the Julian year so year will be correct
		if the day offset is set exactly right so that years -4, 0, 4 are
		leap years, i.e. years -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5 start with
		j =  -1826 -1461, -1095, -730, -365, 0, 366, 731, 1096, 1461, 1827 */
				 
		if(j >= 366) {
			*year = (int) ((double)(j) / 365.25);
			year4 = *year-1;
		} else 
		{
			*year = (int) ((double)(j-365)/ 365.25);
			year4 = *year-4;
		}
    
		*dofy = j - *year * 365 - year4 / 4;
		isLeapyear = (*year%4 == 0);				
	}
	else 
	{
		/* Gregorian Dates and  Gregorian proleptic calendar */
		/* Shift j epoch to 0000-01-01 for the Gregorian proleptic calendar.*/
		j -= MJD_0000G;
				 
		if(j >=366) {
			*year	= (int) ((double)(j) / 365.2425);
			year4 = *year - 1;
			year100 = *year - 1;
			year400 = *year - 1;
		} else {
			*year = (int) ((double)(j-365) / 365.2425);
			year4 = *year - 4;
			year100 = *year - 100;
			year400 = *year - 400;
		}
		 
		*dofy = j - *year * 365 - year4 / 4 + year100 / 100 - year400 / 400;
		isLeapyear = (*year%4 == 0 && *year%100 != 0) || (*year%4 == 0 && *year%400 == 0);		 
		
		
		/* Rare corrections to above average Gregorian relations. */
		if(*dofy < 1) {
			(*year)--;
			incorrect = true;
		} 
		else if(*dofy > 365 && (!isLeapyear || *dofy > 366)) {
			(*year)++;
			incorrect = true;
		} 
		
		if(incorrect) {
			if(j >=366) {
				year4 = *year - 1;
				year100 = *year - 1;
				year400 = *year - 1;
			} 
			else {
				year4 = *year - 4;
				year100 = *year - 100;
				year400 = *year - 400;
			}

			*dofy = j - *year * 365 - year4 / 4 + year100 / 100 - year400 / 400;
			isLeapyear = (*year%4 == 0 && *year%100 != 0) || (*year%4 == 0 && *year%400 == 0);
		}
	}
	return isLeapyear;
}


int DvTime::breakDownMJD(int *year, int *month, int *day, int *hour, int *min, double *sec) const
{ 	
	/* Convert MJD struct into date/time elements */
	/* Note year 0 CE (AD) [1 BCE (BC)] is a leap year */
	/* There are 678943 days from year 0 to MJD(0)   */
	
	int dofy;
	
	DvTime nMJD = this->normalized();

	/* Time part */	
	*sec = nMJD.time_sec;
	if(*sec < -1.e30) {
		// epoch 16 FILLVAL must be trapped 
		*sec = 59.999999999; 
		*min = 59;
		*hour = 23;
		*day = 31;
		*month = 12;
		*year =9999; 
		return 365;
	}
	

	*hour = (int)( *sec / 3600.);
	*sec -= (double) *hour * 3600.;
	*min = (int) ( *sec / 60.);
	*sec -=  (double) *min * 60.;
	
	bool isLeapyear = nMJD.getYandD(year, &dofy);
		
	/* turn day of year into month and day */
	*month = 0;
	if(isLeapyear)
	{
		while(dofy > MonthStartDOY_L[*month])  
		{
			(*month)++;
			if(*month == 12) break;
		}
		*day = dofy - MonthStartDOY_L[*month -1];
	}
	else
	{
		while(dofy > MonthStartDOY[*month]) 
		{
			(*month)++;
			if(*month == 12) break;
		}
		*day = dofy - MonthStartDOY[*month -1];
	}
	
	return dofy;

}

double DvTime::getMJD()
{
	/* Return MJD as a double */
	return  (double) base_day + time_sec / SecInDay ;
}

double DvTime::getJD()
{
	/* Return JD as a double */
	double JD = getMJD() + MJDtoJD; 
	return JD;
}

double DvTime::getDiffDays(const DvTime& MJD) const
{
	/* Return difference this - MJD in days as a double */
	double diff = (double)(this->base_day - MJD.base_day) + (this->time_sec - MJD.time_sec) / SecInDay;
	return diff;
}

double DvTime::getDiffSecs(const DvTime& MJD) const
{
	/* Return difference this - MJD in seconds as a double */
	double diff = (double)(this->base_day - MJD.base_day) * SecInDay + (this->time_sec - MJD.time_sec) ;
	return diff;
}

std::string DvTime::getTAIstring(int delim) const
{
    return string("");
}

std::string DvTime::getTTstring(int delim) const
{
    return string("");
}

std::string DvTime::getISOstring(int delim) const
{
	// ISO time string for UTC 
	char DateTime[50];
	int y, m, d, hour, min;
	int secInt, ysign;
	double sec;
	int slen;

	breakDownMJD(&y, &m, &d, &hour, &min, &sec);

	if(y < 0)
	{
		ysign = 1;
		y=-y;
	}
	else ysign = 0;	
	
	secInt = (int)sec;
	sec -= (double) secInt;
	if(sec < 1.e-10) sec = 0.0;

    secFrac.str("");
    secFrac.clear();
    secFrac.seekp(0);
	secFrac << sec;	
	
	if(delim == 2)
	{
		// T & Z delimiters
		if(ysign == 0)
			sprintf(DateTime, ISO_FORMAT_DZ.c_str(), y, m, d, hour, min, secInt );
		else
			sprintf(DateTime, ISO_FORMAT_DZN.c_str(), y, m, d, hour, min, secInt );
			
		/* remove trailing white space */
		char * ptr;
		while( ( ptr = strrchr(&(DateTime[0]), ' ')) != NULL)  	ptr[0] ='\0';
	}
	else if(delim == 1)
	{
		// T deliiter
		if(ysign == 0)
			sprintf(DateTime, ISO_FORMAT_D.c_str(), y, m, d, hour, min, secInt );
		else
			sprintf(DateTime, ISO_FORMAT_DN.c_str(), y, m, d, hour, min, secInt );
			
		/* remove trailing white space */
		char * ptr;
		while( ( ptr = strrchr(&(DateTime[0]), ' ')) != NULL)  	ptr[0] ='\0';
	}
	else
	{
		// space delimited
		if(ysign == 0)
			sprintf(DateTime, ISO_FORMAT_S.c_str(), y, m, d, hour, min, secInt );
		else
			sprintf(DateTime, ISO_FORMAT_SN.c_str(), y, m, d, hour, min, secInt );
		
		/* remove trailing white space */
		slen = strlen(DateTime)-1;
		while( DateTime[slen] == ' ') 
		{
			DateTime[slen] ='\0';
			slen--;
		}
	}
	std::string res = secFrac.str();
	res.replace(0, 1, DateTime);
	
	return res;
}



size_t DvTime::strfMJD(char * buf, size_t len, const char *format)
{
	/* Format a text string according to the format string.
	   Uses the same syntax as strftime() but does not use current locale. 
	   The null terminator is included in len for safety. */
	
	int year, month, day, hour, min, ysign, sec1, second;
	size_t i;
	int nplaces, slen, count;
	int resolution;
	size_t format_len=strlen(format);
	double shiftPlaces;
	char * ptr;
	double sec;
	const char *dayText;
	const char *monthText;
	char DateTime[80];
	size_t posn = 0;
	size_t last = len -1;
	buf[last] = '\0';
	buf[0] = '\0'; /* force overwrite of old buffer since strncat() used hereafter */
	
	/* Find required resolution */
	resolution = 0;
	i=0;
	while(i<format_len)
	{
		char next = format[i];
		if( next == '%')
		{
			/* find seconds format if used */
			i++;
			next = format[i];
			if( isdigit(next) != 0 )
			{
				nplaces = strtol(&(format[i]), NULL, 10); 
				if(nplaces > resolution) resolution = nplaces;
			}	
			else if( next == '.' )
			{
				resolution = 9; /* maximum resolution allowed */
			}
		}
		i++;
	}
	
	shiftPlaces = std::pow(10.,(double)resolution);
	
	DvTime nMJD = DvTime(*this);
	nMJD.time_sec += 0.5/shiftPlaces;

	int dofy = nMJD.breakDownMJD(&year, &month, &day, &hour, &min, &sec);
	if(year < 0)
	{
		ysign = 1;
		year = -year;
	}
	else ysign = 0;	
	
	second = (int) sec;
	sec1 = (int)sec/10;
	sec -= (double) sec1*10;
	
	/* Read format string, character at a time */
	i=0;
	while(i<strlen(format))
	{
		char next = format[i];
		if( next == '%')
		{
			/* format character or escape */
			i++;
			next =  format[i];
			if(next == '%')
			{
				/* escaped %, pass it on */
				buf[posn] = next;
				posn++;
				if(posn >= last) return posn;
			}
			else if(next == 'a')
			{
				/* short day name */
				dayText = getDayOfWeek();
				strncat(&(buf[posn]), dayText, last - posn);
				posn = strlen(buf);
				if(posn >= last) return posn;
			}
			else if(next == 'A')
			{
				/* long day name */
				dayText = getLongDayOfWeek();
				strncat(&(buf[posn]), dayText, last - posn);
				posn = strlen(buf);
				if(posn >= last) return posn;
			}
			else if(next == 'b' || next == 'h')
			{
				/* short month name */
				monthText = getMonth(month);
				strncat(&(buf[posn]), monthText, last - posn);
				posn = strlen(buf);
				if(posn >= last) return posn;
			}
			else if(next == 'B')
			{
				/* long month name */
				monthText = getLongMonth(month);
				strncat(&(buf[posn]), monthText, last - posn);
				posn = strlen(buf);
				if(posn >= last) return posn;
			}
			else if(next == 'c')
			{
				/* Date and Time with day of week */
				dayText = getDayOfWeek();
				monthText = getMonth(month);
				if(ysign == 0)
					sprintf(DateTime, "%s %s %02d %02d:%02d:%02d %04d", dayText, monthText, day, hour, min, second, year );
				else
					sprintf(DateTime, "%s %s %02d %02d:%02d:%02d -%04d", dayText, monthText, day, hour, min, second, year );
					
				strncat(&(buf[posn]), DateTime, last - posn);
				posn = strlen(buf);
				if(posn >= last) return posn;
			}
			else if(next == 'C')
			{
				/*  year / 100 so, e.g. 1989 is 20th century but comes out as 19 */
				int century = year / 100;
				if(ysign == 0)
					sprintf(DateTime, "%02d", century );
				else
					sprintf(DateTime, "-%02d", century+1 );
					
				strncat(&(buf[posn]), DateTime, last - posn);
				posn = strlen(buf);
				if(posn >= last) return posn;
			}
			else if(next == 'd')
			{
				/* day of month (01 - 31) */
				sprintf(DateTime, "%02d", day);
					
				strncat(&(buf[posn]), DateTime, last - posn);
				posn = strlen(buf);
				if(posn >= last) return posn;
			}
			else if(next == 'D')
			{
				/* month/day/year */
				int y = year %100;
				if(ysign == 0)
					sprintf(DateTime, "%02d/%02d/%02d", month, day, y );
				else
					sprintf(DateTime, "%02d/%02d/-%02d", month, day, y );
					
				strncat(&(buf[posn]), DateTime, last - posn);
				posn = strlen(buf);
				if(posn >= last) return posn;
			}
			else if(next == 'e')
			{
				/* day of month ( 1 - 31) */
				if(day < 10)
					sprintf(DateTime, " %01d", day);
				else
					sprintf(DateTime, "%02d", day);
					
				strncat(&(buf[posn]), DateTime, last - posn);
				posn = strlen(buf);
				if(posn >= last) return posn;
			}
			else if(next == 'F')
			{
				/* year-month-day */
				if(ysign == 0)
					sprintf(DateTime, "%04d-%02d-%02d", year, month, day );
				else
					sprintf(DateTime, "-%04d-%02d-%02d", year, month, day );
					
				strncat(&(buf[posn]), DateTime, last - posn);
				posn = strlen(buf);
				if(posn >= last) return posn;
			}
			else if(next == 'H')
			{
				/* hour, 24 hour clock (00 - 23) */
				sprintf(DateTime, "%02d", hour);
					
				strncat(&(buf[posn]), DateTime, last - posn);
				posn = strlen(buf);
				if(posn >= last) return posn;
			}
			else if(next == 'I')
			{
				/* hour, 12 hour clock (01 - 12) */
				if(hour == 0) sprintf(DateTime, "%02d", hour+12);
				else if(hour > 12) 	 sprintf(DateTime, "%02d", hour-12);
				else  sprintf(DateTime, "%02d", hour);
				
				strncat(&(buf[posn]), DateTime, last - posn);
				posn = strlen(buf);
				if(posn >= last) return posn;
			}
			else if(next == 'j')
			{
				/* day of year */
				sprintf(DateTime, "%03d", dofy);
				
				strncat(&(buf[posn]), DateTime, last - posn);
				posn = strlen(buf);
				if(posn >= last) return posn;
			}
			else if(next == 'k')
			{
				/* hour, 24 hour clock ( 0 - 23) */
				if(hour < 10)
					sprintf(DateTime, " %01d", hour);
				else
					sprintf(DateTime, "%02d", hour);
					
				strncat(&(buf[posn]), DateTime, last - posn);
				posn = strlen(buf);
				if(posn >= last) return posn;
			}
			else if(next == 'l')
			{
				/* hour, 12 hour clock ( 1 - 12) */
				if(hour == 0) sprintf(DateTime, "%02d", hour+12);
				else if(hour < 10) sprintf(DateTime, " %01d", hour);
				else if(hour <= 12) sprintf(DateTime, "%02d", hour);
				else if(hour < 22)  sprintf(DateTime, " %01d", hour-12);
				else sprintf(DateTime, "%02d", hour-12);
				
				strncat(&(buf[posn]), DateTime, last - posn);
				posn = strlen(buf);
				if(posn >= last) return posn;
			}
			else if(next == 'm')
			{
				/* month (01 - 12) */
				sprintf(DateTime, "%02d", month);
				
				strncat(&(buf[posn]), DateTime, last - posn);
				posn = strlen(buf);
				if(posn >= last) return posn;
			}
			else if(next == 'M')
			{
				/* minute (00 - 59) */
				sprintf(DateTime, "%02d", min);
				
				strncat(&(buf[posn]), DateTime, last - posn);
				posn = strlen(buf);
				if(posn >= last) return posn;
			}
			else if(next == 'n')
			{
				/*  newline */
				buf[posn] = '\n';
				posn++;
				if(posn >= last) return posn;
			}
			else if(next == 'p')
			{
				/* am/pm on12 hour clock  */
				if(hour < 0) sprintf(DateTime, "AM");
				else  sprintf(DateTime, "PM");
				
				strncat(&(buf[posn]), DateTime, last - posn);
				posn = strlen(buf);
				if(posn >= last) return posn;
			}
			else if(next == 'r')
			{
					/* hour:min:sec AM, 12 hour clock (01 - 12):(00 - 59):(00 - 59) (AM - PM) */
				if(hour == 0) sprintf(DateTime, "%02d:%02d:%02d AM", hour+12, min, second);
				else if(hour > 12) 	 sprintf(DateTime, "%02d:%02d:%02d PM", hour-12, min, second);
				else if(hour == 12)  sprintf(DateTime, "%02d:%02d:%02d PM", hour, min, second);
				else  sprintf(DateTime, "%02d:%02d:%02d AM", hour, min, second);
				
				strncat(&(buf[posn]), DateTime, last - posn);
				posn = strlen(buf);
				if(posn >= last) return posn;
			}
			else if(next == 'R')
			{
				/* hour:min, 24 hour clock (00 - 23):(00 - 59) */
				sprintf(DateTime, "%02d:%02d", hour, min);
					
				strncat(&(buf[posn]), DateTime, last - posn);
				posn = strlen(buf);
				if(posn >= last) return posn;
			}
			else if(next == 'S')
			{
				/* second (00 - 59) */
				sprintf(DateTime, "%02d", second);
				
				strncat(&(buf[posn]), DateTime, last - posn);
				posn = strlen(buf);
				if(posn >= last) return posn;
			}
			else if(next == 's')
			{
				/* seconds since 01 Jan 1970 Gregorian */
				int s = (int) (time_sec + (base_day - MJD_1970) * SecInDay);
				sprintf(DateTime, "%d", s);
				
				strncat(&(buf[posn]), DateTime, last - posn);
				posn = strlen(buf);
				if(posn >= last) return posn;
			}
			else if(next == 't')
			{
				/*  tab */
				buf[posn] = '\t';
				posn++;
				if(posn >= last) return posn;
			}
			else if(next == 'T')
			{
				/* hour:min:sec, 24 hour clock (00 - 23):(00 - 59):(00 - 59) */
				sprintf(DateTime, "%02d:%02d:%02d", hour, min, second);
					
				strncat(&(buf[posn]), DateTime, last - posn);
				posn = strlen(buf);
				if(posn >= last) return posn;
			}
			else if(next == 'U')
			{
				/* week of year as a number,  (00 - 53) start of week is Sunday */
				int w;
				int days_in_wk1 = (this->base_day - dofy - 4) % 7; 
				
				w = (dofy + 6 - days_in_wk1) / 7;
				
				sprintf(DateTime, "%02d", w);
				
				strncat(&(buf[posn]), DateTime, last - posn);
				posn = strlen(buf);
				if(posn >= last) return posn;
			}
			else if(next == 'u')
			{
				/* weekday as a number,  0 = Monday */
				int d = 1 + (this->base_day - 5) % 7;

				sprintf(DateTime, "%01d", d);
				
				strncat(&(buf[posn]), DateTime, last - posn);
				posn = strlen(buf);
				if(posn >= last) return posn;
			}
			else if(next == 'v')
			{
				/* day-MonthName-year day of month ( 1 - 31) - month (Jan ... Dec) - year (yyyy) */
				
				monthText = getMonth(month);
				
				if(ysign == 0)
				{
					if(day < 10)
						sprintf(DateTime, " %01d-%s-%04d", day, monthText, year);
					else
						sprintf(DateTime, "%02d-%s-%04d", day, monthText, year);
				}
				else
				{
					if(day < 10)
						sprintf(DateTime, " %01d-%s-(-)%04d", day, monthText, year);
					else
						sprintf(DateTime, "%02d-%s-(-)%04d", day, monthText, year);
				}
				
				strncat(&(buf[posn]), DateTime, last - posn);
				posn = strlen(buf);
				if(posn >= last) return posn;
			}
			else if(next == 'V')
			{
				/* week of year as a number,  (01 - 53) start of week is Monday and first week has at least 3 days in year */
				int w;
				int days_in_wk1 = (this->base_day - dofy - 3) % 7; 
				
				if(days_in_wk1 <= 3) w = (dofy +6 - days_in_wk1) / 7; /* ensure first week has at least 3 days in this year */
				else w = 1 + (dofy + 6 - days_in_wk1) / 7;
				
				if(w == 0) w = 53;
				sprintf(DateTime, "%02d", w);
				
				strncat(&(buf[posn]), DateTime, last - posn);
				posn = strlen(buf);
				if(posn >= last) return posn;
			}
			else if(next == 'w')
			{
				/* weekday as a number,  0 = Sunday */
				int d = (this->base_day - 4) % 7;

				sprintf(DateTime, "%01d", d);
				
				strncat(&(buf[posn]), DateTime, last - posn);
				posn = strlen(buf);
				if(posn >= last) return posn;
			}
			else if(next == 'W')
			{
				/* week of year as a number,  (00 - 53) start of week is Monday */
				int w;
				int days_in_wk1 = (this->base_day - dofy - 3) % 7; 
				
				w =  (dofy +6 - days_in_wk1) / 7;
				
				sprintf(DateTime, "%02d", w);
				
				strncat(&(buf[posn]), DateTime, last - posn);
				posn = strlen(buf);
				if(posn >= last) return posn;
			}
			else if(next == 'x')
			{
				/* date string */
				dayText = getDayOfWeek();
				monthText = getMonth(month);
				if(ysign == 0)
					sprintf(DateTime, "%s %s %02d, %04d", dayText, monthText, day, year );
				else
					sprintf(DateTime, "%s %s %02d, -%04d", dayText, monthText, day, year );
					
				strncat(&(buf[posn]), DateTime, last - posn);
				posn = strlen(buf);
				if(posn >= last) return posn;
			}
			else if(next == 'X')
			{
				/* time string */
				sprintf(DateTime, "%02d:%02d:%02d", hour, min, second );
						
				strncat(&(buf[posn]), DateTime, last - posn);
				posn = strlen(buf);
				if(posn >= last) return posn;
			}
			else if(next == 'y')
			{
				/* 2 digit year */
				int y = year %100;
				
				if(ysign == 0)
					sprintf(DateTime, "%02d", y );
				else
					sprintf(DateTime, "-%02d", y );
				
				strncat(&(buf[posn]), DateTime, last - posn);
				posn = strlen(buf);
				if(posn >= last) return posn;
			}
			else if(next == 'Y')
			{
				/* 4 digit year */
				if(ysign == 0)
					sprintf(DateTime, "%04d", year );
				else
					sprintf(DateTime, "-%04d", year );
				
				strncat(&(buf[posn]), DateTime, last - posn);
				posn = strlen(buf);
				if(posn >= last) return posn;
			}
			else if(next == 'Z')
			{
				/* time zone and calendar, alwaus UTC */
				if(year < YearStartGregorian || 
				(year == YearStartGregorian && month < MonthStartGregorian) || 
				(year == YearStartGregorian && month == MonthStartGregorian && day < DayStartGregorian) )
					strncat(&(buf[posn]), "UTC Julian", last - posn);
				else
					strncat(&(buf[posn]), "UTC Gregorian", last - posn);
					
				posn = strlen(buf);
				if(posn >= last) return posn;
			}
			else if(next == 'z')
			{
				/* time zone, always UTC */
				strncat(&(buf[posn]), "+0000", last - posn);
				posn = strlen(buf);
				if(posn >= last) return posn;
			}
			else if(next == '+')
			{
				/* date and time */
				dayText = getDayOfWeek();
				monthText = getMonth(month);
				if(ysign == 0)
					sprintf(DateTime, "%s %s %02d %02d:%02d:%02d UTC %04d",  dayText, monthText, day, hour, min, second, year );
				else
					sprintf(DateTime, "%s %s %02d %02d:%02d:%02d UTC -%04d", dayText, monthText, day, hour, min, second, year );
					
				strncat(&(buf[posn]), DateTime, last - posn);
				posn = strlen(buf);
				if(posn >= last) return posn;
			}
			else if( isdigit(next) != 0 )
			{
				nplaces = strtol(&(format[i]), NULL, 10); 
				/* numeric value is number of decimal places ( > 0 ) */
				double sec_fraction = sec - (double) second;

				for(count=0; count<nplaces; count++) sec_fraction *= 10;
				sprintf(DateTime, ".%d",  (int) sec_fraction);
				
				/* append 0 to pad to length */
				slen = strlen(DateTime);
				while(slen < nplaces+1)
				{
					DateTime[slen] = '0';
					slen++;
					DateTime[slen] = '\0';
				}
				strncat(&(buf[posn]), DateTime, last - posn);
				posn = strlen(buf);
				if(posn >= last) return posn;
			}	
			else if( next == '.' )
			{
				/* fractional part of seconds to maximum available accuracy */
				double sec_fraction = sec - (double) second;
				sprintf(DateTime, "%-11.9f",  sec_fraction);
				while( ( ptr = strrchr(&(DateTime[0]), ' ')) != NULL)  ptr[0] ='\0'; /* remove trailing white space */
				slen = strlen(DateTime) -1;
				while( DateTime[slen] == '0' && DateTime[slen-1] != '.') {
					DateTime[slen] ='\0'; /* remove trailing zeros */
					slen --;
				}
				
				ptr = strchr(DateTime, '.'); /* remove possible lead 0 */
				strncat(&(buf[posn]), ptr, last - posn);
				posn = strlen(buf);
				if(posn >= last) return posn;
			}
		}
		else
		{
			/* regular multi-byte character, pass it on */
			buf[posn] = next;
			posn++;
			if(posn >= last) return posn;
		}
		buf[posn] = '\0';
		i++;
	}
	return posn;
}

void DvTime::setGregorianStartMJD(double GregorianMJD)
{
	int y, dofy, m, d, h, min;
	double s;
	DvTime mjd;
	MJDStartGregorian = (int) GregorianMJD; /* used in breakDownMJD, truncated to start of day */
	mjd.setFromMJD(GregorianMJD);
	dofy = mjd.breakDownMJD(&y, &m, &d, &h, &min, &s);
	YearStartGregorian = y;
	MonthStartGregorian = m;
	DayStartGregorian = d;
	DOYStartGregorian = dofy;

}


// +++++++++++++++++++++++++++++++++++++++++++++++++ class DvEvent
// 
    
DvEvent::DvEvent( const DvTime& st, const DvTime& ed, Sense option )
{
	set(st,ed,option);
}

DvEvent::DvEvent( string &ISOstring ){ 
	if( !this->setFromISOstring(ISOstring)){
		_t1.base_day = MJD2000;
		_t1.time_sec = 0.0;
		_t2.base_day = MJD2000;
		_t2.time_sec = 1.0;
	}
}
      


// --------------------------------------------------------- contains
bool DvEvent::contains( const DvTime& t ) const
{
  bool res;
  Sense sense=get_sense();
  if(sense==TIME_SENSE_NEGATIVE)
    res=(t<=_t1 && _t2<=t);
  else
    res=(_t1<=t && t<=_t2);
  return res;
}

// --------------------------------------------------------- intersects
bool DvEvent::intersects( const DvEvent& tivl ) const
{
	if(this->contains(tivl.start()) ) return true;
	if(this->contains(tivl.end()) ) return true;
	if(tivl.contains(this->start()) ) return true;
	if(tivl.contains(this->end()) ) return true;
	return false;
}

// -------------------------------------------------------- event_union
DvEvent DvEvent::event_union( const DvEvent& tint ) const
{
    if (empty())
        return tint;
    if(tint.empty())
        return *this;    
	DvTime t1 = start_abs();
	DvTime t2 = end_abs();
	if(t1 > tint.start_abs()) t1 = tint.start_abs();
	if(t2 < tint.end_abs()) t2 = tint.end_abs();
	 
    return DvEvent(t1, t2);
}


// ---------------------------------------------------- event_intersect
DvEvent DvEvent::event_intersect( const DvEvent& tint ) const
{
    DvEvent res; // ensure empty interval constructed MJD2000 to MJD2000
    if((tint.start_abs()<end_abs()) && (tint.end_abs() >start_abs()))
	{
		DvTime t1 = start_abs();
		DvTime t2 = end_abs();
		if(t1 < tint.start_abs()) t1 = tint.start_abs();
		if(t2 > tint.end_abs()) t2 = tint.end_abs();

        res.set(t1, t2);
	}
    return res;
}

std::string DvEvent::getISOstring(int delim) const
{
	// ISO time string  
	string res = start().getISOstring(delim);
	res += "/";
	res += end().getISOstring(delim);
	
	return res;
}


