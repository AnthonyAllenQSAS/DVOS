//
//  DvTime.h
//
//  Dvos - time and time interval object classes
//
//  Tony Allen
//  November 2012
//  January 2015 - added TT2000
//  2016 - leap second handling

#ifndef DV_TIME_H
#define DV_TIME_H

#include <string.h>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <time.h>


//#ifdef HAVE_SSTREAM_H
//#include <sstream>
//#else
//#include <strstream>
//#endif

using namespace std;

class DvString;

typedef long long DvInt8;

// MJD for 2000-01-01 used as QSAS default base in GUI feedback
#define MJD2000 51544

// MJD for 0000-01-01 (correctly Jan 01, BCE 1) 
// Julian proleptic calendar value. 
#define MJD_0000J -678943

// Note these routines do not use leap seconds.
// conversion to and from MJD will be correct except
// for the leap second itself which will alias to the first second of the following minute
// Intervals enclosing a leap second will appear to be 1 second short
// Conversions from and to TT and TAI are correct except for the leap second itself.

// NASA CDF TT2000 offset in seconds = 43200 - 32.184
// Terrestrial Time == TAI + 32.184
#define TT2000offset 43167.816

// Atomic clock (TAI) offset in seconds
#define TAI2000offset 43200

/* Gregorian proleptic calendar value.  (At MJD_0000J the Gregorian proleptic
   calendar reads two days behind the Julian proleptic calendar, i.e. - 2 days,
   see http://en.wikipedia.org/wiki/Proleptic_Gregorian_calendar,
   so MJD_0000G = MJD_0000J+2) */
#define MJD_0000G -678941

// MJD for Jan 01, 1970 00:00:00 used for seconds from Unix base epoch
#define MJD_1970 40587

// QSAS version number for using DvTime
#define QMJD_TIME_VERSION_START 20400

#define DEFAULT_DECIMAL_SEC 3
#define DEFAULT_GREGORIAN_START "1582-10-15"

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ class DvTime

// Date and time
// MJD measures from the start of 17 Nov 1858 */
	
//	 Default uses the Gregorian calendar after 4 Oct 1582 (Julian) i.e. from 15 Oct 1582 Gregorian
//	 Note unix libraries use Gregorian only from 14 Sept 1752 
//	 More detailed discussion can be found at http://aa.usno.navy.mil/data/docs/JulianDate.php
//	 These routines have been compared with the results of the US Naval Observatory online converter.
//	 Modified Julian Date (MJD) = Julian Date (JD) - 2400000.5
	 
//	 In all routines, specifying a day, hour, minute or second field greater than would be valid is
//	 handled with modulo arithmetic and safe. 
//	 Thus 2006-12-32 00:62:00.0 will safely, and correctly, be treated as 2007-01-01 01:02:00.0 

//-------------------------------------------------------------------------
class DvTime
{
  	
public:

// Date and time are hld in Modified Julian Date split into int a day (MJD day) and 
// seconds from start of that day. Seconds may exceed a day or be negative

	int base_day; /* integer part of MJD used as default */
	double time_sec; /* seconds from start of base_day */
    int system; /* DvUTC, DvTAI or DvTT */
	
// Constructors
  
	DvTime(double sec=0.0, int day=MJD2000); 	// default day is Jan 01, 2000 = MJD 51544
	DvTime( const char* iso_time);
	DvTime( const DvString& iso_time);
	DvTime( int year, int month=1, int day=1, int hour=0, int min=0, double sec=0.0);
	DvTime( const DvTime& ttime );
	
// Methods
	int setSecResolution(int resolution); // returns previous resolution
	int getSecResolution(); // returns previous resolution
	void setFromUT(int year, int month, int day, int hour, int min, double sec);
	void setFromDOY(int year, int doy, int hour, int min, double sec);
	void setFromBCE(int yearBCE, int month, int day, int hour, int min, double sec);
	void setFromMJD(double ModifiedJulianDate);
	void setFromCDFepoch(double cdfepoch);
	void setFromCDFepoch16(double *cdfepoch16);
	void setFromJD(double JulianDate);
    void setFromEpoch2000(double sec);
    void setFromTT2000(DvInt8 nsec);
    void setFromTAI2000(DvInt8 sec);
	void setFromEpoch1970(double epoch1970);
	bool setFromISOstring(const char* ISOstring);
	bool setFromISOstring(DvString ISOstring);
	int breakDownMJD(int *year, int *month, int *day, int *hour, int *min, double *sec) const;
	double getMJD();
	double getJD();
	double getDiffDays(const DvTime &MJD)const;
	double getDiffSecs(const DvTime &MJD)const;
	double getCDFepoch();
	void getCDFepoch16(double *epoch16);
	double getEpoch2000();
	double getEpoch1970();
    std::string getISOstring(int delim=0) const;
    std::string getTAIstring(int delim=0) const;
    std::string getTTstring(int delim=0) const;
	const char * getDayOfWeek();
	const char * getLongDayOfWeek();
	const char * getMonth(int m);
	const char * getLongMonth(int m);
	bool getYandD(int *year, int* day);
	size_t strfMJD(char * buf, size_t len, const char *format);
    DvTime normalized() const;
	void setGregorianStartMJD(double GregorianMJD);

	DvTime& operator=( const DvTime& ttime )
	{
		// assign from another Time
		base_day = ttime.base_day;
		time_sec = ttime.time_sec;

		return *this;
	}
   
	DvTime& operator=( const double dtime )
	{
		// assign from another Time
		this->setFromEpoch2000(dtime);

		return *this;
	}
   
	DvTime& operator=( const string& isotime )
	{
		// set from string in ISO time format
		if( !setFromISOstring(isotime.c_str()) ) 
		{
			base_day = MJD2000;
			time_sec = 0.0;
		}
		
		return *this;
	}
  
      
	// --- operators needed for intermediate steps (e.g. interpolation, average etc) 
	DvTime& operator*=(const double factor) 
	{
		// keep factor as close to unity as possible to retain accuracy
		double duration=this->getEpoch2000();
		duration *= factor;
		this->setFromEpoch2000(duration);

		return *this;
	} 

	DvTime& operator/=(const double factor) 
	{
		// keep factor as close to unity as possible to retain accuracy
       double duration=this->getEpoch2000();
       duration /= factor;

       this->setFromEpoch2000(duration);
	   return *this;
	} 
   
   
	DvTime& operator-=( DvTime t) 
	{
		// time object temporarily meaningless
		double duration = this->getEpoch2000() - t.getEpoch2000();
		this->setFromEpoch2000(duration);

		return *this;
	} 

	DvTime& operator+=( DvTime t) 
	{
		// time object temporarily meaningless
		double duration = this->getEpoch2000() + t.getEpoch2000();
		this->setFromEpoch2000(duration);
       
       return *this;
	} 
    // -----------------------------
   
	// --- add and subtract seconds 
	DvTime operator+(const double seconds) const
	{
       DvTime tres=*this;
       tres.time_sec += seconds;
       return tres;
	} 
   
	DvTime operator-(const double seconds) const
	{
       DvTime tres=*this;
       tres.time_sec -= seconds;
       return tres;
	}
       
	// --- add and subtract time and update
	DvTime& operator+=(const double seconds)
	{
       time_sec += seconds;
       return *this;
	}
   
	DvTime& operator-=(const double seconds)
	{
       time_sec -= seconds;
       return *this;
	}

	double operator-(const DvTime& ttime) const
	{
		int diffDays = base_day - ttime.base_day;
		double diffSecs = time_sec - ttime.time_sec;
		return diffSecs + diffDays * 86400.;
	}
   
	double operator=(const DvTime& ttime) const
	{
		return time_sec + base_day * 86400.;
	}
   
	// --- logical operators
	bool operator==(const DvTime& ttime) const
	{ 
		if( fabs(this->getDiffSecs(ttime)) < 1.e-10) return true;
        return false;
	} 
   
	bool operator!=(const DvTime& ttime) const
	{ 
   		if( fabs(this->getDiffSecs(ttime)) > 1.e-10) return true;
        return false;
	} 
   
	bool operator<(const DvTime& ttime) const
	{ 
   		if( this->getDiffSecs(ttime) < 0) return true;
        return false;
	}
   
	bool operator<=(const DvTime& ttime) const
	{ 
		double diff = this->getDiffSecs(ttime);
   		if(diff  < 0 || fabs(diff) < 1.e-10 ) return true;
        return false;
	}  
   
	bool operator>(const DvTime& ttime) const
	{ 
   		if( this->getDiffSecs(ttime) > 0) return true;
        return false;
	}
   
	bool operator>=(const DvTime& ttime) const
	{ 
		double diff = this->getDiffSecs(ttime);
   		if(diff  > 0 || fabs(diff) < 1.e-10 ) return true;
        return false;
	}  

   	string iso_srep() const
	{	
		return getISOstring(1);
	}
   
   DvTime& set_today()
   {
     time_t ltime;
     time(&ltime);
     tm *t=localtime(&ltime);
     setFromUT(t->tm_year+1900, t->tm_mon+1, t->tm_mday, 0, 0, 0.0);
     return *this;
   }

};

// ----------- End of class DvTime -------------------

    

// ++++++++++++++++++++++++++++++++++++++++++++++++++++ class DvEvent

class DvEvent 
{
  
public:
    typedef enum Sense_e {TIME_SENSE_NEGATIVE, NO_TIME_SENSE, TIME_SENSE_POSITIVE} Sense;

private:
    DvTime _t1;
    DvTime _t2;

public:

    // --- Constructors
    DvEvent():_t1(DvTime()),_t2(DvTime()+1.0){;}
    DvEvent( const DvEvent& tint): _t1(tint._t1),_t2(tint._t2){;}
    DvEvent( const DvTime& st, const DvTime& ed, Sense option=NO_TIME_SENSE );
    DvEvent( double st, double ed, Sense option=NO_TIME_SENSE ): _t1(DvTime(st)),_t2(DvTime(ed)){;}
    DvEvent( string &ISOstring );
	    
    // --- get start/end 
    
    DvTime start() const { return _t1;}
    DvTime end() const { return _t2;}
    DvTime start_abs() const { 
		if(_t1<_t2) return _t1;
		else return _t2;
	}
    DvTime end_abs() const { 
		if (_t1>_t2) return _t1;
		else return _t2;
	}  
	
	std::string getISOstring(int delim=0) const;

    // --- get duration 
    
    double duration() const { return _t2.getDiffSecs(_t1);}
    double duration_abs() const {return fabs(_t2.getDiffSecs(_t1));}
    DvTime centre() const {
		DvTime tc(_t1);
		tc += 0.5 * _t2.getDiffSecs(_t1);
		return tc;
	}
         
   
    // --- get sense of time direction:  
    // +1 if start<end, -1 if end<start, 0 if start==end
    Sense get_sense() const
    {
      if (_t1==_t2)
		return NO_TIME_SENSE; // does fuzzy test to 0.1 nanosec
      else if(_t1>_t2)
        return TIME_SENSE_NEGATIVE;
      else
        return TIME_SENSE_POSITIVE;
    }         

    // --- set by operator=
    DvEvent& operator=( const DvEvent& tint )
    {
        set(tint.start(),tint.end());
        return *this; 
    }

    // --- set start and end; default option is NO_TIME_SENSE i.e. store
    // start as _t1, end as _t2 regardless of whether _t1>_t2 or _t2<_t1 
    // TIME_SENSE_NEGATIVE: force _t1>=_t2
    // TIME_SENSE_POSITIVE: force _t1<=_t2
    void set( const DvTime& st, const DvTime& ed, Sense option=NO_TIME_SENSE)
    {
        switch(option)
        {
            case NO_TIME_SENSE:
                _t1=st;
                _t2=ed;
                break;
            case TIME_SENSE_NEGATIVE:
                 _t1 = (_t1>_t2) ? _t1:_t2;
                 _t2 = (_t1<_t2) ? _t1:_t2;
                break;
            case TIME_SENSE_POSITIVE:
                 _t1 = (_t1<_t2) ? _t1:_t2;
                 _t2 = (_t1>_t2) ? _t1:_t2;
         }
     }  
	
	bool setFromISOstring(const char* ISOstring)
	{
		string ISO_s1(ISOstring);
		if(ISO_s1.length() < 21) return false;
		
		string ISO_s2(ISOstring);
		if(ISO_s2.length() < 21) return false;
		
		size_t ptr = ISO_s1.find('/');
		if(ptr == string::npos) return false;
		
		ISO_s1 = ISO_s1.substr(0, ptr);
		ISO_s2.replace(0, ptr+1, "");
		
		if( !_t1.setFromISOstring(ISO_s1.c_str()) ) return false;
		if( !_t2.setFromISOstring(ISO_s2.c_str()) ) return false; 
		
		return true;
	}

	bool setFromISOstring(string ISOstring)
	{
		return setFromISOstring(ISOstring.c_str());
	}
	    
    // --- force ordering of start and end to be positive
    void force_sense( Sense option=TIME_SENSE_POSITIVE )
    {
         switch(option)
         {
             case TIME_SENSE_POSITIVE:
             if (_t1>_t2)
                 reverse();
             break;
             case TIME_SENSE_NEGATIVE:
             if (_t1 < _t2)
                 reverse();
             break;
             case NO_TIME_SENSE:
             break;
         }
    }
    
    // --- set start time (_t1) 
    void set_start( const DvTime& st)
    {     
        _t1=st;
    }     
    
     
    //  --- set end time (_t2)
	void set_end( const DvTime& ed)
    {
		_t2=ed;
    }     
    
    // --- reverse start and end
    void reverse()
    {
        DvTime ttemp=_t1;
        _t1=_t2;
        _t2=ttemp;
    }

    // shift interval by given amount of time
    DvEvent& operator+=(double offset)
    {
        _t1.time_sec += offset;
        _t2.time_sec += offset;
        return *this;
    }
        
    DvEvent& operator-=(double offset)
    {
        _t1.time_sec -= offset;
        _t2.time_sec -= offset;
        return *this;
    }

    // --- Set Operations
    DvEvent event_intersect(const DvEvent& ti) const;
    DvEvent event_union(const DvEvent& ti) const;

    bool contains( const DvTime& t ) const;
    bool intersects( const DvEvent& tivl ) const;
    bool empty() const { return _t1 == _t2;}

    double operator*(const double multiplier) const
    {
		double duration = duration_abs();
        return duration * multiplier;
    }
    
    // --- Logical Operations
    bool operator==(const DvEvent& ti) const
    {
        return (start()==ti.start() && end()==ti.end());
    }
    
    bool operator!=(const DvEvent& ti) const
    {
        return !(*this==ti);
    }
    
	bool operator>(const DvEvent& ti) const
    {
		if( start() > ti.start() ) return true;
	
		if( start() == ti.start() && end() > ti.end() ) return true;
		
		return false;

    }

    std::string iso_srep() const { return start().iso_srep() + " - " + end().iso_srep();}

};


// Utilities as DvTime is not just a double

inline DvTime max(DvTime &t1, DvTime &t2)
{
	if (t1 > t2) return t1;
	else return t2;
}
inline DvTime min(DvTime &t1, DvTime &t2)
{
	if (t1 < t2) return t1;
	else return t2;
}


#endif // #ifndef DV_TIME_H



