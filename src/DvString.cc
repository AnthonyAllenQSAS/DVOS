//
//  DvString.cc
//
//  DVOS - string class 
//

#include "DvObj.h"
#include "DvString.h"

  int DvString::dblPrecision=16; // static used for all strings when setting from double

  int DvString::getDblPrecision(){
	return dblPrecision;
  }
	
  void DvString::setDblPrecision(int prec){
	dblPrecision = prec;
  }
	
  void DvString::set(double value){
	std::ostringstream buffer;
  	buffer.precision(dblPrecision);
	buffer << value;
	
	*this = buffer.str();  	
  }

void DvString::set(double value, int dec){
    
    std::ostringstream buffer;
    if(value < 1.e-1){
        buffer.setf (std::ios::scientific, std::ios::floatfield);
        buffer.precision(dec);
    }
    else if(value < 10000){
        buffer.setf (std::ios::fixed, std::ios::floatfield);
        buffer.precision(dec);
    }
    else {
        buffer.setf (std::ios::scientific, std::ios::floatfield);
        buffer.precision(dec);
    }
    buffer << value;
    
    *this = buffer.str();
}

  void DvString::append(double d, const char * spacer){
    *this += spacer;
	ostringstream tmpStr(ostringstream::ate);
	tmpStr.precision(dblPrecision);
	
	tmpStr << d;

  	*this += tmpStr.str();
  }

void DvString::append(double value, int dec, const char * spacer){

    *this += spacer;
    
    std::ostringstream buffer;
    if(value < 1.e-1){
        buffer.setf (std::ios::scientific, std::ios::floatfield);
        buffer.precision(dec);
    }
    else if(value < 10000){
        buffer.setf (std::ios::fixed, std::ios::floatfield);
        buffer.precision(dec);
    }
    else {
        buffer.setf (std::ios::scientific, std::ios::floatfield);
        buffer.precision(dec);
    }
    buffer << value;
    
    *this += buffer.str();
}
