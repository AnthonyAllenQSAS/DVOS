//
//  DvString.h
//
//  DVOS - string class 
//

#ifndef DV_STRING_H
#define DV_STRING_H

#include <iostream>
#include <sstream>
#include <string.h>

using namespace std;

class DvTime;
class DvEvent;

class DvString : public string
{
public:
  DvString() : string() {;}
  DvString( const string& s ) : string(s) {;}
  DvString( const DvString& ds ) : string(ds.c_str()) {;}
  DvString( const char c ) : string(" ") {(*this)[0]=c;}
  DvString( const char* cs ) : string(cs) {;}
  DvString( const int i ) : string() { this->append(i);}
  DvString( const double d ) : string() { this->append(d);}
  DvString( DvTime t ) : string() { this->append(t);}
  DvString( DvEvent e ) : string() { this->append(e);}
  DvString( const char* cs, int len ) : string(cs, (size_t) len) {;}

  ~DvString(){;}

  DvString& operator=( string& s )
    { (string&)(*this) = s; return *this; }
  DvString& operator=( const char* cs )
    { (string&)(*this) = cs; return *this; }
  DvString& operator=( DvString ds )
    { (string&)(*this) = ds.c_str(); return *this; }
    DvString& operator=( const int i )
    { (string&)(*this) = ""; this->append(i); return *this; }
    DvString& operator=( const size_t i )
    { (string&)(*this) = ""; this->append((int)i); return *this; }
  DvString& operator=( const double d )
    { (string&)(*this) = ""; this->append(d); return *this; }
  DvString& operator=( DvTime &t )
    { (string&)(*this) = ""; this->append(t); return *this; }
  DvString& operator=( DvEvent &e )
    { (string&)(*this) = ""; this->append(e); return *this; }

 static int dblPrecision;

  int getDblPrecision();	
  void setDblPrecision(int prec);
  	
    
  DvString substr(size_t i, size_t len){
    return DvString( ((string)*this).substr(i, len));
  }
    
  DvString substr(size_t i){
    return DvString( ((string)*this).substr(i));
  }
    
  void setQuoted(string& s){
	if(s[0]=='\"') (string&)(*this) = s.substr(1, s.size()-2);
	else (string&)(*this) = s;
  }
  void setQuoted(const char * cs){
	if(cs[0]=='\"') (string&)(*this) = string( &(cs[1]), strlen(cs)-2);
	else (string&)(*this) = cs;
  }
  
  bool blank(){
	size_t pos = 0;
  	while ( pos < this->size()){
		if(this->c_str()[pos] == '\"' || this->c_str()[pos] == ' ' || this->c_str()[pos] == '\t' ) pos++;
		else return false;
	}
	return true;
  }
  
  char last(){
  	return this->c_str()[this->length()-1];
  }
  
  DvString before(const char delim){
  	size_t ptr = this->find(delim);
	if (ptr == string::npos) return DvString();
	return this->substr(0, ptr);
  }
  DvString beforeR(const char delim){
  	size_t ptr = this->rfind(delim);
	if (ptr == string::npos) return DvString();
	return this->substr(0, ptr);
  }
  DvString before(const char *delim){
  	size_t ptr = this->find(delim);
	if (ptr == string::npos) return DvString();
	return this->substr(0, ptr);
  }
  DvString before(size_t posn){
	if (posn >= this->length()) return *this;
	return this->substr(0, posn);
  }
  
  DvString afterR(const char delim){
  	size_t ptr = this->rfind(delim);
	if (ptr == string::npos) return DvString();
	return this->substr(ptr+1);
  }
  DvString after(const char delim){
  	size_t ptr = this->find(delim);
	if (ptr == string::npos) return DvString();
	return this->substr(ptr+1);
  }
  DvString after(const char *delim){
  	size_t ptr = this->find(delim);
	if (ptr == string::npos) return DvString();
	ptr += strlen(delim);
	return this->substr(ptr);
  }
  DvString after(size_t posn){
	if (posn >= this->length()) return DvString();
	return this->substr(posn+1);
  }
  DvString between(const char *delim1, const char *delim2){
	return this->after(delim1).before(delim2);
  }
  DvString between(const char delim1, const char delim2){
	return this->after(delim1).before(delim2);
  }
  DvString betweenIncl(const char *delim1, const char *delim2){
	if (this->empty() ) return DvString();

  	DvString ret = this->after(delim1);
	ret = ret.before(delim2);
	if (ret.empty() ) return ret;
	
	DvString txt(delim1);
	txt += ret;
	txt += delim2;
	return txt;
  }
  DvString betweenIncl(const char delim1, const char delim2){
  
  	DvString ret = this->after(delim1).before(delim2);
	if (ret.empty() ) return ret;
	
	ret += delim2;
	return DvString(delim1) + ret;
  }

  DvString front(const char *delim){
  	size_t ptr = this->find(delim);
	if (ptr == string::npos) return *this;

	ptr += strlen(delim)-1;
	return this->substr(0, ptr);
  }
  DvString back(const char *delim){
  	size_t ptr = this->find(delim);
	if (ptr == string::npos) return *this;

	ptr += strlen(delim)-1;
	return this->substr(ptr+1);
  }
  
  DvString common(DvString &str){

  	DvString comStr;
	size_t posn = 0;
	size_t posn1 = 0;
	while(posn < this->size() && posn1 < str.size()){
		if(this->at(posn) == str.at(posn1)) {
			comStr += this->at(posn);
			posn++;
			posn1++;
		}
		else{
		
			while(posn < this->size() && this->at(posn) != '(' && this->at(posn) != '_') 	posn++;
			while(posn1 < str.size()  && posn < this->size()  && str.at(posn1) != this->at(posn)) posn1++;
		}
	}

	return comStr;
  }
  
  DvString replace(const char * str1, const char * str2){
  	size_t posn = this->find(str1);
	DvString ret(*this);
	if (posn != string::npos) ret = ret.string::replace(posn, strlen(str1), str2);
	return ret;
  }
  
  DvString replaceAll(const char *str1, const char * str2){
  	size_t posn = this->find(str1);
	DvString ret(*this);
	while (posn != string::npos){
		ret = ret.string::replace(posn, strlen(str1), str2);
		posn = ret.find(str1);
	}
	return ret;
  }
  
  DvString replaceTokens(const char * str, const char delim='?'){
  	DvString processStr(str);
  	DvString token = processStr.between(delim, delim);
	DvString ret(*this);
	while( !token.empty()){
		DvString newStr = token.after('=');
		DvString oldStr = token.before('=');
		// modify output
		if(!newStr.empty() && !oldStr.empty()) ret = ret.replaceAll(oldStr.c_str(), newStr.c_str());
		// next token (delim occurs twice)
		processStr = processStr.after(delim);
		processStr = processStr.after(delim);
		token = processStr.between(delim, delim);
	}
	
	return ret;
  }

  double toDouble(){
	
	// also handle Latex syntax numbers in some metadata
	
	double value=0.0;
	size_t posn = this->find("10^");
	if (posn != string::npos){
	
		DvString first = this->before(posn);
    	value = first.toDouble();
    	if(value < 1.e-13) value = 1.0; // if vtxt of form "10^6" value is zero since there is no mantissa
		
		DvString second = this->after(posn+2);
		second = second.before('>'); // remove unwanted text
    	double power = strtod( second.c_str(), 0);
		 
    	value *= pow(10., power);
		return value;
  	}

	return strtod( this->c_str(), 0); // normal C formatted numbers
  	
  }
  
  int toInt(){
	return strtol( this->c_str(), 0, 10);
  }  
  
  bool isNumber(){
	size_t i=0;
	while(i < this->size()){
		char c = this->at(i);
		if( !isdigit(c) && c != 'e' && c != 'E' && c != '+' && c != '-' && c != ' ' && c != '.' ) return false;
		i++;
	}	
	return true;
  }
  
  bool isInt(){
	size_t i=0;
	while(i < this->size()){
		char c = this->at(i);
		if( !isdigit(c) && c != '+' && c != '-' && c != ' ' ) return false;
		i++;
	}	
	return true;
  }
  
  int toBool(){
  	if( *this == "TRUE") return true;
  	else if( *this == "FALSE") return false;
	else return (bool) strtol( this->c_str(), 0, 10);
  }
  
  DvTime toTime(){
	return DvTime(*this);
  }
  
  DvEvent toEvent(){
	return DvEvent(*this);
  }
  
    void set(double value);

    void set(double value, int prec);
  
  void set(int value){
	std::ostringstream buffer;
	buffer << value;

	*this = buffer.str();  	
  }
  
    void append(double d, const char * spacer="");

    void append(double d, int prec, const char * spacer="");
  
  void append(int i, const char * spacer=""){
    *this += spacer;
	ostringstream tmpStr(ostringstream::ate);
	tmpStr << i;

  	*this += tmpStr.str();
  }
  void append(size_t posn, const char * spacer=""){
     *this += spacer;
	 unsigned int i = posn;
	ostringstream tmpStr(ostringstream::ate);
	tmpStr << i;

  	*this += tmpStr.str();
  }
  void append(DvTime &t, const char * spacer=""){
    *this += spacer;
  	(string &)(*this) += t.getISOstring();
  }
  void append(DvEvent &e, const char * spacer=""){
    *this += spacer;
  	(string &)(*this) += e.getISOstring();
  }
  
  DvString operator+( double d )
    {  DvString newStr(*this); newStr.append(d); return newStr; }
  
  DvString operator+( int i )
    {  DvString newStr(*this); newStr.append(i); return newStr; }
  
  DvString operator+( DvTime &t )
    {  DvString newStr(*this); newStr.append(t); return newStr; }
  
  DvString operator+( DvEvent &e )
    {  DvString newStr(*this); newStr.append(e); return newStr; }
  
  DvString operator+( DvString &vs )
    {  DvString newStr(*this); newStr += vs; return newStr; }
  
  DvString operator+( string s )
    {  DvString newStr(*this); newStr += s; return newStr; }
  
  DvString operator+( const char *cs )
    {  DvString newStr(*this); newStr += cs; return newStr; }
  
  DvString operator+( const char c )
    {  DvString newStr(*this); newStr += c; return newStr; }
  
  DvString& operator+=( double d )
    {  this->append(d); return *this; }
  
    DvString& operator+=( int i )
    {  this->append(i); return *this; }
    
    DvString& operator+=( size_t i )
    {  this->append((int)i); return *this; }
    
  DvString& operator+=(DvTime &t )
    {  this->append(t); return *this; }
  
  DvString& operator+=( DvEvent &e )
    {  this->append(e); return *this; }
  
  DvString& operator+=( DvString &vs )
    {  (string&)(*this) += (string&)vs; return *this; }
  
  DvString& operator+=( string s )
    {  (string&)(*this) += (string&)s; return *this; }
  
  DvString& operator+=( const char *cs )
    {  (string&)(*this) += string(cs); return *this; }

  DvString& operator+=( const char c )
    {  (string&)(*this) += c; return *this; }

  // get single element as a string
  DvString operator()( int i )
    {  return this->substr(i, 1); }
	  
  bool operator==(  DvString vs )
    {  // case insensitive
		if(this->toUpper().compare(vs.toUpper()) == 0) return true;
		return false; 
	}
  
  bool operator==( const char *cs )
    {  // case insensitive
		return *this == DvString(cs); 
	}

  bool operator==( string ss )
    {  // case insensitive
		return *this == DvString(ss.c_str()); 
	}
  
  bool operator!=(  DvString vs )
    {  // case insensitive
		if(this->toUpper().compare(vs.toUpper()) == 0) return false;
		return true; 
	}
  
  bool operator!=( const char *cs )
    {  // case insensitive
		return *this != DvString(cs); 
	}
  
  bool operator!=( string ss )
    {  // case insensitive
		return *this != DvString(ss.c_str()); 
	}
  
  
  string &str(){return (string &) *this;}
  
  const char * cstr(){
  	// returns c_str() of quoted text without quotes and space padding outside of quotes
  	return this->clean().c_str();
  }
  
  DvString &clean(){
  	// modifies quoted text removing quotes and space padding outside of quotes
	while((*this)[0] == ' ') this->erase(0,1); 
	if( (*this)[0] == '"') this->erase(0,1); 

	while((*this)[this->length()-1] == ' ') this->erase(this->length()-1,1); 
	if( (*this)[this->length()-1] == '"') this->erase(this->length()-1,1); 
	
	return *this;
  }
  
  
  size_t ptr(DvString &vs){
    // case insensitive
	return this->toUpper().find(vs.toUpper());
  }

  bool contains(DvString &vs){
	if(this->toUpper().find(vs.toUpper()) == string::npos) return false;
	return true;
  }

  bool containsCS(DvString &vs){
	if(this->find(vs) == string::npos) return false;
	return true;
  }
  
  bool contains(const char *cs){
  	DvString vs(cs);
  	bool found = contains(vs);
	return found;
  }

  bool containsCS(const char *cs){
  	DvString vs(cs);
  	bool found = containsCS(vs);
	return found;
  }
  
  bool contains(const char c){
	if(this->find(c) == string::npos) return false;
	return true;
  }
  
  bool endsIn(DvString &vs){
  	int vsLen = vs.length();
  	int thisPosn = this->length() - vsLen;
  	if(thisPosn < 0) return false;
	
	for(size_t i=0; i<vs.length(); i++){
		if((*this)[thisPosn+i] != vs[i]) return false;
	}
	
	return true;
  }
  
  bool endsIn(const char *cs){
	DvString vs(cs);
  	return endsIn(vs);
  }

  void strip(const char c){
    size_t pos = this->find(c);
  	while ( pos != string::npos){
		this->erase(pos, 1);
		pos = this->find(c);
	}
  }
  
  DvString & trimmed(){
	
  	// remove leading and trailing space
	while (this->length() != 0 && this->at(0) == ' ' )	{
		this->erase(0, 1);
		
	}  	
	while(length() != 0 && this->at(length()-1) == ' ')	{
			this->erase(length()-1, 1);
	}
	return *this;
  }
  
  DvString toUpper(){
	DvString upper=*this;
	for (size_t i=0; i < upper.length(); i++)  upper[i] = toupper(upper[i]);
	return upper;
  }


};


#endif //  DV_STRING_H

