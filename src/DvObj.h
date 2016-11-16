//
//  DvObj.h
//
//  Dvos -  Data variable object classes
//
//  Tony Allen
//  Nov 2012 
//  Mar 2016
//  July 2016

#ifndef DVOBJ_H
#define DVOBJ_H

#define DVOS_VERSION_STR "1.1.0"
#define DVOS_VERSION_DATE "01 Nov 2016"

#define CACHE_XREF_CHAR '_'

#define DV_JOIN_METHOD "join method"
#define DV_GAP_WIDTH "join gap width"
#define DV_BOX_WIDTH "join boxcar width"
#define DV_GAP_OPTION "join gap option"
#define DV_END_OPTION "join end option"
#define DV_MIN_BOXCAR "minimum boxcar number"
#define DV_LINEAR "Linear"
#define DV_BOXCAR "Boxcar"
#define DV_NN "Nearest Neighbour"
#define DV_REMOVE "Remove"
#define DV_FILL "Fill"
#define DV_ZERO_FILL "Zero Fill"
#define DV_DEF_INTEGRAL "Definite Integral"
#define DV_INDEF_INTEGRAL "Indefinite Integral"
#define DV_DERIV_3PT "3-point estimate"
#define DV_DERIV_5PT "5-point estimate"

// component selection constraints
#define DV_NONE "None"
#define DV_SCALAR "Scalar"
#define DV_TIME "Time"
#define DV_EVENT "Event"
#define DV_ARRAY "Array"
#define DV_INT "Integer"
#define DV_DOUBLE "Double"

// component possibilities
#define DV_METADATA "Metadata..."
#define DV_SEPARATOR "-------"
#define DV_RMAG "R (mag)"
#define DV_X "X"
#define DV_Y "Y"
#define DV_Z "Z"
#define DV_LAT_DEG "Lat (deg)"
#define DV_THETA_DEG "theta (deg)"
#define DV_PHI_DEG  "phi (deg)"
#define DV_LAT_RAD "Lat (rad)"
#define DV_THETA_RAD "theta (rad)"
#define DV_PHI_RAD "phi (rad)"
#define DV_LOCAL_T "Local time"
#define DV_T_START "Start"
#define DV_T_CENTRE "Centre"
#define DV_T_END "End"
#define DV_ARRAY_REDUCE "Reduce Array Dimensions..."
#define DV_SUBSEQUENCE "Subset Sequence..."

#include <iostream>
#include <string.h>
#include <valarray>
#include <vector>
#include <cmath>
#include <complex>
#include <limits>
#include "DvUnit.h"
#include "Xrefs.h"
#include "DvTime.h"
#include "DvString.h"
#include "Dvar.h"

using namespace std;

extern const double DV_CLOSE_ENOUGH;
extern const double DvNaN;
extern const vector <size_t> emptyDim;
extern const vector <size_t> DIM3;
extern const vector <size_t> DIM33;
extern const vector <size_t> DIM4;
extern const vector <size_t> DIM44;


class DvObject;
class DvRecord;

// XML saveset utility calls

char * readXML(DvString fullName);
char* getNextTag(char *xml, size_t &start, DvString &tag);

// var object version

typedef Dvar<DvObject> DvObject_var;

//-----------------------------------------------------------------------------
// class DvMask
//-----------------------------------------------------------------------------

class DvMask  {
	valarray <bool> msk;
	// true means keep record
  public:	
	DvMask(DvMask &mask){ msk.resize(mask.size(), true); msk = mask.msk;}
	DvMask(int i, bool ok=true){ msk.resize(i, ok); }
    DvMask(){;} // empty mask
	
	~DvMask(){;}
	size_t size(){return msk.size();}
	size_t resultSize(){
		size_t res=0;
		for(size_t i=0; i<msk.size(); i++) if(msk[i]) res++;
		return res;
	}
	void resize(size_t n, bool ok=true){msk.resize(n, ok);}
	bool &operator[](size_t i){return msk[i];}
};



//-----------------------------------------------------------------------------
// class DvNode     
//-----------------------------------------------------------------------------

class DvNode {
	DvObject_var _obj;
	DvString _name;
	
	
  public:	
	DvNode(const char *name, DvObject_var &obj);
	DvNode(string name, DvObject_var &obj);
	DvNode(DvString name, DvObject_var &obj);

	DvNode *next;
		
	~DvNode(){if(next) delete next;}
	
//	DvNode operator=(DvNode xref_obj);

	DvString name(){return _name;}
	DvObject_var &obj(){return _obj;}
	void setObj(DvObject_var obj){_obj=obj.ptr();}
	void setName(DvString name){_name=name;}
	
	void toXML(std::ofstream &ofp);
	void fromXML(char *xml);
};

//-----------------------------------------------------------------------------
// class DvList
//-----------------------------------------------------------------------------

class DvList {
	DvNode *nodes; // uses DvNode linked list
		
	public:
		DvList(){nodes = NULL;}

		DvNode *first(){return nodes;}
		
		~DvList(){delete nodes;} // deletes linked list recursively
		
		DvNode *makeNode();
		
		void remove(DvString name){  	
  			DvNode *node = nodes;
			DvNode *prev = NULL;
			while(node){
				if(node->name() == name) {
					if(prev) prev->next = node->next;
					else nodes = node->next;
					node->next = NULL; // stop delete propagating
					delete node;
					return;
				}
				prev = node;
				node = node->next;
			}
		}
		
		DvString append(DvString name, DvObject_var &obj){
			// actually it prepends
			if( find(name) != NULL ){
				name += "A";
				append(name, obj);
			}
			else{
				DvNode *newNode = new DvNode(name, obj);
                DvNode *lastNode = last();
				if(lastNode) lastNode->next = newNode;
                else nodes = newNode;
			}
			return name;
		}
		
		DvNode *find( DvString& name ){
			DvNode *node = nodes;
			while(node){
				if(node->name() == name) return node;
				node = node->next;
			}
			return NULL;
  		}

		DvNode *find( const char * name ){
			DvNode *node = nodes;
			while(node){
				if(node->name() == name) return node;
				node = node->next;
			}
			return NULL;
  		}
		
        DvNode *last(){
            DvNode *node = nodes;
            while(node){
                if(node->next == NULL) return node;
                node=node->next;
            }
            return nodes;
        }
    
		int size(){
			DvNode *node = nodes;
			int count = 0;
			while(node){
				count++;
				node = node->next;
			}
			return count;
		}
		
		void clear(){
			DvNode *node = nodes;
			while(node){ 
				DvNode *child = node->next;
				delete node;
				node = child;
			}
			nodes = NULL;
		}
	
		void toXML(std::ofstream &ofp, DvString ListName);
		void fromXML(char *xml);							
};



//-----------------------------------------------------------------------------
// class DvObject
//-----------------------------------------------------------------------------
//
// Not template class to avoid trouble with creating var pointer of template


class DvObject  {

  static int ID_counter;
  unsigned int ref_count; // used by Dvar to know when safe to delete
  
  
  friend class DvRecord;

  protected:

  	// sequence of arrays unpacked
  	valarray <double> Ddata;
  	valarray <int> Idata;
  	valarray <DvString> Sdata;
  	valarray <DvTime> Tdata;
  	valarray <DvEvent> Edata;	
	
	DvList xrefs; // xrefs attached to this object (singly linked list)
	
	
  public:


  int object_id; // unique Data Variable object identifier
  vector <size_t> dims;  // array dimensions
  size_t seqLen;
  vector <DvString> errList;
  int secRes; // number of valid digits after decimal point for time

  // var counters for use by Dvar
  unsigned int get_ref_count() const { return ref_count; }
  unsigned int inc_ref_count()       { return ++ref_count; }
  unsigned int dec_ref_count()       { return --ref_count; }

  // enquiries
  int get_id() const { return object_id; }

  bool is_ok(){ return (object_id != -1); }
  bool is_nil(){ return (object_id == -1); }
  bool is_dbl(){ return (Ddata.size() > 0);}
  bool is_int(){ return (Idata.size() > 0);}
  bool is_str(){ return (Sdata.size() > 0);}
  bool is_time(){ return (Tdata.size() > 0);}
  bool is_event(){ return (Edata.size() > 0);}
  bool not_dbl(){ return !is_dbl();}
  bool not_int(){ return !is_int();}
  bool not_str(){ return !is_str();}
  bool not_time(){ return !is_time();}
  bool not_event(){ return !is_event();}
  
  size_t seqSize() const { return seqLen; } // length of sequence for data sequence
  vector <size_t> & Dims() {return dims;} // dimensions of array in sequence
  size_t nDims(){return dims.size();} // not normally needed as dims itself is needed too
  
  size_t arraySize(){ // length of array in sequence
	size_t len=1;
  	for(size_t i=0; i<dims.size(); i++) len *= dims[i];
  	return len;
  }
  size_t totalSize(){
  	if(Ddata.size() > 0) return Ddata.size(); 
	else if(Idata.size() > 0) return Idata.size();
	else if(Sdata.size() > 0) return Sdata.size();
	else if(Tdata.size() > 0) return Tdata.size();
	else if(Edata.size() > 0) return Edata.size();
	else return 0;
  }
  
  size_t maxStrLen(){
  	size_t len = 1;
  	if(this->is_str()){
		for(size_t i=0; i<Sdata.size(); i++) if(Sdata[i].length() > len) len = Sdata[i].length();		
	}
	return len;
  }
  

  // empty constructors
  DvObject();

  // create sequence nRecs of value
  // Defaults are single value constructors
  DvObject(double value, size_t nRecs=1);
  DvObject(int value, size_t nRecs=1);
  DvObject(DvString value, size_t nRecs=1);
  DvObject(const char * value, size_t nRecs=1);
  DvObject(DvTime value, size_t nRecs=1);
  DvObject(DvEvent value, size_t nRecs=1);

  // create sequence nRecs of values in valarray 
  DvObject(valarray <double> &from);
  DvObject(valarray <int> &from);
  DvObject(valarray <DvString> &from);
  DvObject(valarray <DvTime> &from);
  DvObject(valarray <DvEvent> &from);

  // create object of dimensions dims and length nRecs filled with value 
  DvObject(double value, const vector <size_t> &dims, size_t nRecs=1);
  DvObject(int value, const vector <size_t> &dims, size_t nRecs=1);
  DvObject(const char *value, const vector <size_t> &dims, size_t nRecs=1);
  DvObject(DvString value, const vector <size_t> &dims, size_t nRecs=1);
  DvObject(DvTime value, const vector <size_t> &dims, size_t nRecs=1);
  DvObject(DvEvent value, const vector <size_t> &dims, size_t nRecs=1);

  // copy constructors
  DvObject(DvObject_var &inObj, bool withXrefs=true);
  DvObject(DvObject &inObj, bool withXrefs=true);
  DvObject(DvObject *inObj, bool withXrefs=true);
  DvObject(DvObject_var &inObj, size_t nRecs, bool withXrefs=true);
  DvObject(DvObject *inObj, size_t nRecs, bool withXrefs=true);
  DvObject(valarray <double> &from, vector <size_t> &dims, size_t nRecs=1);
  DvObject(valarray <int> &from, vector <size_t> &dims, size_t nRecs=1);
  DvObject(valarray <DvString> &from, vector <size_t> &dims, size_t nRecs=1);
  DvObject(valarray <DvTime> &from, vector <size_t> &dims, size_t nRecs=1);
  DvObject(valarray <DvEvent> &from, vector <size_t> &dims, size_t nRecs=1);
  
  // assign regular values in steps of h (for Time or double)
    void assign_regular(double h);
    void assign_regular(int h);

  // copy single record to new object (with correct xrefs)
  DvObject_var subSequence(size_t nRec);

  // copy record range to new object
  DvObject_var subSequence(size_t fromRec, size_t toRec);
  
  // copy data within event range to new object
  DvObject_var subSequence(DvEvent &range);
  
  // copy data for slice to new object
  DvObject_var subSequence(slice &sl, size_t len=0);
  
  // strip out records from data and xrefs
  void apply_mask(DvMask &msk);

  
  // Locate bounding records for interval and time
  bool getEventRecBound(int &recL, int &recU, DvEvent &range);
  bool getRecBounds(int &recL, int &recU, double start, double end);
  void getBounds(int &recL, int &recU, const DvTime &t);   // recursive algorithm, so gets both records 
  void getBounds(int &recL, int &recU, double d);   // recursive algorithm, so gets both records 

  // create object of same type, different dimensions (fill with default element)
  DvObject *create(const vector <size_t> &dims=emptyDim, size_t nRecs=1);
  DvObject *create(size_t nRecs=1){ return this->create(this->dims, nRecs);} // scalar seq
  DvObject_var resize(size_t nRecs);
    
  // replace data valarray
  void replaceData(valarray <int> &data, size_t len=0){
		if(len == 0) len = data.size();
  		Idata.resize(len, 0);
		Idata = data[std::slice(0, len, 1)];
		if(is_dbl()) Ddata.resize(0, 0.0);
		else if(is_str()) Sdata.resize(0, "");
		else if(is_time()) Tdata.resize(0, DvTime());
		else if(is_event()) Edata.resize(0, DvEvent());
  }
  void replaceData(valarray <double> &data, size_t len=0){
		if(len == 0) len = data.size();
  		Ddata.resize(len, 0.0);
		Ddata = data[std::slice(0, len, 1)];
		if(is_int()) Idata.resize(0, 0);
		else if(is_str()) Sdata.resize(0, "");
		else if(is_time()) Tdata.resize(0, DvTime());
		else if(is_event()) Edata.resize(0, DvEvent());
  }
  void replaceData(valarray <DvString> &data, size_t len=0){
		if(len == 0) len = data.size();
  		Sdata.resize(len, DvString(""));
		Sdata = data[std::slice(0, len, 1)];
		if(is_int()) Idata.resize(0, 0);
		else if(is_dbl()) Ddata.resize(0, 0.0);
		else if(is_time()) Tdata.resize(0, DvTime());
		else if(is_event()) Edata.resize(0, DvEvent());
  }
  void replaceData(valarray <DvTime> &data, size_t len=0){
		if(len == 0) len = data.size();
  		Tdata.resize(len, DvTime());
		Tdata = data[std::slice(0, len, 1)];
		if(is_int()) Idata.resize(0, 0);
		else if(is_dbl()) Ddata.resize(0, 0.0);
		else if(is_str()) Sdata.resize(0, "");
		else if(is_event()) Edata.resize(0, DvEvent());
  }
  void replaceData(valarray <DvEvent> &data, size_t len=0){
		if(len == 0) len = data.size();
  		Edata.resize(len, DvEvent());
		Edata = data[std::slice(0, len, 1)];
		if(is_int()) Idata.resize(0, 0);
		else if(is_dbl()) Ddata.resize(0, 0.0);
		else if(is_str()) Sdata.resize(0, "");
		else if(is_time()) Tdata.resize(0, DvTime());
  }
  
  ~DvObject(){
  }

  
  int getSecResolution(){
	int len = 0;
	DvString value;
	
	if(this->is_time()){
		// not perfect since seconds is a double
		value = Tdata[0].time_sec;
	}
	else if(this->is_event()){
		// not perfect since seconds is a double
		value = Edata[0].start().time_sec;
	}
	size_t posn = value.rfind('.');
	if(posn == string::npos) return 0;
	posn++;
	while( isdigit(value[posn+len]) ) len++;
	return len;
  }


  // Fast access to individual elements (allows modification of value)
  
  valarray<double>& d_valarray(){return Ddata;}
  
  double & dbl(size_t nRec, size_t i, size_t j, size_t k, size_t m);
  double & dbl(size_t nRec, size_t i, size_t j, size_t k);
  double & dbl(size_t nRec, size_t i, size_t j);
  double & dbl(size_t nRec, size_t i);
  double & dbl(size_t posn=0); // element at this location in valarray, single scalar can use dbl()
  double & dbl(size_t nRec, vector <size_t> &index); // element at record and index posn
  
  int & itg(size_t nRec, size_t i, size_t j, size_t k, size_t m);
  int & itg(size_t nRec, size_t i, size_t j, size_t k);
  int & itg(size_t nRec, size_t i, size_t j);
  int & itg(size_t nRec, size_t i);
  int & itg(size_t posn=0); // element at this location in valarray, single scalar can use itg()
  int & itg(size_t nRec, vector <size_t> &index); // element at record and index posn

  DvString & str(size_t nRec, size_t i, size_t j, size_t k, size_t m);
  DvString & str(size_t nRec, size_t i, size_t j, size_t k);
  DvString & str(size_t nRec, size_t i, size_t j);
  DvString & str(size_t nRec, size_t i);
  DvString & str(size_t posn=0); // element at this location in valarray, single string can use str()
  DvString & str(size_t nRec, vector <size_t> &index); // element at record and index posn

  DvTime & time(size_t nRec, size_t i, size_t j, size_t k, size_t m);
  DvTime & time(size_t nRec, size_t i, size_t j, size_t k);
  DvTime & time(size_t nRec, size_t i, size_t j);
  DvTime & time(size_t nRec, size_t i);
  DvTime & time(size_t posn=0); // element at this location in valarray, single time can use time()
  DvTime & time(size_t nRec, vector <size_t> &index); // element at record and index posn

  DvEvent & event(size_t nRec, size_t i, size_t j, size_t k, size_t m);
  DvEvent & event(size_t nRec, size_t i, size_t j, size_t k);
  DvEvent & event(size_t nRec, size_t i, size_t j);
  DvEvent & event(size_t nRec, size_t i);
  DvEvent & event(size_t posn=0); // element at this location in valarray, single event can use event()
  DvEvent & event(size_t nRec, vector <size_t> &index); // element at record and index posn


  // access record (used with DvRecord operators)
  DvRecord operator[](size_t nRec);
  
  // Copy back double array at record
  valarray <double> record(size_t nRec);
  
    // SAFE versions of fast access to individual elements (allows modification of value)
  double & dblS(size_t nRec, size_t i, size_t j, size_t k, size_t m, double safe=0.);
  double & dblS(size_t nRec, size_t i, size_t j, size_t k, double safe=0.);
  double & dblS(size_t nRec, size_t i, size_t j, double safe=0.);
  double & dblS(size_t nRec, size_t i, double safe=0.);
  double & dblS(size_t posn=0, double safe=0.); // element at this location in valarray, single scalar can use dbl()

  int & itgS(size_t nRec, size_t i, size_t j, size_t k, size_t m, int safe); // cannot set default, overloading gets confused
  int & itgS(size_t nRec, size_t i, size_t j, size_t k, int safe);
  int & itgS(size_t nRec, size_t i, size_t j, int safe);
  int & itgS(size_t nRec, size_t i, int safe);
  int & itgS(size_t posn=0, int safe=0); // element at this location in valarray, single scalar can use itg()

  DvString & strS(size_t nRec, size_t i, size_t j, size_t k, size_t m, DvString safe=DvString());
  DvString & strS(size_t nRec, size_t i, size_t j, size_t k, DvString safe=DvString());
  DvString & strS(size_t nRec, size_t i, size_t j, DvString safe=DvString());
  DvString & strS(size_t nRec, size_t i, DvString safe=DvString());
  DvString & strS(size_t posn=0, DvString safe=DvString()); // element at this location in valarray, single string can use str()

  DvTime & timeS(size_t nRec, size_t i, size_t j, size_t k, size_t m, DvTime safe=DvTime());
  DvTime & timeS(size_t nRec, size_t i, size_t j, size_t k, DvTime safe=DvTime());
  DvTime & timeS(size_t nRec, size_t i, size_t j, DvTime safe=DvTime());
  DvTime & timeS(size_t nRec, size_t i, DvTime safe=DvTime());
  DvTime & timeS(size_t posn=0, DvTime safe=DvTime()); // element at this location in valarray, single time can use time()

  DvEvent & eventS(size_t nRec, size_t i, size_t j, size_t k, size_t m, DvEvent safe=DvEvent());
  DvEvent & eventS(size_t nRec, size_t i, size_t j, size_t k, DvEvent safe=DvEvent());
  DvEvent & eventS(size_t nRec, size_t i, size_t j, DvEvent safe=DvEvent());
  DvEvent & eventS(size_t nRec, size_t i, DvEvent safe=DvEvent());
  DvEvent & eventS(size_t posn=0, DvEvent safe=DvEvent()); // element at this location in valarray, single event can use event()

  // operate on records for any type, dims and type MUST match (fast, no checking done)
  void setRecord(size_t nRec, size_t oRec, DvObject* obj); // set rec nRec of this to record oRec of obj

  // conversions
  valarray <double> toDouble(double safeDefault=0.0);
  valarray <int> toInt(int safeDefault=0);
  valarray <DvString> toStr(DvString safeDefault=DvString(""));
  valarray <DvTime> toTime(DvTime safeDefault=DvTime());
  
  
  // single element conversions (first element of first record, convenience call)
  double asDouble(size_t posn=0, double safeDefault=0.0);
  double asDouble(size_t rec, size_t i, double safeDefault=0.0){return asDouble(i + dims[0]*rec, safeDefault);}
  double asDouble(size_t rec, size_t i, size_t j, double safeDefault=0.0){return asDouble(j + dims[1]*(i + dims[0]*rec), safeDefault);}
  double asDouble(size_t rec, size_t i, size_t j, size_t k, double safeDefault=0.0){return asDouble(k + dims[2]*(j + dims[1]*(i + dims[0]*rec)), safeDefault);}
  double asDouble(size_t rec, size_t i, size_t j, size_t k, size_t m, double safeDefault=0.0){return asDouble(m + dims[3]*(k + dims[2]*(j + dims[1]*(i + dims[0]*rec))), safeDefault);}
  int asInt(size_t posn=0, int safeDefault=0);
  DvString asStr(size_t posn=0, DvString safeDefault=DvString(""));
  DvString asText(size_t posn=0, DvString safeDefault=DvString(""));
  DvTime asTime(size_t posn=0, DvTime safeDefault=DvTime());
  
  DvObject_var convertToDbl();
  
  // operators
  
  DvString getTypeAsText();
  
  
  // special scaling
  DvObject &operator*=( double arg);
  DvObject &operator*=( int arg);
  DvObject &operator+=( double arg);
  DvObject &operator+=( int arg);
  DvObject &operator+=( DvEvent arg);
  DvObject &operator+=( DvString arg);
  DvObject &operator-=( double arg);
  DvObject &operator-=( int arg);
  DvObject &operator/=( double arg);
  DvObject &operator/=( int arg);
  
  // general arithmetic
  DvObject &operator+=( DvObject &arg);
  DvObject &operator-=( DvObject &arg);
  DvObject &operator*=( DvObject &arg);
  DvObject &operator/=( DvObject &arg);
  
  DvObject *operator+( DvObject &arg); 
  DvObject *operator-( DvObject &arg); 
  DvObject *operator*( DvObject &arg); 
  DvObject *operator/( DvObject &arg); 


  DvObject *multiplyInt( DvObject &arg); // internal use called by * if both ints


  // error flags
  void clearErr(){errList.clear();}
  void error(const char* flag){ 
  	cout << " ::: " << flag << endl;fflush(stdout);
  	errList.push_back( DvString(flag) );
  }
  DvString getError(size_t i=1000){
  	if(i == 1000) {
		if(errList.size() > 0) return errList[errList.size()-1]; // last error
		else return DvString();
	}
  	else if (errList.size() > i) return errList[i];
	else return DvString();
  }
  int nError(){
  	return errList.size();
  }
    
    void addErr(DvObject &obj) {for(size_t n=0; n<obj.nError(); n++) errList.push_back(obj.getError(n));}
    
  
  // xref methods
  // ------------

  
  DvNode *first_xref( ){
    if(xrefs.size() < 1) return NULL;
      
    return xrefs.first();
  }
  
  DvNode *find_xref( const char * xref_name ){
	return xrefs.find(xref_name);
  }
  
  DvNode *find_xref( DvString& xref_name ){
	return xrefs.find(xref_name);
  }

  DvObject_var get_xref( const char * xref_name ){

  	DvNode *xref = xrefs.find(xref_name);
  	if(xref) {
		return xref->obj();
	}
	else return new DvObject(); // empty object
  }
  
  DvObject_var get_xref(  DvString xref_name ){ 
  	DvNode *xref = xrefs.find(xref_name);
  	if(xref) return xref->obj();
	else return new DvObject(); // empty object
  }


  void delete_xref( const char * xref_name ){ 
	xrefs.remove(DvString(xref_name));
  }
  
  void delete_xref( DvString xref_name ){ 	
	xrefs.remove(xref_name);
  }

  bool xref_exists( const char * xref_name){
  	return (xrefs.find(xref_name) != NULL );
  }
  
  bool xref_exists( DvString xref_name){
  	return (xrefs.find(xref_name) != NULL );
  }

  void list_xref_names( vector<DvString>& sl ){
  	DvNode *xref = xrefs.first();
	while(xref){
		sl.push_back(xref->name());
		xref = xref->next;		
	}
  }

  bool xref_in_graph(DvObject &xref_obj){

	// used to prevent circular xrefs ever being created
	if(this->get_id() == xref_obj.get_id()) return true; 

	// check this is not in xref object list of xrefs
  	DvNode *xref = xref_obj.first_xref();
	while(xref){	
	
		if( this->xref_in_graph(*(xref->obj())) ) return true;
		xref = xref->next;		
	}
	
	return false;
  }

	
  void change_xref(  const char * xref_name, DvObject_var &xref_obj){
	DvString name(xref_name);
	change_xref(name, xref_obj);
  }
   
  void change_xref(  const char * xref_name, DvString &xref){
	DvObject_var xref_var = new DvObject(xref);
	DvString name(xref_name);
	change_xref(name, xref_var);
  }
  
  void change_xref( DvString & xref_name, DvString &xref){
	change_xref(xref_name.c_str(), xref);
  }
  
  void change_xref(  const char * xref_name, const char *xref){
    DvString xref_str(xref);
	DvObject_var xref_var = new DvObject(xref_str);
	DvString name(xref_name);
	change_xref(name, xref_var);
  }
  
  void change_xref(  const char * xref_name, int xref){
	DvObject_var xref_var = new DvObject(xref);
	DvString name(xref_name);
	change_xref(name, xref_var);
  }
   
  void change_xref(  const char * xref_name, double xref){
	DvObject_var xref_var = new DvObject(xref);
	DvString name(xref_name);
	change_xref(name, xref_var);
  }
   
  void change_xref(  DvString xref_name, DvObject_var &xref_obj)
  { 
	if(this->get_id() == xref_obj->get_id()) return; // forbidden to attach object to itself
	
    if(xref_in_graph(*xref_obj)) {
		// do not permit circular linkage
		DvObject_var copy = new DvObject(xref_obj, false); // copy without xrefs
		change_xref(xref_name, copy);
		return;
	}
	
    if( xref_exists(xref_name) ){
		DvNode *oldXref = find_xref(xref_name);
		oldXref->setObj(xref_obj);
	}
	else{
		xrefs.append(xref_name, xref_obj);
	}

  }
  
  void copy_xrefs_from(DvObject_var &from_obj){
  	DvNode *xref = from_obj->first_xref();
	while(xref){	
		
		change_xref(xref->name(), xref->obj());
		xref = xref->next;
	}
  }

  void copy_xrefs_from(DvObject &from_obj){
  	DvNode *xref = from_obj.first_xref();
	while(xref){	
		change_xref(xref->name(), xref->obj());
		xref = xref->next;
	}
  }

  void copy_xref_from( const char *xref_name, DvObject &from_obj){

	copy_xref_from(DvString(xref_name), from_obj);
  }
  
  void copy_xref_from( DvString xref_name, DvObject &from_obj){
  	DvObject_var xref = from_obj.get_xref(xref_name);
	
	if( xref->is_ok()) change_xref(xref_name, xref);
  }

  void copy_xref_from( const char *xref_name, DvObject_var &from_obj){
	copy_xref_from(DvString(xref_name), from_obj);
  }
  
  void copy_xref_from( DvString xref_name, DvObject_var &from_obj){
  	DvObject_var xref = from_obj->get_xref(xref_name);
	
	if( xref->is_ok()) change_xref(xref_name, xref);
  }


  // General Attribute Utilities
  
  DvString getXrefText(DvString &name);
  DvString getXrefText(const char *name);
  int get_iFILL();
  double get_dFILL();
  DvObject_var getDep0();
  DvObject_var getTimeTags();
  bool hasTimeTags();
  DvEvent getTimeRange();
  bool getDataRange(double &min, double &max);
  bool getDataRange(size_t recFrom, size_t recTo, double &min, double &max);
   
  // Join Utilities
  bool isRegular(double &spacing);
  bool isJoined(DvObject &arg);
  bool same(DvObject &arg);
  DvObject_var Join(DvObject &Target, bool withXrefs=true);
  DvObject_var Join(DvObject_var &Target, bool withXrefs=true);
  DvObject_var linearJoin(DvObject &Target, bool withXrefs=true, DvMask *Gmsk=NULL);
  DvObject_var boxcarJoin(DvObject &Target, bool withXrefs=true, DvMask *Gmsk=NULL);
  DvObject_var nnJoin(DvObject &Target, bool withXrefs=true,  DvMask *Gmsk=NULL);
  double get_spacing();
  void joinXrefs(DvObject &obj, DvObject &Target);
  DvObject_var interpAt(DvTime &target);
  valarray <double> getCentres();
  void applyLinearFill(DvMask &msk);


  // Units Utilities
  bool sameBaseSI(DvObject &arg);
  bool sameBaseSI(const char * SIstr);
  bool sameUnits(DvObject &arg);
  bool sameUnits(const char * SIstr, int i=0);
  DvString getSIC(size_t i=0);
  int getSICdim();
  double convFactor();
  double convFactorToBaseSI();
  double convFactorFrom(DvString &SIC);
  DvString getBaseSI(int i=0);
  DvString getUnitsProduct(DvObject &arg);
  DvString getUnitsRatio(DvObject &arg);
  DvString getUnitsPower(double pow);
  DvString getUnitsInverse();
  DvString getSICProduct(DvObject &arg);
  DvString getSICInverse();
  DvString getSICRatio(DvObject &arg);
  DvString getSICPower(double pow);
  bool hasUnits();
  bool isDeg();
  bool isRad();
  bool isVectDeg(); // checks all components
  bool isVectRad();
  bool isAngle();
  bool isPhiAngle(bool isDeg, int index, double &minVal, int &kStart, int &kEnd);
  DvString angleUnitStr();
  DvObject_var angleMod(int angMax, DvString rad_deg);
  DvString getCommonText();


  // Units conversion
  DvObject_var changeUnitsTo(DvString SIC, DvString Units=DvString(""));

  // Coordinate System Utilities
  bool sameFrame(DvObject &arg);
  bool sameFrame(DvString &arg);
  void setFrameAttr(DvString frameAttr);
  void setFrameAttr(DvString frm, DvString rep, int order);
  void setAttrsProduct(DvObject &arg);
  void setToFrameAttr(const char *frame);
  DvString getFrameAttr();
  DvString getFrame();
  DvString getRep();
  int getOrder();
  bool isThreeVector();
  bool isVectorXYZ();
  void ensureVectorXYZ();
  
  // Matrix Operations
  
  // slice access to valarray (no bounds protection)
  std::slice_array <double> getDSlice(slice sl){ return Ddata[sl];}
  std::slice_array <int> getISlice(slice sl){ return Idata[sl];}
  std::slice_array <DvString> getSSlice(slice sl){ return Sdata[sl];}
  std::slice_array <DvTime> getTSlice(slice sl){ return Tdata[sl];}
  std::slice_array <DvEvent> getESlice(slice sl){ return Edata[sl];}
  
  void setSlice(slice sl, std::slice_array <double> val){ Ddata[sl] = val;}
  void setSlice(slice sl, std::slice_array <int> val){ Idata[sl] = val;}
  void setSlice(slice sl, std::slice_array <DvString> val){ Sdata[sl] = val;}
  void setSlice(slice sl, std::slice_array <DvTime> val){ Tdata[sl] = val;}
  void setSlice(slice sl, std::slice_array <DvEvent> val){ Edata[sl] = val;}
  
  void setSliceDim(size_t d, size_t at, std::slice_array <double> val);
  void setSliceDim(size_t d, size_t at, std::slice_array <int> val);
  void setSliceDim(size_t d, size_t at, std::slice_array <DvString> val);
  void setSliceDim(size_t d, size_t at, std::slice_array <DvTime> val);
  void setSliceDim(size_t d, size_t at, std::slice_array <DvEvent> val);

  valarray <double> detValarray();
  std::slice_array <double> elementValarray(size_t i);

  DvObject_var matrixInverse();
  DvObject_var det();
  DvObject_var adjoint();
  DvObject_var element(size_t i);
  DvObject_var transpose(bool withXref=true);
  DvObject_var subMatrix(size_t i, size_t j, bool withXref=true);
  DvObject_var subDimension(size_t i, bool withXref=true);
  DvObject_var trace();
  DvObject_var subsetDim(size_t d, size_t from, size_t to);
  DvObject_var sliceDim(size_t d, size_t at);
  DvObject_var sumDim(size_t d, size_t from, size_t to, bool stripBad=true);
  DvObject_var averageDim(size_t d, size_t from, size_t to, bool stripBad=true);
  DvObject_var unwrapDim(size_t d, size_t from, size_t to);
  DvObject_var reverseDim(size_t d, size_t from, size_t to);
  void unwrapXrefs(size_t oldSeqLen,  size_t nElem);
  void adjustDeltas(size_t nElem);
  DvObject_var outerProduct(DvObject_var &obj);
 
  void cycleMod(int nstart, int n, int modn);
  DvObject_var correctColToRow();
  
  // Components and subset
  vector <DvString> available(DvString constraint=DV_NONE);
  DvObject_var getObjComponent(DvString comp);
  DvObject_var getVectorComponent(DvString comp);
  DvObject_var getTimeComponent(DvString comp);
  DvObject_var getIntervalComponent(DvString comp);
  DvObject_var getArrayComponent(DvString comp);
  DvObject_var getTimeSubset(DvString comp);
  DvObject_var getTimeAve(DvString comp);
  DvEvent getTimeInterval(DvString comp=DvString(""));
  DvObject_var sum();
  DvObject_var sumRecEvery(size_t stride);
  DvObject_var aveRecEvery(size_t stride);
  DvObject_var sliceRecEvery(size_t stride);
  
  DvObject_var getVecFromRLP();
  DvObject_var getVecFromRTP();


  // algorithms
  DvObject_var minusStart();
  DvObject_var integrate(DvString gapMethod, DvString integralMethod);
  DvObject_var differentiate(DvString derivMethod=DV_DERIV_3PT, double gapSize=-1.);
  DvObject_var regression();
  DvObject_var regression(DvObject_var &xobj);
  DvObject_var mergeWith(DvObject_var &dobj);
  DvObject_var mergeByMask(size_t nMask, vector <int> &pairMask, DvObject_var &dobj);
  void stripDuplicates(double tolerance);
  DvObject_var cleanObject();

  DvObject_var makeMonotonic();
  DvObject_var applyOrder(valarray <size_t> &order);
  void ensureSIC();

  // Dimension Utilities
  bool okDims(DvObject &arg);
  bool sameDims(DvObject &arg);
  bool conformalDims(DvObject &arg);
  bool squareMat();

  DvString stringArray(size_t n=0){
    
    DvString strArr("");
      
    if(arraySize() == 1) strArr = this->asStr(n);
    
    else if(this->nDims() == 1){
        for(size_t i=0; i<this->dims[0]; i++){
            strArr += this->asStr(n*arraySize()+i);
            if(i < this->dims[0] - 1)  strArr += ", ";    
        }
    }
      
    else if(this->nDims() == 2){
        for(size_t i=0; i<this->dims[0]; i++){
            for(size_t j=0; j<this->dims[1]; j++){
                strArr += this->asStr(n*arraySize() +i*this->dims[1] + j);
                if(j < this->dims[1] - 1)  strArr += ", ";    
            }
            strArr += "\n";
        }
    }
      
    else  if(this->nDims() == 3){
        for(size_t i=0; i<this->dims[0]; i++){
            for(size_t j=0; j<this->dims[1]; j++){
                for(size_t k=0; k<this->dims[2]; k++){
                    strArr += this->asStr(((n*dims[0] + i)*dims[1] + j)*dims[2] + k);
                    if(k < this->dims[2] - 1)  strArr += ", ";   
                }
                strArr += "\n";
            }
            strArr += "\n";
        }
    }
      
    else{
        for(size_t i=0; i<this->dims[0]; i++){
            for(size_t j=0; j<this->dims[1]; j++){
                for(size_t k=0; k<this->dims[2]; k++){
                    for(size_t m=0; m<this->dims[3]; m++){
                        strArr += this->asStr((((n*dims[0] + i)*dims[1] + j)*dims[2] + k)*dims[3] + m);
                        if(m < this->dims[3] - 1)  strArr += ", ";   
                    }
                    strArr += "\n";
                }
                strArr += "\n";
            }
            strArr += "\n";
        }
    }
      
    strArr += "\n";

    return strArr;
  }
    
  // operators
  DvObject_var sqrt();
  void      sqrtThis();
  
  DvObject_var abs();
  void      absThis();
  
  DvObject_var log();
  bool      logThis();
  
  DvObject_var log10();
  bool      log10This();
  
  DvObject_var inverse();
  void      inverseThis();

  DvObject_var chgSign();
  void      chgSignThis();
  
  DvObject_var power(double p);
  void      powerThis(double p);
  DvObject_var power(int p);
  void      powerThis(int p);

  DvObject_var remainder(double d);
  void      remainderThis(double d);
  DvObject_var remainder(DvObject_var &obj);
  bool      remainderThis(DvObject_var &obj);

  DvObject_var exp();
  bool      expThis();
  
  double 	max();
  double 	min();

  DvObject_var toDeg();
  void 		toDegThis();

  DvObject_var toRad();
  void 		toRadThis();

  
  DvObject_var XYtoRP();
  void 		XYtoRPThis();
  DvObject_var XYtoRP_deg();
  void 		XYtoRP_degThis();
  DvObject_var RPtoXY();
  void 		RPtoXYThis();

  DvObject_var XYZtoRPZ();
  void 		XYZtoRPZThis();
  DvObject_var XYZtoRPZ_deg();
  void 		XYZtoRPZ_degThis();
  DvObject_var RPZtoXYZ();
  void 		RPZtoXYZThis();
  
  DvObject_var XYZtoRTP();
  void 		XYZtoRTPThis();
  DvObject_var XYZtoRTP_deg();
  void 		XYZtoRTP_degThis();
  DvObject_var RTPtoXYZ();
  void 		RTPtoXYZThis();
  
  DvObject_var XYZtoRLP();
  void 		XYZtoRLPThis();
  DvObject_var XYZtoRLP_deg();
  void 		XYZtoRLP_degThis();
  DvObject_var RLPtoXYZ();
  void 		RLPtoXYZThis();
  
  DvObject_var sin();
  void      sinThis();
  
  DvObject_var cos();
  void      cosThis();
  
  DvObject_var tan();
  void      tanThis();
  
  DvObject_var cot();
  void      cotThis();
  
  DvObject_var sinh();
  void      sinhThis();
  
  DvObject_var cosh();
  void      coshThis();
  
  DvObject_var tanh();
  void      tanhThis();
  
  DvObject_var coth();
  void      cothThis();

  DvObject_var asin(bool indeg=false);
  void      asinThis(bool indeg=false);

  DvObject_var acos(bool indeg=false);
  void      acosThis(bool indeg=false);

  DvObject_var atan(bool indeg=false);
  void      atanThis(bool indeg=false);

  DvObject_var atan2(DvObject_var &x, bool indeg=false);
  void      atan2This(DvObject_var &x, bool indeg=false);

  DvObject_var acot(bool indeg=false);
  void      acotThis(bool indeg=false);
 
  DvObject_var asinh();
  void      asinhThis();
  
  DvObject_var acosh();
  void      acoshThis();
  
  DvObject_var atanh();
  void      atanhThis();
  
  DvObject_var acoth();
  void      acothThis();


  DvObject_var dot(DvObject_var& obj);
  DvObject_var vec(DvObject_var& obj);
  DvObject_var normalize();
  
  DvObject_var multElements(DvObject_var& obj);
  DvObject_var divideElements(DvObject_var& obj);
		
  void      toXML(std::ofstream &ofp);
  void      fromXML(char *xml);
  

};



//-----------------------------------------------------------------------------
// class DvJoinList
//-----------------------------------------------------------------------------

class DvJoinList : public DvList {
		
	public:

		DvString MultiJoin( DvObject_var     &target,
							DvJoinList       &retList); // DvString return is status message
							
};

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ class DvMatrixEig
// Class for calculating eigenvalues/egigenvectors for a rank 2 square matrix
// -----------------------------------------------------------------------------
	void d_eigen(valarray<double> &a, valarray<double> &v, size_t n);
	void c_eigen(std::vector<std::complex<double> > &array, std::vector<double> &eigvalues, size_t size);

	class DvMatrixEig
	{
		protected:
			size_t _n;               // no. rows/cols (also used as error flag = 0)
			DvObject_var _eigvals;   // _n elt 
			DvObject_var _eigvecs;   // _n x _n elt 

		public:
     // --- Constructors
     // Empty constructor
			DvMatrixEig(){;}
			
     // Construct eigenvalues/eigenvectors from Matrix
			DvMatrixEig(DvObject_var &m)
			{
				// verify that matrix argument is rank 2
				if (m.is_nil() || m->nDims() !=2 ) _n=0;
				else _n = m->Dims()[0];
				
				if(m->seqSize() > 1) _n=0;
				
				//  verify that matrix is symmetric
				if (_n != 0 && m->Dims()[1] == _n)
				{
					vector<size_t> dims;
					dims.push_back(_n);
					_eigvals = new DvObject(0.0, dims); // dim _n
					
					dims.push_back(_n);
					_eigvecs = new DvObject(m,false); // dim _n x _n copy of m initially

					d_eigen(_eigvecs->d_valarray(), _eigvals->d_valarray(), _n);
					
					_eigvecs = _eigvecs->transpose(false); // convert to eigenvectors as rows (no metadata)
					
				}
			}

     // --- Misc functions
			DvObject_var  eigenvalues()
			{
				return _eigvals;
			}

			DvObject_var  eigenvectors()
			{
				return _eigvecs;
			}

	};
	


class DvMinVariance
{
    DvObject_var _eigvals;       // eigenvalues
    DvObject_var _eigvecs;       // eigenvectors 
    DvObject_var _eigqual;       // Matrix Dim(3) eigenvalue quality
    DvObject_var _covariance;    // covariance matrix
    DvObject_var P;              // projection matrix P (see ISSI book, p 195)
    double _orgqual;             // orthogonality quality
    
public:
    DvMinVariance(const DvObject_var &ms, const DvObject_var &constraint, const int Qsel=0)	
	{
		_covariance = new DvObject((double)0.0, DIM33);
		P           = new DvObject((double)0.0, DIM33);
		
		if (ms.is_nil()){
			_covariance->error("MV: nil input sequence");
			return;
		}
		if (ms->nDims() != 1){
			_covariance->error("MV: matrix rank not 1");
			return;
		}
		if (ms->Dims()[0] != 3){
			_covariance->error("MV: not a vector");
			return;
		}

		bool constrained =false;
		if(constraint.is_ok()) constrained = true;
		
		// compute the variance matrix
		int sz = ms->seqSize();
		DvObject_var sum_ms  = ms->sum();
		double sumxx = 0.0;
		double sumxy = 0.0;
		double sumxz = 0.0;
		double sumyy = 0.0;
		double sumyz = 0.0;
		double sumzz = 0.0;
		for(int i=0; i<sz; i++){
			sumxx += ms->dbl(i,0) * ms->dbl(i,0);
			sumxy += ms->dbl(i,0) * ms->dbl(i,1);
			sumxz += ms->dbl(i,0) * ms->dbl(i,2);
			sumyy += ms->dbl(i,1) * ms->dbl(i,1);
			sumyz += ms->dbl(i,1) * ms->dbl(i,2);
			sumzz += ms->dbl(i,2) * ms->dbl(i,2);
		}
		

		switch (Qsel) {
		case 2: // Weimer type covariance matrix
			//printf("MVA matrix with Weimer type covariance matrix\n");
			_covariance->dbl(0,0,0) = (sumxx - sum_ms->dbl(0) * sum_ms->dbl(0)) / sz;
			_covariance->dbl(0,0,1) = (sumxy - sum_ms->dbl(0) * sum_ms->dbl(1)) / sz;
			_covariance->dbl(0,0,2) = (sumxz - sum_ms->dbl(0) * sum_ms->dbl(2)) / sz;
			_covariance->dbl(0,1,1) = (sumyy - sum_ms->dbl(1) * sum_ms->dbl(1)) / sz;
			_covariance->dbl(0,1,2) = (sumyz - sum_ms->dbl(1) * sum_ms->dbl(2)) / sz;
			_covariance->dbl(0,2,2) = (sumzz - sum_ms->dbl(2) * sum_ms->dbl(2)) / sz;
			constrained = false;
			break;
		case 3: // Siscoe type covariance matrix
			//printf("MVA matrix with Siscoe type covariance matrix\n");
			_covariance->dbl(0,0,0) = sumxx / sz;
			_covariance->dbl(0,0,1) = sumxy / sz;
			_covariance->dbl(0,0,2) = sumxz / sz;
			_covariance->dbl(0,1,1) = sumyy / sz;
			_covariance->dbl(0,1,2) = sumyz / sz;
			_covariance->dbl(0,2,2) = sumzz / sz;
			constrained = false;
			break;
		case 4: // extreme Weimer (for N -> infinity) NB! for record only - gives nonsense
			//printf("MVA matrix with extreme Weimer type covariance matrix\n");
			_covariance->dbl(0,0,0) = -sum_ms->dbl(0) * sum_ms->dbl(0) / sz;
			_covariance->dbl(0,0,1) = -sum_ms->dbl(0) * sum_ms->dbl(1) / sz;
			_covariance->dbl(0,0,2) = -sum_ms->dbl(0) * sum_ms->dbl(2) / sz;
			_covariance->dbl(0,1,1) = -sum_ms->dbl(1) * sum_ms->dbl(1) / sz;
			_covariance->dbl(0,1,2) = -sum_ms->dbl(1) * sum_ms->dbl(2) / sz;
			_covariance->dbl(0,2,2) = -sum_ms->dbl(2) * sum_ms->dbl(2) / sz;
			constrained = false;
			break;
		default:
			//printf("MVA matrix with classic covariance matrix\n");
			_covariance->dbl(0,0,0) = (sumxx - sum_ms->dbl(0) * sum_ms->dbl(0) / sz) / sz;
			_covariance->dbl(0,0,1) = (sumxy - sum_ms->dbl(0) * sum_ms->dbl(1) / sz) / sz;
			_covariance->dbl(0,0,2) = (sumxz - sum_ms->dbl(0) * sum_ms->dbl(2) / sz) / sz;
			_covariance->dbl(0,1,1) = (sumyy - sum_ms->dbl(1) * sum_ms->dbl(1) / sz) / sz;
			_covariance->dbl(0,1,2) = (sumyz - sum_ms->dbl(1) * sum_ms->dbl(2) / sz) / sz;
			_covariance->dbl(0,2,2) = (sumzz - sum_ms->dbl(2) * sum_ms->dbl(2) / sz) / sz;
			break;
		}

		_covariance->dbl(0,1,0) = _covariance->dbl(0,0,1);
		_covariance->dbl(0,2,0) = _covariance->dbl(0,0,2);
		_covariance->dbl(0,2,1) = _covariance->dbl(0,1,2);

		//
		// Constrained MVA : find eigenvec of [P # C # P] * n = lamda * n
		// instead of [C] * n = lambda * n.
		//
		//   0 - standard MVA
		//   1 - constrain so that <ts2> * dot n == 0 (see p 209 in blue ISSI book)
		//
		// SEH, 14 May 2002
		//      24 Aug 2002 (constrain against ts2)
		//

		if (constrained) {
		
			// constrained MV
			DvObject_var sum_c = constraint->sum();
			double sz2 = sum_c->dbl(0) * sum_c->dbl(0) + sum_c->dbl(1) * sum_c->dbl(1) + sum_c->dbl(2) * sum_c->dbl(2);

			P->dbl(0,0,0) = 1.0 - (sum_c->dbl(0) * sum_c->dbl(0)) / sz2;
			P->dbl(0,0,1) = -(sum_c->dbl(0) * sum_c->dbl(1)) / sz2;
			P->dbl(0,0,2) = -(sum_c->dbl(0) * sum_c->dbl(2)) / sz2;

			P->dbl(0,1,0) = -(sum_c->dbl(1) * sum_c->dbl(0)) / sz2;
			P->dbl(0,1,1) = 1.0 - (sum_c->dbl(1) * sum_c->dbl(1)) / sz2;
			P->dbl(0,1,2) = -(sum_c->dbl(1) * sum_c->dbl(2)) / sz2;

			P->dbl(0,2,0) = -(sum_c->dbl(2) * sum_c->dbl(0)) / sz2;
			P->dbl(0,2,1) = -(sum_c->dbl(2) * sum_c->dbl(1)) / sz2;
			P->dbl(0,2,2) = 1.0 - (sum_c->dbl(2) * sum_c->dbl(2)) / sz2;

			DvObject_var temp = _covariance * P;
			_covariance = P * temp;

			// force symmetric covariance matrix again
			_covariance->dbl(0,1,0) = _covariance->dbl(0,0,1);
			_covariance->dbl(0,2,0) = _covariance->dbl(0,0,2);
			_covariance->dbl(0,2,1) = _covariance->dbl(0,1,2);
		}

		// get eigenvectors and values
		DvMatrixEig eig(_covariance);
		_eigvecs = eig.eigenvectors();
		_eigvals = eig.eigenvalues();

		// If constrained : eigval 0 == 0. Permute eigenvecs and eigenvals.
		// For e.g., finding boundary normals, one has to look at the
		// eigenvector corresponding to the intermediate eigenvalue
		// to get a proxy for the normal.
		// See blue ISSI book, ch8 for more info.

		if (constrained) {
			_eigvals->dbl(0) = 0.0;
		}

		// check the positive sense of the eigenvector associated with
		// the minimum eigenvalue and if necessary reverse the signs
		// of all three components i.e. vint x vmax = vector perpendicular
		// to both vint and vmax. If vmin is perpendicular to vint
		// and vmax (i.e. system is ortogonal), then (vinx x vmax).vmin=1
		// If < 0 then not a right-handed system, so reverse signs of vmin

		DvObject_var eigV0 = _eigvecs->sliceDim(0,0);
		DvObject_var eigV1 = _eigvecs->sliceDim(0,1);
		DvObject_var eigV2 = _eigvecs->sliceDim(0,2);
		
		DvObject_var rh = eigV1->vec(eigV2);
		rh = rh->dot(eigV0);
		double rh_factor = rh->dbl();
		if (rh_factor < 0.) *eigV0 *= -1.0;
		_orgqual = fabs(fabs(rh_factor) - 1.0);

		// ensure that the sign of the average minimum variance component
		// is positive by computing average of sum of B-field data.min
		// eigenvector. If less than zero negate both vmin and vint (so
		// preserving right-handed system)
		// double avg_minvar_comp=(ms*_eigvecs[0]).sum();
		//        if(avg_minvar_comp<0)

		// SEH 9 May 2005 : ennforce positive x-dir for minvar direction


		if (constrained && eigV1->dbl(0) < 0.0) {
			*eigV0 *= -1.0;
			*eigV1 *= -1.0;
		}
		else if ( !constrained && eigV0->dbl(0) < 0.0) {
			*eigV0 *= -1.0;
			*eigV1 *= -1.0;
		}
		
		// replace values in _eigvecs 2D array
		slice ds(0, eigV0->totalSize(), 1); // array_slice of all data in vector
		
		_eigvecs->setSliceDim(0, 0, eigV0->getDSlice(ds));
		_eigvecs->setSliceDim(0, 1, eigV1->getDSlice(ds));
		_eigvecs->setSliceDim(0, 2, eigV2->getDSlice(ds));
		
        _eigqual  = new DvObject(_eigvals);
        *_eigqual /= _eigvals->dbl(2);

	}

    // methods for returning members and calculating variance
    DvObject_var covariance() 
    { return _covariance; }

    DvObject_var projection() 
    { return P; }
    
    DvObject_var eigenvalues() 
    { return _eigvals;}
     
    DvObject_var eigenvectors()
    { return _eigvecs;}   
    
    DvObject_var eigenquality() 
    { return _eigqual;}
    
    double orthogonality_quality()
    { return _orgqual;}    // orthogonality quality     

        
};
 
// ------------------ internal utilities -----------------

valarray <double> atan_yx(valarray <double> &y, valarray <double> &x);
double atan_yx( double y,  double x);
bool isFill(double val, double fill);
bool notFill(double val, double fill);

// ----------------  version utility calls --------------

inline DvString getDVOS_versionStr(){
	return DVOS_VERSION_STR;
}
inline DvString getDVOS_versionDate(){
	return DVOS_VERSION_DATE;
}
inline int getDVOS_version_major(){
	DvString ver = DVOS_VERSION_STR;
	int val = ver.before('.').toInt();
	return val;
}
inline int getDVOS_version_minor(){
	DvString ver = DVOS_VERSION_STR;
	ver = ver.after('.');
	int val = ver.before('.').toInt();
	return val;
}
inline int getDVOS_version_sub(){
	DvString ver = DVOS_VERSION_STR;
	ver = ver.after('.');
	ver = ver.after('.');
	int val = ver.toInt();
	return val;
}

#endif 

