//
//  DvObj.h
//
//  Dvos -  Data variable object classes
//
//  Tony Allen
//  Nov 2012
//

#include <iostream>

#include "DvObj.h"
#include "DvRecord.h"
	
	
//   -------- DvNode class ---------



	DvNode::DvNode(const char *name, DvObject_var &obj):_obj(obj),_name(name){
		next = NULL;
	}
	
	DvNode::DvNode(string name, DvObject_var &obj):_obj(obj),_name(name.c_str()){
		next = NULL;
	}
	
	DvNode::DvNode(DvString name, DvObject_var &obj):_obj(obj),_name(name.c_str()){
		next = NULL;
	}

//   --------- DvList Class  --------

DvNode *DvList::makeNode(){
		DvObject_var nilObj = new DvObject();
		DvNode *newNode = new DvNode("", nilObj);
        DvNode *lastNode = last();
        if(lastNode) lastNode->next = newNode;
        else nodes = newNode;
    
		return newNode;
	}



// ------------

int DvObject::ID_counter = 0;

// constructors

DvObject::DvObject() 
{
  // empty object used as place holder in empty DvObject_var
  object_id=-1; // id of all empty objects
  ref_count=0; // used by Dvar to know when safe to delete
  
  seqLen=1;
  secRes = 3;
  dims.push_back(1);
  Ddata.resize(1, 0.0); // ensure something exists so safe to de-reference

}

DvObject::DvObject(double value, size_t nRecs) 
{
  object_id=++ID_counter; // unique object id
  ref_count=0; // used by Dvar to know when safe to delete
 
  seqLen=nRecs;
  Ddata.resize(nRecs, value);
  secRes = 3;
  dims.clear();
  dims.push_back(1);

}
DvObject::DvObject(int value, size_t nRecs) 
{
  object_id=++ID_counter; // unique object id
  ref_count=0; // used by Dvar to know when safe to delete
 
  seqLen=nRecs;
  Idata.resize(nRecs, value);
  secRes = 3;
  dims.clear();
  dims.push_back(1);
}
DvObject::DvObject(DvString value, size_t nRecs) 
{
  object_id=++ID_counter; // unique object id
  ref_count=0; // used by Dvar to know when safe to delete
 
  seqLen=nRecs;
  Sdata.resize(nRecs, value);
  secRes = 3;
  dims.clear();
  dims.push_back(1);
}
DvObject::DvObject(const char * value, size_t nRecs) 
{
  object_id=++ID_counter; // unique object id
  ref_count=0; // used by Dvar to know when safe to delete
 
  seqLen=nRecs;
  Sdata.resize(nRecs, DvString(value));
  secRes = 3;
  dims.clear();
  dims.push_back(1);
}
DvObject::DvObject(DvTime value, size_t nRecs) 
{
  object_id=++ID_counter; // unique object id
  ref_count=0; // used by Dvar to know when safe to delete
 
  seqLen=nRecs;
  Tdata.resize(nRecs, value);
  secRes = value.getSecResolution();
  dims.clear();
  dims.push_back(1);
}
DvObject::DvObject(DvEvent value, size_t nRecs) 
{
  object_id=++ID_counter; // unique object id
  ref_count=0; // used by Dvar to know when safe to delete
 
  seqLen=nRecs;
  Edata.resize(nRecs, value);
  secRes = value.start().getSecResolution();
  dims.clear();
  dims.push_back(1);
}

// construct sequence of arrays and fill with value
// allows elements to be changed one at a time later 
DvObject::DvObject(double value, const vector <size_t> &dims_, size_t nRecs) 
{
  object_id=++ID_counter; // unique object id
  ref_count=0; // used by Dvar to know when safe to delete

  seqLen = nRecs;
  dims = dims_;
  
  size_t len = nRecs;
  for(size_t i=0; i<dims.size(); i++) len *= dims[i];
  Ddata.resize(len, value);
  secRes = 3;
}
DvObject::DvObject(int value, const vector <size_t> &dims_, size_t nRecs) 
{
  object_id=++ID_counter; // unique object id
  ref_count=0; // used by Dvar to know when safe to delete

  seqLen = nRecs;
  dims = dims_;
  
  size_t len = nRecs;
  for(size_t i=0; i<dims.size(); i++) len *= dims[i];
  Idata.resize(len, value);
  secRes = 3;
}
DvObject::DvObject(const char *value, const vector <size_t> &dims_, size_t nRecs)
{
    object_id=++ID_counter; // unique object id
    ref_count=0; // used by Dvar to know when safe to delete
    
    seqLen = nRecs;
    dims = dims_;
    
    DvString strVal = value;
    size_t len = nRecs;
    for(size_t i=0; i<dims.size(); i++) len *= dims[i];
    Sdata.resize(len, strVal);
    secRes = 3;
}
DvObject::DvObject(DvString value, const vector <size_t> &dims_, size_t nRecs)
{
    object_id=++ID_counter; // unique object id
    ref_count=0; // used by Dvar to know when safe to delete
    
    seqLen = nRecs;
    dims = dims_;
    
    size_t len = nRecs;
    for(size_t i=0; i<dims.size(); i++) len *= dims[i];
    Sdata.resize(len, value);
    secRes = 3;
}
DvObject::DvObject(DvTime value, const vector <size_t> &dims_, size_t nRecs)
{
  object_id=++ID_counter; // unique object id
  ref_count=0; // used by Dvar to know when safe to delete

  seqLen = nRecs;
  dims = dims_;
  
  size_t len = nRecs;
  for(size_t i=0; i<dims.size(); i++) len *= dims[i];
  Tdata.resize(len, value);
  secRes = value.getSecResolution();
}
DvObject::DvObject(DvEvent value, const vector <size_t> &dims_, size_t nRecs) 
{
  object_id=++ID_counter; // unique object id
  ref_count=0; // used by Dvar to know when safe to delete

  seqLen = nRecs;
  dims = dims_;
  
  size_t len = nRecs;
  for(size_t i=0; i<dims.size(); i++) len *= dims[i];
  Edata.resize(len, value);
  secRes = value.start().getSecResolution();
}

// Copy Constructors 

DvObject::DvObject(DvObject &inObj, bool withXrefs) 
{
  object_id=++ID_counter; // unique object id
  ref_count=0; // used by Dvar to know when safe to delete

  seqLen = inObj.seqLen;
  dims = inObj.dims;
  
  if (withXrefs) copy_xrefs_from(inObj);

  Ddata.resize(inObj.Ddata.size());
  Ddata = inObj.Ddata;
  Idata.resize(inObj.Idata.size());
  Idata = inObj.Idata;
  Sdata.resize(inObj.Sdata.size());
  Sdata = inObj.Sdata;
  Tdata.resize(inObj.Tdata.size());
  Tdata = inObj.Tdata;
  Edata.resize(inObj.Edata.size());
  Edata = inObj.Edata;
  secRes = inObj.secRes;
    
  this->addErr(inObj);
    

}

DvObject::DvObject(DvObject_var &inObj, bool withXrefs) 
{
  object_id=++ID_counter; // unique object id
  ref_count=0; // used by Dvar to know when safe to delete

  seqLen = inObj->seqLen;
  dims = inObj->dims;
  
  if (withXrefs)  copy_xrefs_from(inObj);
  
  Ddata.resize(inObj->Ddata.size());
  Ddata = inObj->Ddata;
  Idata.resize(inObj->Idata.size());
  Idata = inObj->Idata;
  Sdata.resize(inObj->Sdata.size());
  Sdata = inObj->Sdata;
  Tdata.resize(inObj->Tdata.size());
  Tdata = inObj->Tdata;
  Edata.resize(inObj->Edata.size());
  Edata = inObj->Edata;
  secRes = inObj->secRes;
  
    this->addErr(*inObj);
}

DvObject::DvObject(DvObject *inObj, bool withXrefs) 
{

  object_id=++ID_counter; // unique object id
  ref_count=0; // used by Dvar to know when safe to delete

  seqLen = inObj->seqLen;
  dims = inObj->dims;
  
  if (withXrefs) copy_xrefs_from(*inObj);
  
  Ddata.resize(inObj->Ddata.size());
  Ddata = inObj->Ddata;
  Idata.resize(inObj->Idata.size());
  Idata = inObj->Idata;
  Sdata.resize(inObj->Sdata.size());
  Sdata = inObj->Sdata;
  Tdata.resize(inObj->Tdata.size());
  Tdata = inObj->Tdata;
  Edata.resize(inObj->Edata.size());
  Edata = inObj->Edata;
  secRes = inObj->secRes;
  
  this->addErr(*inObj);

}


DvObject::DvObject(DvObject_var &inObj, size_t nRecs, bool withXrefs) 
{
  // copy constructor repeating to nRecs
  
  object_id=++ID_counter; // unique object id
  ref_count=0; // used by Dvar to know when safe to delete
 
  seqLen = nRecs;
  dims = inObj->dims;
  
  if(withXrefs) copy_xrefs_from(*inObj);
  
  size_t arrayLen = inObj->arraySize();
  size_t len = arrayLen*nRecs;
  if(inObj->is_dbl()){
  	Ddata.resize(len);
  	for(size_t i=0; i<nRecs; i++)
  		for(size_t j=0; j<arrayLen; j++) 
			Ddata[i*arrayLen + j] = inObj->Ddata[j];
  }
  else if(inObj->is_int()){  
  	Idata.resize(len);
  	for(size_t i=0; i<nRecs; i++)
  		for(size_t j=0; j<arrayLen; j++) 
			Idata[i*arrayLen + j] = inObj->Idata[j];
  }
  else if(inObj->is_str()){  
  	Sdata.resize(len);
  	for(size_t i=0; i<nRecs; i++)
  		for(size_t j=0; j<arrayLen; j++)
  			Sdata[i*arrayLen + j] = inObj->Sdata[j];
  }
  else if(inObj->is_time()){ 
    Tdata.resize(len);
  	for(size_t i=0; i<nRecs; i++)
  		for(size_t j=0; j<arrayLen; j++)
  			Tdata[i*arrayLen + j] = inObj->Tdata[j];
  }
  else if(inObj->is_event()){ 
	Edata.resize(len);
  	for(size_t i=0; i<nRecs; i++)
  		for(size_t j=0; j<arrayLen; j++)
  			Edata[i*arrayLen + j] = inObj->Edata[j];
  }
  secRes = inObj->secRes;

  this->addErr(*inObj);

}

DvObject::DvObject(DvObject *inObj, size_t nRecs, bool withXrefs) 
{
  // copy constructor repeating to nRecs
  
  object_id=++ID_counter; // unique object id
  ref_count=0; // used by Dvar to know when safe to delete
 
  seqLen = nRecs;
  dims = inObj->dims;
  
  if(withXrefs) copy_xrefs_from(*inObj);
  
  size_t arrayLen = inObj->arraySize();
  size_t len = arrayLen*nRecs;
  if(inObj->is_dbl()){
  	Ddata.resize(len);
  	for(size_t i=0; i<nRecs; i++)
  		for(size_t j=0; j<arrayLen; j++) 
			Ddata[i*arrayLen + j] = inObj->Ddata[j];
  }
  else if(inObj->is_int()){  
  	Idata.resize(len);
  	for(size_t i=0; i<nRecs; i++)
  		for(size_t j=0; j<arrayLen; j++) 
			Idata[i*arrayLen + j] = inObj->Idata[j];
  }
  else if(inObj->is_str()){  
  	Sdata.resize(len);
  	for(size_t i=0; i<nRecs; i++)
  		for(size_t j=0; j<arrayLen; j++)
  			Sdata[i*arrayLen + j] = inObj->Sdata[j];
  }
  else if(inObj->is_time()){ 
    Tdata.resize(len);
  	for(size_t i=0; i<nRecs; i++)
  		for(size_t j=0; j<arrayLen; j++)
  			Tdata[i*arrayLen + j] = inObj->Tdata[j];
  }
  else if(inObj->is_event()){ 
	Edata.resize(len);
  	for(size_t i=0; i<nRecs; i++)
  		for(size_t j=0; j<arrayLen; j++)
  			Edata[i*arrayLen + j] = inObj->Edata[j];
  }
  secRes = inObj->secRes;

  this->addErr(*inObj);

}



// constructors for simple sequences
DvObject::DvObject(valarray <double> &from) 
{
  object_id=++ID_counter; // unique object id
  ref_count=0; // used by Dvar to know when safe to delete
 
    // default to sequence of single elements
  	seqLen = from.size();
  
  Ddata.resize(from.size());
  Ddata = from;
  secRes = 3;
  dims.clear();
  dims.push_back(1);
}
DvObject::DvObject(valarray <int> &from) 
{
  object_id=++ID_counter; // unique object id
  ref_count=0; // used by Dvar to know when safe to delete
 
    // default to sequence of single elements
  	seqLen = from.size();
  
  Idata.resize(from.size());
  Idata = from;
  secRes = 3;
  dims.clear();
  dims.push_back(1);
}
DvObject::DvObject(valarray <DvString> &from) 
{
  object_id=++ID_counter; // unique object id
  ref_count=0; // used by Dvar to know when safe to delete
 
    // default to sequence of single elements
  	seqLen = from.size();
  
  Sdata.resize(from.size());
  Sdata = from;
  secRes = 3;
  dims.clear();
  dims.push_back(1);
}
DvObject::DvObject(valarray <DvTime> &from) 
{
  object_id=++ID_counter; // unique object id
  ref_count=0; // used by Dvar to know when safe to delete

    // default to sequence of single elements
  	seqLen = from.size();
  
  Tdata.resize(from.size());
  Tdata = from;
  if(from.size() > 0) secRes = from[0].getSecResolution();
  dims.clear();
  dims.push_back(1);
}
DvObject::DvObject(valarray <DvEvent> &from) 
{
  object_id=++ID_counter; // unique object id
  ref_count=0; // used by Dvar to know when safe to delete
 
  // default to sequence of single elements
  seqLen = from.size();
  
  Edata.resize(from.size());
  Edata = from;
  if(from.size() > 0) secRes = from[0].start().getSecResolution();
  dims.clear();
  dims.push_back(1);
}

// constructors from valarrays 
DvObject::DvObject(valarray <double> &from, vector <size_t> &dims_, size_t nRecs) 
{
  object_id=++ID_counter; // unique object id
  ref_count=0; // used by Dvar to know when safe to delete

  seqLen = nRecs;
  dims = dims_;
  
  Ddata.resize(from.size());
  Ddata = from;
  secRes = 3;
}

DvObject::DvObject(valarray <int> &from, vector <size_t> &dims_, size_t nRecs) 
{
  object_id=++ID_counter; // unique object id
  ref_count=0; // used by Dvar to know when safe to delete

  seqLen = nRecs;
  dims = dims_;
  
  Idata.resize(from.size());
  Idata = from;
  secRes = 3;
}
DvObject::DvObject(valarray <DvString> &from, vector <size_t> &dims_, size_t nRecs) 
{
  object_id=++ID_counter; // unique object id
  ref_count=0; // used by Dvar to know when safe to delete

  seqLen = nRecs;
  dims = dims_;
  
  Sdata.resize(from.size());
  Sdata = from;
  secRes = 3;
}
DvObject::DvObject(valarray <DvTime> &from, vector <size_t> &dims_, size_t nRecs) 
{
  object_id=++ID_counter; // unique object id
  ref_count=0; // used by Dvar to know when safe to delete
 
  seqLen = nRecs;
  dims = dims_;
  
  Tdata.resize(from.size());
  Tdata = from;
  if(from.size() > 0) secRes = from[0].getSecResolution();
}
DvObject::DvObject(valarray <DvEvent> &from, vector <size_t> &dims_, size_t nRecs) 
{
  object_id=++ID_counter; // unique object id
  ref_count=0; // used by Dvar to know when safe to delete
 
  seqLen = nRecs;
  dims = dims_;
  
  Edata.resize(from.size());
  Edata = from;
  if(from.size() > 0) secRes = from[0].start().getSecResolution();
}

//  regular sequence in steps of h
void DvObject::assign_regular(double h){
    // starts from first value and increments by h
    if(is_dbl()){
        for(size_t i=1; i<Ddata.size(); i++) Ddata[i] = Ddata[i-1] + h;
    }
    else if(is_time()){
        for(size_t i=1; i<Tdata.size(); i++) Tdata[i] = Tdata[i-1] + h;
    }
    else{
        this->error("[assign regular] Not valid for data type");
    }
}

//  regular sequence in steps of h
void DvObject::assign_regular(int h){
    // starts from first value and increments by h
    if(is_int()){
        for(size_t i=1; i<Idata.size(); i++) Idata[i] = Idata[i-1] + h;
    }
    else{
        this->error("[assign regular (int)] Not valid for data type");
    }
}


// record copy constructor

DvObject_var DvObject::subSequence(size_t nRec){

  // makes a copy of data at record nRec
  DvObject_var res;
		
  if(nRec >= seqLen) {
	return new DvObject(); // return empty object
  }
  
  size_t Len = 1;
  
  int arrLen = this->arraySize();
  slice sl(nRec*arrLen, arrLen, 1);
	
  if (is_dbl()) {
    valarray <double> data = Ddata[sl];
  	res = new DvObject(data, this->Dims(), Len);
  }
  else if(is_int()) {
    valarray <int> data = Idata[sl];
  	res = new DvObject(data, this->Dims(), Len);
  }
  else if(is_str()) {
  	res = new DvObject(DvString(), this->Dims(), Len);
	res->Sdata = Sdata[sl];
  }
  else if(is_time()) {
    valarray <DvTime> data = Tdata[sl];
  	res = new DvObject(data, this->Dims(), Len);
  }
  else if(is_event()) {
    valarray <DvEvent> data = Edata[sl];
  	res = new DvObject(data, this->Dims(), Len);
  }
  else res = new DvObject();
  
  // xrefs
  DvNode *xref = this->first_xref();
  while(xref){	
  	// take xref for these records	
	if(seqSize() == xref->obj()->seqSize()){
		DvObject_var newXref = xref->obj()->subSequence(nRec);
		res->change_xref(xref->name(), newXref);
	}
	else res->change_xref(xref->name(), xref->obj());
	
	xref = xref->next;	
  }
	
  	res->secRes = secRes;
	
	return res;
}

// -----------------

// record range copy constructor
DvObject_var DvObject::subSequence(DvEvent &range){
 
      // Subset  ts object on interval
	  
	  int start, end;
	  
	  this->getEventRecBound(start, end, range);	  	  
	  
	  return this->subSequence(start, end);

}

// -----------------

// record interval copy constructor

DvObject_var DvObject::subSequence(size_t fromRec, size_t toRec){
  // makes a copy of data in record range inclusively
	
  if(fromRec == 0 && toRec == seqLen-1) return this;
  DvObject_var res;

  if(toRec < fromRec || toRec >= seqLen) {
	return new DvObject(); // return empty object
  }
  
  size_t Len = toRec - fromRec + 1;
  
  int arrLen = arraySize();
  slice sl(fromRec*arrLen, arrLen*Len, 1);

  res = this->subSequence(sl, Len);

  // xrefs
  DvNode *xref = this->first_xref();
  while(xref){	
	
  	// take xref for these records	
	if(seqSize() == xref->obj()->seqSize()){
		DvObject_var newXref = xref->obj()->subSequence(fromRec, toRec);
		res->change_xref(xref->name(), newXref);
	}
	else res->change_xref(xref->name(), xref->obj());
	
	xref = xref->next;	
  }

  return res;
}

// -----------------

// slice copy constructor
DvObject_var DvObject::subSequence(slice &sl, size_t len){

	DvObject_var res;
	if(len == 0) len = sl.size() / arraySize();
	
	if (is_dbl()) {
		valarray <double> data(Ddata[sl]);
		res = new DvObject(data, this->Dims(), len);
	}
	else if(is_int()) {
		valarray <int> data(Idata[sl]);
		res = new DvObject(data, this->Dims(), len);
	}
	else if(is_str()) {
		res = new DvObject(DvString(), this->Dims(), len);
		res->Sdata = Sdata[sl];
	}
	else if(is_time()) {
		valarray <DvTime> data(Tdata[sl]);
		res = new DvObject(data, this->Dims(), len);
	}
	else if(is_event()) {
		valarray <DvEvent> data(Edata[sl]);
		res = new DvObject(data, this->Dims(), len);
	}
	else res = new DvObject();
  
    res->secRes = secRes;

	return res;
}
// -------------------------------------------------------------------

// record range bounds
bool DvObject::getEventRecBound(int &recStart, int &recEnd, DvEvent &range){
 
	// Find records bounding specified time interval using binary search on two ends
	  	  	
	DvTime tStart = range.start_abs();
	DvTime tEnd = range.end_abs();
	DvObject_var tt = this->getTimeTags();
	if(tt.is_nil()){
		recStart = -1;
		recEnd = -1;
		return false;
	}
	
	int recL=0;
	int recU=tt->seqSize()-1;
	
	if(tStart > tt->time(recU) || tEnd < tt->time(recL)) {
		recStart = -1;
		recEnd = -1;		
		return false; // ivl outside tt range
	}
	
	if(tStart < tt->time(recL) ) recStart = 0; // outside tt range
	else{
		tt->getBounds(recL, recU, tStart);
		if(tt->time(recL) == tStart ) recStart = recL;
		else recStart = recU;
	}

	recL=0;
	recU=tt->seqSize()-1;
	if(tEnd > tt->time(recU) ) recEnd = recU; // outside tt range
	else{
		tt->getBounds(recL, recU, tEnd);
		if(tt->time(recU) == tEnd ) recEnd = recU;
		else recEnd = recL;
	}
	return true;
}

// -------------------------------------------------------------------

// record range bounds
bool DvObject::getRecBounds(int &recStart, int &recEnd, double start, double end){
 
	// Find records bounding specified scalar interval using binary search on two ends
	if( start > end) {
		double hold = start;
		start = end;
		end = hold;
	}
	
	DvObject_var ss = this->getDep0();
	if(ss.is_nil() || !ss->is_dbl() ){
		recStart = -1;
		recEnd = -1;
		return false;
	}
	
	int recL=0;
	int recU=ss->seqSize()-1;
	
	if(start > ss->dbl(recU) || end < ss->dbl(recL)) {
		recStart = -1;
		recEnd = -1;		
		return false; //  outside range
	}
	
	if(start < ss->dbl(recL) ) recStart = 0; // outside range
	else{
		ss->getBounds(recL, recU, start);
		if(ss->dbl(recL) == start ) recStart = recL;
		else recStart = recU;
	}

	recL=0;
	recU=ss->seqSize()-1;
	if(end > ss->dbl(recU) ) recEnd = recU; // outside tt range
	else{
		ss->getBounds(recL, recU, end);
		if(ss->dbl(recU) == end ) recEnd = recU;
		else recEnd = recL;
	}
	
	return true;
}

// -------------------------------------------------------------------

void DvObject::getBounds(int &recL, int &recU, const DvTime& t){

	// Find the records that bound a time value using recursive binary search

	if( recU == recL+1 || recL == recU) return; // found bounding records
	
	int recH =  (recU + recL) / 2;
	if( this->asTime(recH) < t) recL = recH;
	else recU = recH;
	
	return this->getBounds(recL, recU, t);
}

// -------------------------------------------------------------------

void DvObject::getBounds(int &recL, int &recU, double d){

	// Find the records that bound a value using recursive binary search

	if( recU == recL+1 || recL == recU) return; // found bounding records
	
	int recH =  (recU + recL) / 2;
	if( this->asDouble(recH) < d) recL = recH;
	else recU = recH;
	
	return this->getBounds(recL, recU, d);
}



DvObject * DvObject::create(const vector <size_t> &dims, size_t nRecs){
    
    if(is_dbl()) {
        return new DvObject((double)0.0, dims, nRecs);
    }
    else if(is_int()) {
        return new DvObject((int)0, dims, nRecs);
    }
    else if(is_str()) {
        DvString nullStr = DvString();
        return new DvObject(nullStr, dims, nRecs);
    }
    else if(is_time()) {
        DvTime nullTime = DvTime();
        return new DvObject(nullTime, dims, nRecs);
    }
    else if(is_event()) {
        DvEvent nullEvent = DvEvent();
        return new DvObject(nullEvent, dims, nRecs);
    }
    
    return new DvObject();
    
    
}

DvObject_var DvObject::resize(size_t nRecs) 
{
    size_t copySize = this->seqSize() * this->arraySize();
    if(this->seqSize() > nRecs) copySize = nRecs * this->arraySize();
    
    DvObject_var ret = this->create(nRecs);
    
    ret->copy_xrefs_from(*this);
    
    // copy data as far as we have
    if(this->is_dbl()){
        for(size_t i=0; i<copySize; i++) ret->dbl(i) = this->Ddata[i];
    }
    else if(this->is_int()){  
        for(size_t i=0; i<copySize; i++) ret->itg(i) = this->Idata[i];
    }
    else if(this->is_str()){  
        for(size_t i=0; i<copySize; i++) ret->str(i) = this->Sdata[i];
    }
    else if(this->is_time()){ 
        for(size_t i=0; i<copySize; i++) ret->time(i) = this->Tdata[i];
        // for time and event repeat last value to avoid non-monotonic result
        for(size_t i=copySize; i<ret->totalSize(); i++) ret->time(i) = this->Tdata[copySize-1];        
    }
    else if(this->is_event()){ 
        for(size_t i=0; i<copySize; i++) ret->event(i) = this->Edata[i];
        for(size_t i=copySize; i<ret->totalSize(); i++) ret->event(i) = this->Edata[copySize-1];
    }
    
    ret->secRes = this->secRes;
    
    ret->addErr(*this);
    
    DvNode *xref = ret->xrefs.first();
    while(xref){
        if(xref->obj()->seqSize() == this->seqSize()) xref->obj() = xref->obj()->resize(nRecs);
        xref = xref->next;
    }
    
    return ret;
    
}



// fast element access (allows modification of element)
// use after dimensionality and type known
//---------------------------------------
// THESE ARE THE ONLY UNPROTECTED METHODS


double & DvObject::dbl(size_t nRec, size_t i, size_t j, size_t k, size_t m){
	return Ddata[m + dims[3]*(k + dims[2]*(j + dims[1]*(i + dims[0]*nRec)))];
}

double & DvObject::dbl(size_t nRec, size_t i, size_t j, size_t k){
	return Ddata[k + dims[2]*(j + dims[1]*(i + dims[0]*nRec))];
}

double & DvObject::dbl(size_t nRec, size_t i, size_t j){
	return Ddata[j + dims[1]*(i + dims[0]*nRec)];
}

double & DvObject::dbl(size_t nRec, size_t i){
	// may be 1D or counter over all array elements
	return Ddata[i + arraySize()*nRec];
}

double & DvObject::dbl(size_t posn){
	return Ddata[posn];
}

double & DvObject::dbl(size_t nRec, vector <size_t> &index){

	if(index.size() != dims.size() || index.size() == 0) return Ddata[nRec];
	
	size_t posn = nRec*arraySize();
	size_t element = 0;
	for(size_t i=0; i<index.size(); i++){
		element += index[i];
		if(i+1 < index.size()) element *= dims[i+1];
	}
	posn += element;
	return Ddata[posn];
}

// integer type
int & DvObject::itg(size_t nRec, size_t i, size_t j, size_t k, size_t m){
	return Idata[m + dims[3]*(k + dims[2]*(j + dims[1]*(i + dims[0]*nRec)))];
}

int & DvObject::itg(size_t nRec, size_t i, size_t j, size_t k){
	return Idata[k + dims[2]*(j + dims[1]*(i + dims[0]*nRec))];
}

int & DvObject::itg(size_t nRec, size_t i, size_t j){
	return Idata[j + dims[1]*(i + dims[0]*nRec)];
}

int & DvObject::itg(size_t nRec, size_t i){
	// may be 1D or counter over all array elements
	return Idata[i + arraySize()*nRec];
}

int & DvObject::itg(size_t posn){
	return Idata[posn];
}

int & DvObject::itg(size_t nRec, vector <size_t> &index){

	if(index.size() != dims.size() || index.size() == 0) return Idata[nRec];
	
	size_t posn = nRec*arraySize();
	size_t element = 0;
	for(size_t i=0; i<index.size(); i++){
		element += index[i];
		if(i+1 < index.size()) element *= dims[i+1];
	}
	posn += element;
	return Idata[posn];
}

// string type
DvString & DvObject::str(size_t nRec, size_t i, size_t j, size_t k, size_t m){
	return Sdata[m + dims[3]*(k + dims[2]*(j + dims[1]*(i + dims[0]*nRec)))];
}

DvString & DvObject::str(size_t nRec, size_t i, size_t j, size_t k){
	return Sdata[k + dims[2]*(j + dims[1]*(i + dims[0]*nRec))];
}

DvString & DvObject::str(size_t nRec, size_t i, size_t j){
	return Sdata[j + dims[1]*(i + dims[0]*nRec)];
}

DvString & DvObject::str(size_t nRec, size_t i){
	// may be 1D or counter over all array elements
	return Sdata[i + arraySize()*nRec];
}

DvString & DvObject::str(size_t posn){
	return Sdata[posn];
}

DvString & DvObject::str(size_t nRec, vector <size_t> &index){

	if(index.size() != dims.size() || index.size() == 0) return Sdata[nRec];
	
	size_t posn = nRec*arraySize();
	size_t element = 0;
	for(size_t i=0; i<index.size(); i++){
		element += index[i];
		if(i+1 < index.size()) element *= dims[i+1];
	}
	posn += element;
	return Sdata[posn];
}

// time type
DvTime & DvObject::time(size_t nRec, size_t i, size_t j, size_t k, size_t m){
	return Tdata[m + dims[3]*(k + dims[2]*(j + dims[1]*(i + dims[0]*nRec)))];
}

DvTime & DvObject::time(size_t nRec, size_t i, size_t j, size_t k){
	return Tdata[k + dims[2]*(j + dims[1]*(i + dims[0]*nRec))];
}

DvTime & DvObject::time(size_t nRec, size_t i, size_t j){
	return Tdata[j + dims[1]*(i + dims[0]*nRec)];
}

DvTime & DvObject::time(size_t nRec, size_t i){
	// may be 1D or counter over all array elements
	return Tdata[i + arraySize()*nRec];
}

DvTime & DvObject::time(size_t posn){
	return Tdata[posn];
}

DvTime & DvObject::time(size_t nRec, vector <size_t> &index){

	if(index.size() != dims.size() || index.size() == 0) return Tdata[nRec];
	
	size_t posn = nRec*arraySize();
	size_t element = 0;
	for(size_t i=0; i<index.size(); i++){
		element += index[i];
		if(i+1 < index.size()) element *= dims[i+1];
	}
	posn += element;
	return Tdata[posn];
}

// event type (time interval)
DvEvent & DvObject::event(size_t nRec, size_t i, size_t j, size_t k, size_t m){
	return Edata[m + dims[3]*(k + dims[2]*(j + dims[1]*(i + dims[0]*nRec)))];
}

DvEvent & DvObject::event(size_t nRec, size_t i, size_t j, size_t k){
	return Edata[k + dims[2]*(j + dims[1]*(i + dims[0]*nRec))];
}

DvEvent & DvObject::event(size_t nRec, size_t i, size_t j){
	return Edata[j + dims[1]*(i + dims[0]*nRec)];
}

DvEvent & DvObject::event(size_t nRec, size_t i){
	// may be 1D or counter over all array elements
	return Edata[i + arraySize()*nRec];
}

DvEvent & DvObject::event(size_t posn){
	return Edata[posn];
}

DvEvent & DvObject::event(size_t nRec, vector <size_t> &index){

	if(index.size() != dims.size() || index.size() == 0) return Edata[nRec];
	
	size_t posn = nRec*arraySize();
	size_t element = 0;
	for(size_t i=0; i<index.size(); i++){
		element += index[i];
		if(i+1 < index.size()) element *= dims[i+1];
	}
	posn += element;
	return Edata[posn];
}


// SAFE versions of Fast Access routines
// fast element access (allows modification of element)
// use after dimensionality and type known
// Safe return value must be set.
//---------------------------------------


DvRecord DvObject::operator[](size_t nRec){     
    return DvRecord(this, nRec);
}

valarray <double>  DvObject::record(size_t nRec){

  // returns double data at record nRec
  // The calling module should first test that object is double.
  
  int arrLen = this->arraySize();

  if(Ddata.size() < (nRec+1)*arrLen) {
	error("[]");
	return valarray<double>(0.0, arrLen);
  }
  
  slice sl(nRec*arrLen, arrLen, 1);
	  	
  return Ddata[sl];

}

double & DvObject::dblS(size_t nRec, size_t i, size_t j, size_t k, size_t m, double safe){
	static double _safe;
	_safe = safe;
	
	if(nDims() < 4 ) return _safe;
	size_t posn = m + dims[3]*(k + dims[2]*(j + dims[1]*(i + dims[0]*nRec)));
	
	if(Ddata.size() <= posn) return _safe;
	
	return Ddata[posn];
}

double & DvObject::dblS(size_t nRec, size_t i, size_t j, size_t k, double safe){
	static double _safe;
	_safe = safe;
	
	if(nDims() < 3 ) return _safe;
	size_t posn = k + dims[2]*(j + dims[1]*(i + dims[0]*nRec));
	
	if(Ddata.size() <= posn) return _safe;
	
	return Ddata[posn];
}

double & DvObject::dblS(size_t nRec, size_t i, size_t j, double safe){
	static double _safe;
	_safe = safe;
	
	if(nDims() < 2 ) return _safe;
	size_t posn = j + dims[1]*(i + dims[0]*nRec);
	
	if(Ddata.size() <= posn) return _safe;
	
	return Ddata[posn];
}

double & DvObject::dblS(size_t nRec, size_t i, double safe){
	static double _safe;
	_safe = safe;
	
	if(nDims() < 1 ) return _safe;
	size_t posn = i + dims[0]*nRec;
	
	if(Ddata.size() <= posn) return _safe;
	
	return Ddata[posn];
}

double & DvObject::dblS(size_t posn, double safe){
	static double _safe;
	_safe = safe;
	
	if(Ddata.size() <= posn) return _safe;
 
	return Ddata[posn];
}

// integer type
int & DvObject::itgS(size_t nRec, size_t i, size_t j, size_t k, size_t m, int safe){
	static int _safe;
	_safe = safe;
	
	if(nDims() < 4 ) return _safe;
	size_t posn = m + dims[3]*(k + dims[2]*(j + dims[1]*(i + dims[0]*nRec)));
	
	if(Idata.size() <= posn) return _safe;
	
	return Idata[posn];
}

int & DvObject::itgS(size_t nRec, size_t i, size_t j, size_t k, int safe){
	static int _safe;
	_safe = safe;
	
	if(nDims() < 3 ) return _safe;
	size_t posn = k + dims[2]*(j + dims[1]*(i + dims[0]*nRec));
	
	if(Idata.size() <= posn) return _safe;
	
	return Idata[posn];
}

int & DvObject::itgS(size_t nRec, size_t i, size_t j, int safe){
	static int _safe;
	_safe = safe;
	
	if(nDims() < 2 ) return _safe;
	size_t posn = j + dims[1]*(i + dims[0]*nRec);
	
	if(Idata.size() <= posn) return _safe;
	
	return Idata[posn];
}

int & DvObject::itgS(size_t nRec, size_t i, int safe){
	static int _safe;
	_safe = safe;
	
	if(nDims() < 1 ) return _safe;
	size_t posn = i + dims[0]*nRec;
	
	if(Idata.size() <= posn) return _safe;
	
	return Idata[posn];
}

int & DvObject::itgS(size_t posn, int safe){
	static int _safe;
	_safe = safe;
		
	if(Idata.size() <= posn) return _safe;
	
	return Idata[posn];
}

// string type
DvString & DvObject::strS(size_t nRec, size_t i, size_t j, size_t k, size_t m, DvString safe){
	static DvString _safe;
	_safe = safe;
	
	if(nDims() < 4 ) return _safe;
	size_t posn = m + dims[3]*(k + dims[2]*(j + dims[1]*(i + dims[0]*nRec)));
	
	if(Sdata.size() <= posn) return _safe;
	
	return Sdata[posn];
}

DvString & DvObject::strS(size_t nRec, size_t i, size_t j, size_t k, DvString safe){
	static DvString _safe;
	_safe = safe;
	
	if(nDims() < 3 ) return _safe;
	size_t posn = k + dims[2]*(j + dims[1]*(i + dims[0]*nRec));
	
	if(Sdata.size() <= posn) return _safe;
	
	return Sdata[posn];
}

DvString & DvObject::strS(size_t nRec, size_t i, size_t j, DvString safe){
	static DvString _safe;
	_safe = safe;
	
	if(nDims() < 2 ) return _safe;
	size_t posn = j + dims[1]*(i + dims[0]*nRec);
	
	if(Sdata.size() <= posn) return _safe;
	
	return Sdata[posn];
}

DvString & DvObject::strS(size_t nRec, size_t i, DvString safe){
	static DvString _safe;
	_safe = safe;
	
	if(nDims() < 1 ) return _safe;
	size_t posn = i + dims[0]*nRec;
	
	if(Sdata.size() <= posn) return _safe;
	
	return Sdata[posn];
}

DvString & DvObject::strS(size_t posn, DvString safe){
	static DvString _safe;
	_safe = safe;
		
	if(Sdata.size() <= posn) return _safe;
	
	return Sdata[posn];
}

// time type
DvTime & DvObject::timeS(size_t nRec, size_t i, size_t j, size_t k, size_t m, DvTime safe){
	static DvTime _safe;
	_safe = safe;
	
	if(nDims() < 4 ) return _safe;
	size_t posn = m + dims[3]*(k + dims[2]*(j + dims[1]*(i + dims[0]*nRec)));
	
	if(Tdata.size() <= posn) return _safe;
	
	return Tdata[posn];
}

DvTime & DvObject::timeS(size_t nRec, size_t i, size_t j, size_t k, DvTime safe){
	static DvTime _safe;
	_safe = safe;
	
	if(nDims() < 3 ) return _safe;
	size_t posn = k + dims[2]*(j + dims[1]*(i + dims[0]*nRec));
	
	if(Tdata.size() <= posn) return _safe;
	
	return Tdata[posn];
}

DvTime & DvObject::timeS(size_t nRec, size_t i, size_t j, DvTime safe){
	static DvTime _safe;
	_safe = safe;
	
	if(nDims() < 2 ) return _safe;
	size_t posn = j + dims[1]*(i + dims[0]*nRec);
	
	if(Tdata.size() <= posn) return _safe;
	
	return Tdata[posn];
}

DvTime & DvObject::timeS(size_t nRec, size_t i, DvTime safe){
	static DvTime _safe;
	_safe = safe;
	
	if(nDims() < 1 ) return _safe;
	size_t posn = i + dims[0]*nRec;
	
	if(Tdata.size() <= posn) return _safe;
	
	return Tdata[posn];
}

DvTime & DvObject::timeS(size_t posn, DvTime safe){
	static DvTime _safe;
	_safe = safe;
		
	if(Tdata.size() <= posn) return _safe;
	
	return Tdata[posn];
}

// event type (time interval)
DvEvent & DvObject::eventS(size_t nRec, size_t i, size_t j, size_t k, size_t m, DvEvent safe){
	static DvEvent _safe;
	_safe = safe;
	
	if(nDims() < 4 ) return _safe;
	size_t posn = m + dims[3]*(k + dims[2]*(j + dims[1]*(i + dims[0]*nRec)));
	
	if(Edata.size() <= posn) return _safe;
	
	return Edata[posn];
}

DvEvent & DvObject::eventS(size_t nRec, size_t i, size_t j, size_t k, DvEvent safe){
	static DvEvent _safe;
	_safe = safe;
	
	if(nDims() < 3 ) return _safe;
	size_t posn = k + dims[2]*(j + dims[1]*(i + dims[0]*nRec));
	
	if(Edata.size() <= posn) return _safe;
	
	return Edata[posn];
}

DvEvent & DvObject::eventS(size_t nRec, size_t i, size_t j, DvEvent safe){
	static DvEvent _safe;
	_safe = safe;
	
	if(nDims() < 2 ) return _safe;
	size_t posn = j + dims[1]*(i + dims[0]*nRec);
	
	if(Edata.size() <= posn) return _safe;
	
	return Edata[posn];
}

DvEvent & DvObject::eventS(size_t nRec, size_t i, DvEvent safe){
	static DvEvent _safe;
	_safe = safe;
	
	if(nDims() < 1 ) return _safe;
	size_t posn = i + dims[0]*nRec;
	
	if(Edata.size() <= posn) return _safe;
	
	return Edata[posn];
}

DvEvent & DvObject::eventS(size_t posn, DvEvent safe){
	static DvEvent _safe;
	_safe = safe;
		
	if(Edata.size() <= posn) return _safe;
	return Edata[posn];
}

void DvObject::setRecord(size_t nRec, size_t oRec, DvObject* obj){
	// No defensive programming, types, dimensions must match and records in range 
	size_t nArr = arraySize();
	if( is_dbl() ){
		for(size_t i=0; i<nArr; i++) Ddata[nArr*nRec+i] = obj->Ddata[nArr*oRec+i];
	}
	else if( is_time() ){
		for(size_t i=0; i<nArr; i++) Tdata[nArr*nRec+i] = obj->Tdata[nArr*oRec+i];
	}
	else if( is_int() ){
		for(size_t i=0; i<nArr; i++) Idata[nArr*nRec+i] = obj->Idata[nArr*oRec+i];
	}
	else if( is_event() ){
		for(size_t i=0; i<nArr; i++) Edata[nArr*nRec+i] = obj->Edata[nArr*oRec+i];
	}
	else if( is_str() ){
		for(size_t i=0; i<nArr; i++) Sdata[nArr*nRec+i] = obj->Sdata[nArr*oRec+i];
	}
	  
}
