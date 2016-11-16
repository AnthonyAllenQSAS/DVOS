//
//  DvObj.h
//
//  Dvos -  Data variable object classes
//
//  Utility calls 
//
//  Tony Allen
//  Nov 2012
//

#include <iostream>

#include "DvObj.h"
#include "DvRecord.h"

using namespace std;

const double DV_CLOSE_ENOUGH = 2.e-7; 
const double DvNaN = std::numeric_limits<double>::quiet_NaN();
const vector <size_t> emptyDim;
const vector <size_t> DIM3(1, 3); // (3) vector
const vector <size_t> DIM33(2,3); // (3,3)
const vector <size_t> DIM4(1, 4); // (4)
const vector <size_t> DIM44(2,4); // (4,4)

double twopi = 8.*std::atan(1.0);

// ----------------  General translation ---------------- 

DvString DvObject::getXrefText(DvString &name){

	if( !xref_exists(name) ) return DvString("");
	
	DvObject_var xref_ptr = this->get_xref(name);
	DvString ret;
	if( xref_ptr.is_ok()) ret = xref_ptr->asStr(0);
	for(size_t i=1; i<xref_ptr->totalSize(); i++) {
		ret += ", ";
		ret += xref_ptr->asStr(i);
	}
	return ret;
}

DvString DvObject::getXrefText(const char *name){

	if( !xref_exists(name) ) return DvString("");

	DvObject_var xref_ptr = this->get_xref(name);
	
	DvString ret;
	// convert each element to string
	if( xref_ptr.is_ok()) ret = xref_ptr->asStr(0);
	for(size_t i=1; i<xref_ptr->totalSize(); i++) {
		ret += ", ";
		ret += xref_ptr->asStr(i);
	}
	
	return ret;
}

int DvObject::get_iFILL(){

	if( !xref_exists(FILL_VALUE) ) return -1;

	DvObject_var xref_ptr = this->get_xref(FILL_VALUE);
	
	int iFILL = xref_ptr->asInt(0);
	
	return iFILL;
}

double DvObject::get_dFILL(){

    if( !xref_exists(FILL_VALUE) ) return -1.;

	DvObject_var xref_ptr = this->get_xref(FILL_VALUE);
	
	double dFILL = xref_ptr->asDouble(0);
	
	return dFILL;
}



// ---------------- DEPEND_0 handling (Time Tags or Scalar Sequence) -----------------

DvObject_var DvObject::getDep0(){
	return get_xref(DEPEND_0);
}


DvObject_var DvObject::getTimeTags(){

	if(this->is_time()) return DvObject_var(this);

	DvObject_var ret = getDep0();
	if(ret->is_time()) return ret;
	
	return DvObject_var();
}

bool DvObject::hasTimeTags(){
	
	DvObject_var D0 = get_xref(DEPEND_0);
	if(D0->is_time()) return true;

	return false;
}


// ---------------- 
DvEvent DvObject::getTimeRange(){

	if(this->is_event()) {
		DvEvent ret(this->event(0).start(), this->event(seqLen-1).end());
		return ret;
	}
	else if(this->is_time()) {
		DvEvent ret(this->time(0), this->time(seqLen-1));
		return ret;
	}
	
	DvObject_var tt = this->getTimeTags();

	if(this->is_time()){
		tt = this;
	}
	
	if(tt->is_nil()) return DvEvent();
	
	DvEvent ret( tt->time(0), tt->time(tt->seqSize()-1) );

	return ret;
}

// ---------------- 

bool DvObject::getDataRange(size_t recFrom, size_t recTo, double &min, double &max){
	
	if(this->is_nil()) return false;
	
	if(recTo > seqSize() ) return false;

	size_t len = arraySize();
	
	if(this->is_dbl()) {
		// find first valid data
		while(recFrom < recTo && isnan(Ddata[recFrom]) ) recFrom++;
		if( isnan(Ddata[recFrom]) ) return false;
		
		min = Ddata[recFrom*len];
		max = min;
		for(size_t i=(size_t)recFrom; i<=(size_t)recTo; i++){
			for(size_t j=0; j<len; j++){
				double d = Ddata[i*len + j];
				if( !isnan(d) && d < min) min = d;
				if( !isnan(d) && d > max) max = d;
			}
		}
		return true;
	}
	if(this->is_int()) {
		min = (double) Idata[recFrom*len];
		max = min;
		for(size_t i=(size_t)recFrom; i<=(size_t)recTo; i++){
			for(size_t j=0; j<len; j++){
				double d = (double) Idata[i*len + j];
				if(d < min) min = d;
				if(d > max) max = d;
			}
		}
		return true;
	}
	
	if(this->is_time()){
		min = (double) Tdata[recFrom*len].getEpoch2000();
		max = min;
		for(size_t i=(size_t)recFrom; i<=(size_t)recTo; i++){
			for(size_t j=0; j<len; j++){
				double d = (double) Tdata[i*len + j].getEpoch2000();
				if(d < min) min = d;
				if(d > max) max = d;
			}
		}
		return true;
	}
	
	return false;

}

bool DvObject::getDataRange(double &min, double &max){
	
	if(this->is_nil()) return false;
	
	if(this->is_dbl() || this->is_int()) {
		min = this->min();
		max = this->max();
		return true;
	}
	
	if(this->is_time()){
		DvEvent ivl = this->getTimeRange();
		min = ivl.start().getEpoch2000();
		max = ivl.end().getEpoch2000();
		return true;
	}
	
	return false;

}


// ----------------  Coordinate system ---------------- 

bool DvObject::sameFrame(DvObject &arg){

	DvString thisFrame = this->getFrameAttr();
	if(thisFrame.empty()) return true; // no idea, so proceed

	DvString argFrame = arg.getFrameAttr();
	if(argFrame.empty()) return true; // no idea, so proceed
	
    if(thisFrame.contains("scalar") || argFrame.contains("scalar") ) return true; // scalar, so safe

    if(thisFrame.contains("any") || argFrame.contains("any") ) return true; // any frame, generic vector safe

	if(thisFrame.after('>') == argFrame.after('>')) return true; // operator== is case insensitive
	return false; // frames exist and differ
}


bool DvObject::sameFrame(DvString &argFrame){

	DvString thisFrame = this->getFrameAttr();
	if(thisFrame.empty()) return true; // no idea, so proceed

	if(argFrame.empty()) return true; // no idea, so proceed
	
	if(thisFrame.contains("scalar") || argFrame.contains("scalar") ) return true; // scalar, so safe
 
    if(thisFrame.contains("any") || argFrame.contains("any") ) return true; // any frame, generic vector safe

	if(thisFrame.after('>') == argFrame.after('>')) return true; // operator== is case insensitive
		
	return false; // frames exist and differ
}

DvString DvObject::getFrameAttr(){
	// get the CSDS Frame syntax, e.g. vector>gse_xyz

	// if CSDS Frame attribute exists get it
	DvString frm = this->getXrefText(FRAME);
	if( !frm.empty() ) return frm;
	
	// else construct from CAA syntax
	DvString frameAttr("");
	
	frm = getFrame();
	if( frm.empty() ) return frameAttr;

	DvString rep   = getRep();
	if( rep.empty() ) return frameAttr;	 
	
	int order = getOrder();

	if(order == 1 ){
		frameAttr = "vector>";
		frameAttr += frm;
		frameAttr += "_";
		frameAttr += rep;
	}
	else if(order == 2 ){
		frameAttr = "tensor>";
		frameAttr += frm;
		frameAttr += "_";
		frameAttr += rep;
	}
	else if(order == 0 ){
		frameAttr = "scalar>na";
	}
	return frameAttr;
}

void DvObject::setFrameAttr(DvString frameAttr){
	// set the CSDS and CAA Frame from CSDS syntax, e.g. vector>gse_xyz

	change_xref(FRAME, frameAttr );
	
	// FRAME attribute now exists, so can get components as normal and set CAA attributes

	DvString frm   = getFrame();
	DvString rep   = getRep();
	int order = getOrder();

	change_xref(REPRESENTATION, rep);
	if(rep.length() > 1) {
		vector <size_t> dim;
		dim.push_back(rep.length());
		DvObject_var rep_ptr = new DvObject(rep, dim);
		for(size_t i=0; i<rep.length(); i++){
			rep_ptr->str(i) = rep(i);
		}
		change_xref(REPRESENTATION_1, rep_ptr);
		if(order == 2) change_xref(REPRESENTATION_2, rep_ptr);

	}
	
	change_xref(TENSOR_ORDER, order );

	if(frm == "GSE" ) frm = "GSE>Geocentric Solar Ecliptic";
	else if(frm == "GSM" ) frm = "GSM>Geocentric Solar Magnetic";
	else if(frm == "SR2" ) frm = "SR2>Spacecraft Despun";
	else if(frm == "ISR2" ) frm = "ISR2>Inverted Spacecraft Despun";
	
	change_xref(COORDINATE_SYSTEM, frm );
	
}

void DvObject::setFrameAttr(DvString frm, DvString rep, int order){
	
	// set the CSDS and CAA Frame from CAA syntax

	change_xref(REPRESENTATION, rep);
	if(rep.length() > 1) {
		vector <size_t> dim;
		dim.push_back(rep.length());
		DvObject_var rep_ptr = new DvObject(rep, dim);
		for(size_t i=0; i<rep.length(); i++){
			rep_ptr->str(i) = rep(i);
		}
		change_xref(REPRESENTATION_1, rep_ptr);

	}
	
	change_xref(TENSOR_ORDER, order );

	if(frm == "GSE" ) frm = "GSE>Geocentric Solar Ecliptic";
	else if(frm == "GSM" ) frm = "GSM>Geocentric Solar Magnetic";
	else if(frm == "SR2" ) frm = "SR2>Spacecraft Despun";
	else if(frm == "ISR2" ) frm = "ISR2>Inverted Spacecraft Despun";
	
	change_xref(COORDINATE_SYSTEM, frm );
	
	// delete existing attribute if it exists (it may be out of date)
	delete_xref(FRAME);
	
	// get new version (from CAA attrs which now exist)
	DvString frameAttr = getFrameAttr();
	change_xref(FRAME, frameAttr );
	
}

void DvObject::setToFrameAttr(const char *frame){

	DvString Attribute_str = "vector>";
 	Attribute_str += frame;
  	Attribute_str += "_xyz";
   
	this->change_xref("ToFrame", Attribute_str);
   
} 



 void DvObject::setAttrsProduct(DvObject &arg){
	
	// sets frame for this * arg
	
	DvString frameAttr("");
 
	// This is used for non-vector operation multiply
	// input may be rotation matrix 
 
	DvString argFrm = arg.getFrameAttr();
	bool argIsRotn = false;
	if (argFrm.contains("ROT")) argIsRotn = true;
 
	DvString thisFrm = this->getFrameAttr();
	bool thisIsRotn = false;
	if (thisFrm.contains("ROT")) thisIsRotn = true;
 
  
	bool thisIsVec = this->isThreeVector();
	bool argIsVec = arg.isThreeVector();

	// if both 3 vectors use scalar (cross product separate operator)

	if(thisIsVec && argIsRotn){
		frameAttr = arg.getXrefText("ToFrame");
	}
	else if(thisIsVec)
    	frameAttr = thisFrm;
	else if(thisIsRotn && argIsVec)
		frameAttr = this->getXrefText("ToFrame");
  	else if(argIsVec)
		frameAttr = argFrm;
	else if(thisIsRotn && argIsRotn){
	
		// This is usual form Vec = Rot * Vec 
		frameAttr = arg.getFrameAttr();   // R2 is applied 1st to vector, so from Frame from R2
    	DvString toFrame = this->getXrefText("ToFrame");// ToFrame will be from R1 
		change_xref("ToFrame", toFrame);
	} 
  	setFrameAttr(frameAttr);


 	// DEPEND_i
 	// special case of one operand is a single scalar
 	if ( arg.seqSize() == 1 )
 	{
 		// do nothing dims are already correct
 	}
 	else if ( arg.arraySize() == 1 )
 	{
        // 2nd obj is scalar use depends from first
 		for( size_t i=1; i<= arg.nDims(); i++ ){
 			DvString dep_str = "DEPEND_";
 			dep_str += (int) i;
 			DvObject_var depi = arg.get_xref(dep_str);
 			if(depi.is_ok()) this->change_xref(dep_str, depi);
 		}
 	}
 	else{
 		// handle matrix
 		
        // this has nDims 1st + nDims 2nd - 2
        
 		size_t iDim = this->nDims()-arg.nDims()+2; // will replace last dim with first dim of obj2
        DvString dep_str = "DEPEND_";
        DvString rep_str = "REPRESENTATION_";
        dep_str += (int) iDim;
        rep_str += (int) iDim;
        this->delete_xref(dep_str); // in case both objects have only 1 dim
        this->delete_xref(rep_str); 

        size_t i=2;
 		while( i <= arg.nDims() ){ // this ignores first dim of obj2

            DvString dep_str = "DEPEND_";
            DvString dep_strm1 = "DEPEND_";
            DvString rep_str = "REPRESENTATION_";
            DvString rep_strm1 = "REPRESENTATION_";
            dep_str += (int) i;
            dep_strm1 += (int) iDim;
            rep_str += (int) i;
            rep_strm1 += (int) iDim;
            DvObject_var depi = arg.get_xref(dep_str);
            DvObject_var repi = arg.get_xref(rep_str);
            if(depi.is_ok()) this->change_xref( dep_strm1, depi  );
            if(repi.is_ok()) this->change_xref( rep_strm1, repi  );
 			iDim++;
 			i++;
 		}
 	}
	
     if(this->nDims() != 2) delete_xref(TOFRAME);
	return;

}

DvString DvObject::getFrame(){
	// just the actual frame element, eg gse

	DvString attrPtr = getXrefText( FRAME );
  
	if( attrPtr.empty() ){
		// didn't find frame, look for COORDINATE_SYSTEM
		attrPtr = getXrefText( COORDINATE_SYSTEM );
		attrPtr = attrPtr.front(">"); // gets, eg gse from "gse>geocentric solar ecliptic"

	}
	else{
	 	    
		attrPtr = attrPtr.after('>'); // gets, eg gse_xyz, may be empty
		attrPtr = attrPtr.before('_'); // gets, eg gse, may be empty
	 
	}

	return attrPtr;
}

DvString DvObject::getRep(){
	// just the actual rep element, eg xyz
	   
	DvString attrPtr = getXrefText( FRAME );

	if( attrPtr.empty() ){
		// didn't find frame, look for REPRESENTATION
		attrPtr = getXrefText( REPRESENTATION );
		if( attrPtr.empty()){
			attrPtr = getXrefText(REPRESENTATION_1);
			attrPtr.strip(' ');
			attrPtr.strip(',');
		}
		// catch Themis syntax
		attrPtr = attrPtr.replace("Rep_","");
		attrPtr = attrPtr.replace("rep_","");
		return  attrPtr; // may be empty
	}
	else{
	 	    
		attrPtr = attrPtr.after('_'); // gets, eg xyz, may be empty
   }

   return attrPtr;
}

int DvObject::getOrder(){
	// just the actual order element, eg 1 for vector
      
	// use Frame attribute if present
	DvString attrPtr = getXrefText( FRAME );
	if( !attrPtr.empty()){
		DvString order = attrPtr.before('>');
		if(order == "VECTOR") return 1;
		else if(order == "ROT" || order == "TENSOR") return 2;
		else return 0; // default, also empty frame attribute
	}
	else {
		DvObject_var xref = get_xref(TENSOR_ORDER);
   		if( xref.is_ok() ) return xref->asInt(); // get integer value	
	}
  
   return 0;
}


// ----------------  UNITS and SI_conversion ---------------- 

bool DvObject::sameUnits(DvObject &arg){
	// same unit and conversion factor
	DvString base_str1 = this->getBaseSI();  
    DvString base_str2 = arg.getBaseSI();
	
	if( base_str1.empty() || base_str2.empty() ) return true; // accept numeric value w/out xref
	
    if(base_str1.after('>') != base_str2.after('>')) return false;
	
	double conv1 = base_str1.before('>').toDouble();
	double conv2 = base_str2.before('>').toDouble();
   
	if( fabs( (conv1 - conv2) / (conv1 + conv2) ) > DV_CLOSE_ENOUGH ) return false; // fuzzy equality test
	
	return true;

}

bool DvObject::sameBaseSI(DvObject &arg){
	// same base form, may need converting
	DvString base_str1 = this->getBaseSI().after('>');  
    DvString base_str2 = arg.getBaseSI().after('>');
	if( base_str1.empty() || base_str2.empty() ) return true; // accept numeric value w/out xref
	
    if(base_str1 != base_str2) return false;
   
	return true;

}
bool DvObject::sameUnits(const char * SIstr, int i){
	// same unit and conversion factor
	// int i allows specification of which dimension in, e.g. RPL vector
	
	DvString base_str1 = this->getBaseSI(i);  
    DvString base_str2 = DvUnitGetBaseSI(SIstr);
	
	if( base_str1.empty() || base_str2.empty() ) return true; // accept numeric value w/out xref
	
    if(base_str1.after('>') != base_str2.after('>')) return false;
	
	double conv1 = base_str1.before('>').toDouble();
	double conv2 = base_str2.before('>').toDouble();
   
	if( fabs( (conv1 - conv2) / (conv1 + conv2) ) > DV_CLOSE_ENOUGH ) return false; // fuzzy equality test
	
	return true;

}

bool DvObject::sameBaseSI(const char *  SIstr){
	// same base form, may need converting
	DvString base_str1 = this->getBaseSI().after('>');  
    DvString base_str2 = DvString(DvUnitGetBaseSI(SIstr)).after('>');
	
	if( base_str1.empty() || base_str2.empty() ) return true; // accept numeric value w/out xref
	
    if(base_str1 != base_str2) return false;
   
	return true;

}

DvString DvObject::getSIC(size_t i){
	
	// SI_conversion may be an array e.g. (RLP vectors)

	DvObject_var SIC = this->get_xref(SI_CONVERSION);
	
    if( SIC->is_str() ) {
        if(i<arraySize()){
            return SIC->str(i);
        }
        else{
            return SIC->str();
        }
    }

	return DvString("");
}

int DvObject::getSICdim(){
	
	// SI_conversion may be an array e.g. (RLP vectors)

	DvObject_var SIC = this->get_xref(SI_CONVERSION);
	
	return SIC->arraySize();
}

double DvObject::convFactor(){
	
	return getSIC().before('>').toDouble();
}


DvString DvObject::getBaseSI(int i){
	
	DvString SIstr = this->getSIC(i);
	if(SIstr.empty()) return SIstr;
	
	return DvString(DvUnitGetBaseSI(SIstr.c_str()));
}

double DvObject::convFactorToBaseSI(){
	
	return getBaseSI().before('>').toDouble();
}



double DvObject::convFactorFrom(DvString &SIC){
	
	DvString thisConvStr = getBaseSI().before('>');
	DvString argConvStr = SIC.before('>');
	
	if(thisConvStr.empty() || argConvStr.empty() ) {
		error("unit conv, missing SI_conversion");
		return 1.;
	}
	
	double argConv = argConvStr.toDouble() / thisConvStr.toDouble();	

	return argConv;
	
}


DvString DvObject::getUnitsProduct(DvObject &arg){
	
	// units for this * arg
	if(!xref_exists(UNITS) && !(arg.xref_exists(UNITS))) return DvString("");

	DvString units = this->getXrefText(UNITS);
	units += " ";
	units += arg.getXrefText(UNITS);
	
	return units;
}

DvString DvObject::getUnitsRatio(DvObject &arg){
	
	// units for this * arg
	if(!xref_exists(UNITS) && !(arg.xref_exists(UNITS))) return DvString("");

	DvString units = this->getXrefText(UNITS);
	units += " / (";
	units += arg.getXrefText(UNITS);
	units += " )";

	return units;
}

DvString DvObject::getUnitsPower(double pow){
	
	// units for this ^ pow
	
	if(!xref_exists(UNITS) ) return DvString("");

	DvString units = "(";
	units += this->getXrefText(UNITS);
	units += ")^";
	units += pow;
	
	return units;
}

DvString DvObject::getUnitsInverse(){
		
	if(!xref_exists(UNITS) ) return DvString("");

	DvString units = "(";
	units += this->getXrefText(UNITS);
	units += ")^-1";
	
	return units;
}


DvString DvObject::getSICProduct(DvObject &arg){
 
	if(!xref_exists(SI_CONVERSION) && !(arg.xref_exists(SI_CONVERSION))) return DvString("");

	DvString null_unit = "1 > ()";

	DvString siconv1 = this->getXrefText(SI_CONVERSION);
  	DvString siconv2 = arg.getXrefText(SI_CONVERSION);
	if( siconv1.empty()) {  
    	siconv1 = null_unit;
	}
	if(siconv2.empty()) {
  		siconv2 = null_unit;
	}
  
	double value1 = siconv1.before('>').toDouble();
	double value2 = siconv2.before('>').toDouble();
  
	double value = value1 * value2;
   
	DvString siconv; 

	siconv += value;
	siconv += ">";
  	siconv += siconv1.after('>');  
  	siconv += " "; 
	siconv += siconv2.after('>'); 

	siconv =  DvUnitGetTidySI(siconv.c_str());  
  
	return siconv;

}

DvString DvObject::getSICInverse(){
 
	if(!xref_exists(SI_CONVERSION) ) return DvString("");

	DvString siconv = getXrefText(SI_CONVERSION);
	siconv =  DvUnitGetPowerSI(siconv.c_str(), -1);
    
	return siconv;

}

DvString DvObject::getSICRatio(DvObject &arg){

	if(!xref_exists(SI_CONVERSION) && !(arg.xref_exists(SI_CONVERSION))) return DvString("");

	DvString null_unit = "1 > ()";

	DvString siconv1 = this->getXrefText(SI_CONVERSION);

	if( siconv1.empty()) {  
    	siconv1 = null_unit;
	}
	
	DvString siconv2 = arg.getSICInverse();
	
	// the rest is as per product
	
	double value1 = siconv1.before('>').toDouble();
	double value2 = siconv2.before('>').toDouble();
  
	double value = value1 * value2;
   
	DvString siconv; 

	siconv += value;
	siconv += ">";
  	siconv += siconv1.after('>');  
  	siconv += " "; 
	siconv += siconv2.after('>'); 
    
	siconv =  DvUnitGetTidySI(siconv.c_str());  
  
	return siconv;
}

DvString DvObject::getSICPower(double pow){
 
	DvString siconv = getXrefText(SI_CONVERSION);
	siconv =  DvUnitGetPowerSI(siconv.c_str(), pow);
    
	return siconv;

}

bool DvObject::hasUnits(){
 
  DvString siconv = getXrefText(SI_CONVERSION);
  if (siconv.empty() ) return false;
  
  return DvUnitHasUnits(siconv.c_str());
  
}

bool DvObject::isVectDeg(){
	int nSIC = getSICdim();
	for(int i=0; i<nSIC; i++){
		if( sameUnits("1>deg", i) ) return true;
	}
	return false;
}

bool DvObject::isVectRad(){
	int nSIC = getSICdim();
	for(int i=0; i<nSIC; i++){
		if( sameUnits("1>rad", i) ) return true;
	}
	return false;
}

bool DvObject::isDeg(){
	if( sameUnits("1>deg") ) return true;
	return false;
}

bool DvObject::isRad(){
	if( sameUnits("1>rad") ) return true;
	return false;
}

bool DvObject::isAngle(){

	if( !xref_exists(SI_CONVERSION) && !xref_exists(UNITS)) return false;
	
	if( sameBaseSI("1>radian") ) { return true;}
	else if( sameBaseSI("1>degree") ) { return true;}
	else {
		DvString unit_str = getXrefText(UNITS);
		if(unit_str.contains("deg") || unit_str.contains("rad")) return true;
	}
	
	return false;
}

bool DvObject::isPhiAngle(bool isDeg, int index, double &minVal, int &kStart, int &kEnd){
	
	// used on DEPEND_i object, so rank is 1
	
	double conv2rad=1.;
	double pi = 4.*std::atan(1.);
	if(isDeg) conv2rad=pi/180.;
	
	minVal=1.e30;
    double maxVal=-1.e30;
	
	bool isPhi = false;
	
	int n_dim = dims[0];
	for(int k=0; k<n_dim; k++)
	{
		if( !isnan(asDouble(k))  )
		{
			if(kStart == 0) kStart = k+1; // slots use range (1, n)
			kEnd = k+1;
		}
	}
	
	minVal = this->min() * conv2rad;
	maxVal = this->max() * conv2rad;
	
	if( fabs(maxVal-minVal)>pi ) isPhi = true;
	
	if(!isPhi) // Let's see if the variable name gives us a hint...
	{
		DvString name=getXrefText("FIELDNAM");
		isPhi=(name.contains("phi") || name.contains("Azimut"));
	}

	if(!isPhi)
	{
		DvString representation=getXrefText("REPRESENTATION");
		isPhi=(representation=="rtp" || representation=="rlp" || representation=="rp")&&(index==3);
	}

	return isPhi;
}


DvString DvObject::angleUnitStr(){

	// unit object is an angle, else ""
	DvString strUnits = ""; 
  
	DvObject_var xref;

	if(this->xref_exists(SI_CONVERSION)) xref = this->get_xref(SI_CONVERSION);
	else if (this->xref_exists(UNITS)) xref = this->get_xref(UNITS);
	
	if( xref.is_ok() ) {
    	
		strUnits = xref->asStr();

  		if (strUnits.contains("DEG"))  strUnits = DvString("deg");
  		else if (strUnits.contains("RAD"))  strUnits = DvString("rad");
		else strUnits = DvString("");
	}
	
  	return strUnits;

}

DvObject_var DvObject::angleMod(int angMax, DvString rad_deg){
	
	DvObject_var newObj = new DvObject(*this);
	double one_eighty = 180.;
	if(rad_deg == "rad") one_eighty = 4.*std::atan(1.);

	double three_sixty = 2.*one_eighty;
	double ninety = 0.5*one_eighty;
	
	// Use dblS as it is safe as l or r value even if data is not double (though nothing changes if not double)
		
	for(size_t j=0; j<newObj->seqSize(); j++) {
		for(size_t i=0; i<newObj->arraySize(); i++){
			if(angMax == 360){
				 if(newObj->dblS(j,i) < 0) newObj->dblS(j,i) += three_sixty; //  0,360
			}
			else if(angMax == 180 ){
				 if( newObj->dblS(j,i) > one_eighty) newObj->dblS(j,i) -= three_sixty; //  -180,180
			}
			else  {
				newObj->dblS(j,i) = ninety - newObj->dblS(j,i); // 0,180 or -90,90
			}
		}
	}
	
	return newObj;
}




// ----------------  Dimension handling ---------------- 

bool DvObject::isThreeVector(){

	if(this->is_nil()) return false;
	
	// first test for frame (CSDS or CAA syntax)
	DvString framePtr = this->getFrameAttr();

	// if no frame accept object, if Frame exists test it is a vector
	if( !framePtr.empty() && !framePtr.contains("vector") ) return false;
 
    if( arraySize() != 3 ) return false; // 3 is prime, so must have 1 dim

	return true;

}

bool DvObject::isVectorXYZ(){

	if(!isThreeVector()) return false;

	if(getRep() == "xyz") return true;
	return false;
}

bool DvObject::okDims(DvObject &arg){
	// always allow operation with single number
	if(arg.arraySize() == 1) return true;
	
	return sameDims(arg);
}

bool DvObject::sameDims(DvObject &arg){
	
	if(dims.size() != arg.dims.size() ) return false;
	for(size_t i=0; i<dims.size(); i++) if(dims[i] != arg.dims[i]) return false;
	
	return true;
}

bool DvObject::conformalDims(DvObject &arg){
	
	if(arg.arraySize() == 1) return true; // multiply by scalar
	if(this->arraySize() == 1) return true; // multiply by scalar

	if(dims[dims.size()-1] == arg.dims[0] ) return true;
	
	return false;
}

bool DvObject::squareMat(){
	
	if(dims.size() == 0) return true; // scalar is square matrix
	if(dims.size() == 1 && dims[0] == 1) return true; // scalar is square matrix
	if(dims.size()!= 2) return false;
	if(dims[1] != dims[0] ) return false;
	
	return true;
}




DvObject_var DvObject::cleanObject(){

	// Removes records containing DvNaN from double type data
	// Removes records containing FILLVAL from Int type data
	// Other types are not cleaned 
	
	if(this->is_nil()) return this;

	bool objChanged = false;

	int iFill = -1;
	double dFill = DvNaN;

	if (this->xref_exists(FILL_VALUE)){
			DvObject_var fill_xref = this->get_xref(FILL_VALUE);
			iFill = fill_xref->asInt(0);
			dFill = fill_xref->asDouble(0);
	}

	if(this->is_dbl()){

		// initialise mask
		size_t seqLen = this->seqSize();
      	DvMask msk(seqLen);
		size_t arrLen = this->arraySize();
		for(size_t n=0; n<seqLen; n++){
			slice sl(n*arrLen, arrLen, 1);
			valarray <double> val = this->Ddata[sl];
			for(size_t i = 0; i<val.size(); i++){
				if( isFill(val[i], dFill) ) {
					objChanged = true;
					msk[n] = false;
					break;
				}
			}
		}
		if(objChanged) {
			DvObject_var ret = new DvObject(this); // attaches xrefs automatically
			ret->apply_mask(msk); //  will also apply mask to xrefs (after copying them)
			return ret;
		}
		else return this;
		
	}
	else if(this->is_int()){
		
		// initialise mask
		size_t seqLen = this->seqSize();
      	DvMask msk(seqLen);
		size_t arrLen = this->arraySize();
		for(size_t n=0; n<seqLen; n++){
			slice sl(n*arrLen, arrLen, 1);
			valarray <int> val = this->Idata[sl];
			for(size_t i = 0; i<val.size(); i++){
				if(val[i] == iFill) {
					objChanged = true;
					msk[n] = false;
					break;
				}
			}
		}
		if(objChanged) {
			DvObject_var ret = new DvObject(this); // attaches xrefs automatically
			ret->apply_mask(msk); //  will also apply mask to xrefs (after copying them)
			return ret;
		}		
		else return this;

	}


   // string and other input types not cleaned - return input object

   return this;

}

DvObject_var DvObject::correctColToRow(){
    // corrects data read from Column Major (Fortran) record to Row Major (C)
    // dimensions are assumed correct
    
    DvObject_var ret = new DvObject(*this); // sets size, types and metadata
    
    if(nDims() < 2) return ret;
    
    if(nDims() == 2) {
        for(size_t n=0; n<seqSize(); n++){
            for(size_t i=0; i<dims[0]; i++){
                for(size_t j=0; j<dims[0]; j++){
                    // not efficient having ifs here, but rarely used
                    if(is_dbl()) ret->dbl(n,j,i) = this->dbl(n,i,j);
                    else if(is_int()) ret->itg(n,j,i) = this->itg(n,i,j);
                    else if(is_str()) ret->str(n,j,i) = this->str(n,i,j);
                    else if(is_time()) ret->time(n,j,i) = this->time(n,i,j);
                    else  ret->event(n,j,i) = this->event(n,i,j);
                }
            }
        }
    }
    else if(nDims() == 3) {
        for(size_t n=0; n<seqSize(); n++){
            for(size_t i=0; i<dims[0]; i++){
                for(size_t j=0; j<dims[1]; j++){
                    for(size_t k=0; k<dims[2]; k++){
                        if(is_dbl()) ret->dbl(n,k,j,i) = this->dbl(n,i,j,k);
                        else if(is_int()) ret->itg(n,k,j,i) = this->itg(n,i,j,k);
                        else if(is_str()) ret->str(n,k,j,i) = this->str(n,i,j,k);
                        else if(is_time()) ret->time(n,k,j,i) = this->time(n,i,j,k);
                        else ret->event(n,k,j,i) = this->event(n,i,j,k);
                    }
                }
            }
        }
    }
    else if(nDims() == 4) {
        for(size_t n=0; n<seqSize(); n++){
            for(size_t i=0; i<dims[0]; i++){
                for(size_t j=0; j<dims[1]; j++){
                    for(size_t k=0; k<dims[2]; k++){
                        for(size_t m=0; m<dims[3]; m++){
                            if(is_dbl()) ret->dbl(n,m,k,j,i) = this->dbl(n,i,j,k,m);
                            else if(is_int())   ret->itg(n,m,k,j,i) = this->itg(n,i,j,k,m);
                            else if(is_str())   ret->str(n,m,k,j,i) = this->str(n,i,j,k,m);
                            else if(is_time())   ret->time(n,m,k,j,i) = this->time(n,i,j,k,m);
                            else    ret->event(n,m,k,j,i) = this->event(n,i,j,k,m);
                        }
                    }
                }
            }
        }
    }
    else {
        ret = new DvObject();
        ret->error(" dimension > 4 not handled for column conversion");
    }

    return ret;

}


// -------------- utilities-------------------

valarray <double> atan_yx(valarray <double> &y, valarray <double> &x){

	valarray <double> val = std::atan2(y, x);
	// 
	for(size_t i=0; i<val.size(); i++) {
		if(y[i] < 0.0) val[i] += twopi;
	}
	return val;
}

double atan_yx( double y,  double x){

	double val = std::atan2(y, x);
	if(y < 0.0) val += twopi;
	return val;
}

bool isFill(double val, double fill){
	if( isnan(fill) ) {
		if( isnan(val) ) return true;
		return false;
	}
	if( fabs(val - fill)  < DV_CLOSE_ENOUGH * fabs(val) ) return true;
	return false;
}

bool notFill(double val, double fill){
	if( isnan(fill) ) {
		if( isnan(val) ) return false;
		return true;
	}
	if( fabs(val - fill)  > DV_CLOSE_ENOUGH * fabs(val)) return true;
	return false;
}

