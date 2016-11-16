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


// This is here as DvRecord is incomplete when Dvar class is defined
template <> DvRecord Dvar<DvObject>::operator[](size_t rec) { return (*(DvObject*)p)[rec]; }

double DV_PI = 4.0*atan(1.0);

// operators
// ---------

// Single value, no metadata testing
//-----------------------------------

DvObject &DvObject::operator*=(double arg)
{
	// special operator to allow rescaling (e.g. units changes)
	
	if(Ddata.size() > 0 )  Ddata *= arg; 
	else if(Idata.size() > 0 ) Idata *= (int) arg; 
	else error("*= double not allowed"); // meaningless for events

	return *this;
	
}
DvObject &DvObject::operator*=(int arg)
{
	// special operator to allow rescaling (e.g. units changes)
	
	if(Ddata.size() > 0 )  Ddata *= (double) arg; 
	else if(Idata.size() > 0 ) Idata *= arg; 
	else error("*= int not allowed"); // meaningless for events

	return *this;
	
}

DvObject &DvObject::operator+=(double arg)
{
	// single value
	if(Ddata.size() > 0 )Ddata += arg; 
	else if(Idata.size() > 0 ) Idata += (int) arg; 
	else if(Sdata.size() > 0 ) { for(size_t i=0; i<Sdata.size(); i++) Sdata[i] += arg; }
	else if(Tdata.size() > 0 ) { for(size_t i=0; i<Tdata.size(); i++) Tdata[i] += arg; }
	
	else error("+= not defined for types"); // meaningless for events

	return *this;
} 
DvObject &DvObject::operator+=(int arg)
{
	// single value
	if(Ddata.size() > 0 )Ddata += (double) arg; 
	else if(Idata.size() > 0 ) Idata += arg; 
	else if(Sdata.size() > 0 ) { for(size_t i=0; i<Sdata.size(); i++) Sdata[i] += arg; }
	else if(Tdata.size() > 0 ) { for(size_t i=0; i<Tdata.size(); i++) Tdata[i] += (double) arg; }
	
	else error("+= not defined for types"); // meaningless for events

	return *this;
} 
DvObject &DvObject::operator+=(DvEvent arg)
{
	// special case
	if(Edata.size() > 0 ) {
		valarray <DvEvent> hold(Edata.size());
		hold = Edata;
		Edata.resize(hold.size()+1, DvEvent());
		for(size_t i=0; i<hold.size(); i++){
			Edata[i] = hold[i];
		}
		Edata[hold.size()] = arg;
		seqLen = hold.size()+1;
	}
	else error("+= not defined for type"); 

	return *this;
} 
DvObject &DvObject::operator+=(DvString arg)
{
	// special case
	if(Sdata.size() > 0 ) {
		valarray <DvString> hold(Sdata.size());
		hold = Sdata;
		Sdata.resize(hold.size()+1, DvString());
		for(size_t i=0; i<hold.size(); i++){
			Sdata[i] = hold[i];
		}
		Sdata[hold.size()] = arg;
		seqLen = hold.size()+1;
	}
	else error("+= not defined for type"); 

	return *this;
} 

DvObject &DvObject::operator-=(double arg)
{
	// single value
	if(Ddata.size() > 0 ) Ddata -= arg; 
	else if(Idata.size() > 0 ) Idata -= (int) arg; 
	else if(Tdata.size() > 0 ) { for(size_t i=0; i<Tdata.size(); i++) Tdata[i] -= arg; }
	else error("-= not defined for types"); 

	return *this;
} 
DvObject &DvObject::operator-=(int arg)
{
	// single value
	if(Ddata.size() > 0 ) Ddata -= (double) arg; 
	else if(Idata.size() > 0 ) Idata -=  arg; 
	else if(Tdata.size() > 0 ) { for(size_t i=0; i<Tdata.size(); i++) Tdata[i] -= (double) arg; }
	else error("-= not defined for types"); 

	return *this;
} 

DvObject &DvObject::operator/=(double arg)
{
	// special operator to allow rescaling (e.g. units changes)
	
	if(Ddata.size() > 0 )  Ddata /= arg; 
	else if(Idata.size() > 0 ) Idata /= (int) arg; 
//	else if(Tdata.size() > 0 ) Tdata /= arg; 
	else if(Tdata.size() > 0 ) for(size_t i=0; i<Tdata.size(); i++) Tdata[i] /= arg; 
	else error("/= not defined for types"); // meaningless for events

	return *this;
	
}
DvObject &DvObject::operator/=(int arg)
{
	// special operator to allow rescaling (e.g. units changes)
	
	if(Ddata.size() > 0 )  Ddata /= arg; 
	else if(Idata.size() > 0 ) Idata /= arg; 
	else error("/= not defined for types"); // meaningless for events

	return *this;
	
}


// double is always a special case as large data arrays are usually double, so speed matters

// +=    
//

DvObject &DvObject::operator+=(DvObject &arg)
{
    this->addErr(arg);
    
	if(seqLen == 0 || arg.seqLen == 0 ){
		error("+= empty object");
		return *this;
	}
	
	// needed for vectors, but benign for scalars or missing Frame info
	if( !sameFrame(arg) ){
		error("+= frames differ");
		return *this;
	}
	
	// test dimensions (must match or arg must be single value
	if(!okDims(arg) ){ 
		error("+= dimension mis-match");
		return *this;
	}
	
	// convert units if needed
	
	if( !sameBaseSI(arg)) {
		error("+= base units differ");
		return *this;
	}
	
	// join if needed (new object) arg_var now has no metadata
	DvObject_var arg_var;
	if(arg.seqSize() != 1) {
		arg_var = arg.Join(*this);
	}
	else arg_var = arg;
	
	if( sameUnits(arg) == false ) {
	 	// not sameconversion factor
		DvString argSIC = arg.getSIC();
		double argConv = convFactorFrom(argSIC);
		
		// copy if unchanged after join to avoid modifying input arg
		
		if(arg_var->get_id() == arg.get_id()) arg_var = new DvObject(arg);
		
		*arg_var *= argConv; // argConv is a double, no metadata testing
	}

	// Do the actual operation on valarray
	
	if(Ddata.size() > 0 ){
		if(arg_var->is_dbl()) {
			// special case of double with double for speed
			if (arg_var->totalSize()==totalSize() ) Ddata += arg_var->Ddata; // conformal sequences
			else if (arg_var->totalSize()==1) Ddata += arg_var->Ddata[0]; // single value
			else if (arg_var->seqSize() == 1 && arg_var->arraySize()==arraySize()) { // single matrix
				// NOTE: actual dimensions are not tested since array may have been unpacked
				size_t block = arraySize();
				for(size_t i=0; i<seqLen; i++){
					Ddata[slice(i*block, block, 1)] += arg_var->Ddata;  
				}
			}
            else this->error("+= dimension or length not conformal");
		}
		else{
			if (arg_var->totalSize()==totalSize() ) Ddata += arg_var->toDouble(); // conformal sequences
			else if (arg_var->totalSize()==1) Ddata += (arg_var->toDouble())[0]; // single value
			else if (arg_var->seqSize() == 1 && arg_var->arraySize()==arraySize()) { // single matrix
				size_t block = arraySize();
				for(size_t i=0; i<seqLen; i++){
					Ddata[slice(i*block, block, 1)] += arg_var->toDouble();  
				}
            }
            else this->error("+= dimension or length not conformal");
		}
	}
	
	else if(Idata.size() > 0 ){
		
		if (arg_var->totalSize()==totalSize() ) Idata += arg_var->toInt(); // conformal sequences
		else if (arg_var->totalSize()==1) Idata += (arg_var->toInt())[0]; // single value
		else if (arg_var->seqSize() == 1 && arg_var->arraySize()==arraySize()) { // single matrix
			size_t block = arraySize();
			for(size_t i=0; i<seqLen; i++){
				Idata[slice(i*block, block, 1)] += arg_var->toInt();  
			}
        }
        else this->error("+= dimension or length not conformal");
	}
	
	else if(Sdata.size() > 0 ){
		
		if (arg_var->totalSize()==totalSize() ) Sdata += arg_var->toStr(); // conformal sequences
		else if (arg_var->totalSize()==1) Sdata += (arg_var->toStr())[0]; // single value
		else if (arg_var->seqSize() == 1 && arg_var->arraySize()==arraySize()) { // single matrix
			size_t block = arraySize();
			for(size_t i=0; i<seqLen; i++){
				Sdata[slice(i*block, block, 1)] += arg_var->toStr();  
			}
        }
        else this->error("+= dimension or length not conformal");
	}
	
	else if(Tdata.size() > 0 && arg_var->is_dbl()){
		// only doubles can be added to time
		if (arg_var->totalSize()==totalSize() ) { // conformal sequences
			for(size_t i=0; i<Tdata.size(); i++) Tdata[i] += arg_var->Ddata[i];
		}
		else if (arg_var->totalSize()==1) { // single value
			for(size_t i=0; i<Tdata.size(); i++) Tdata[i] += arg_var->Ddata[0];
		} 
		else if (arg_var->seqSize() == 1 && arg_var->arraySize()==arraySize()) { // single matrix
			size_t block = arraySize();
			for(size_t i=0; i<seqLen; i++){
				for(size_t j=0; j<block; j++) Tdata[i*block + j] += arg_var->Ddata[j];  
			}
        }
        else this->error("+= dimension or length not conformal");
	}
	
	else if(Edata.size() > 0 && arg_var->is_event()){
		// add event to sequence
		if( arraySize()==1 ) { // sequence of single events (event list)
			// copy this event list to safety
			size_t thisLen = Edata.size();
			size_t argLen = arg_var->Edata.size();
			valarray<DvEvent> copyE(thisLen);
			copyE = Edata;
			Edata.resize(thisLen+argLen);
			for(size_t i=0; i<thisLen; i++){
				Edata[i] = copyE[i];
			}
			for(size_t i=0; i<argLen; i++){
				Edata[i+thisLen] = arg_var->Edata[i];
			}
        }
        else this->error("+= dimension or length not conformal");
		
	}
	else error("+= not defined for types"); // meaningless for events

	return *this;
} 

// -= 
//  

DvObject &DvObject::operator-=(DvObject &arg)
{
    this->addErr(arg);

    if(seqLen == 0 || arg.seqLen == 0 ){
		error("-= empty object");
		return *this;
	}
	
	// needed for vectors, but benign for scalars or missing Frame info
	if( !sameFrame(arg) ){
		error("-= frames differ");
		return *this;
	}
	
	// test dimensions (must match or arg must be single value
	if(!okDims(arg) ){ 
		error("-= dimension mis-match");
		return *this;
	}
	
	// convert units if needed
	
	if( !sameBaseSI(arg)) {
		error("-= base units differ");
		return *this;
	}
	
	// join if needed (new object) arg_var now has no metadata
	DvObject_var arg_var;
	if(arg.seqSize() != 1) {
		arg_var = arg.Join(*this);
	}
	else arg_var = arg;
	
	if( sameUnits(arg) == false ) {
	 	// not sameconversion factor
		DvString argSIC = arg.getSIC();
		double argConv = convFactorFrom(argSIC);
		
		// copy if unchanged after join to avoid modifying input arg
		
		if(arg_var->get_id() == arg.get_id()) arg_var = new DvObject(arg);
		
		*arg_var *= argConv; // argConv is a double, no metadata testing
	}

	
	if(Ddata.size() > 0 ){
		if(arg_var->is_dbl()) {
			// special case of double with double for speed
			if (arg_var->totalSize()==totalSize() ) Ddata -= arg_var->Ddata; // conformal sequences
			else if (arg_var->totalSize()==1) Ddata -= arg_var->Ddata[0]; // single value
			else if (arg_var->seqSize() == 1 && arg_var->arraySize()==arraySize()) { // single matrix
				// NOTE: actual dimensions are not tested since array may have been unpacked
				size_t block = arraySize();
				for(size_t i=0; i<seqLen; i++){
					Ddata[slice(i*block, block, 1)] -= arg_var->Ddata;  
				}
            }
            else this->error("-= dimension or length not conformal");
		}
		else{
			if (arg_var->totalSize()==totalSize() ) Ddata -= arg_var->toDouble(); // conformal sequences
			else if (arg_var->totalSize()==1) Ddata -= (arg_var->toDouble())[0]; // single value
			else if (arg_var->seqSize() == 1 && arg_var->arraySize()==arraySize()) { // single matrix
				size_t block = arraySize();
				for(size_t i=0; i<seqLen; i++){
					Ddata[slice(i*block, block, 1)] -= arg_var->toDouble();  
				}
            }
            else this->error("-= dimension or length not conformal");
		}
	}
	
	else if(Idata.size() > 0 ){
		
		if (arg_var->totalSize()==totalSize() ) Idata -= arg_var->toInt(); // conformal sequences
		else if (arg_var->totalSize()==1) Idata -= (arg_var->toInt())[0]; // single value
		else if (arg_var->seqSize() == 1 && arg_var->arraySize()==arraySize()) { // single matrix
			size_t block = arraySize();
			for(size_t i=0; i<seqLen; i++){
				Idata[slice(i*block, block, 1)] -= arg_var->toInt();  
			}
        }
        else this->error("-= dimension or length not conformal");
	}
		
	else if(Tdata.size() > 0 && arg_var->is_dbl()){
		// only doubles can be added to time
		if (arg_var->totalSize()==totalSize() ) { // conformal sequences
			for(size_t i=0; i<Tdata.size(); i++) Tdata[i] -= arg_var->Ddata[i];
		}
		else if (arg_var->totalSize()==1) { // single value
			for(size_t i=0; i<Tdata.size(); i++) Tdata[i] -= arg_var->Ddata[0];
		} 
		else if (arg_var->seqSize() == 1 && arg_var->arraySize()==arraySize()) { // single matrix
			size_t block = arraySize();
			for(size_t i=0; i<seqLen; i++){
				for(size_t j=0; j<block; j++) Tdata[i*block + j] -= arg_var->Ddata[j];  
			}
        }
        else this->error("-= dimension or length not conformal");
	}
	else error("-="); // meaningless for events and strings

	return *this;
} 

// *= 
//  

DvObject &DvObject::operator*=(DvObject &arg)
{	
    this->addErr(arg);
    
	if(seqLen == 0 || arg.seqLen == 0 ){
		error("*= empty object");
		return *this;
	}	
	
	if( !arg.is_dbl() || !this->is_dbl() ) {
		error("*= not defined for data type ");
		return *this;
	}
	
	// needed for vectors, but benign for scalars or missing Frame info
	if( !sameFrame(arg) ){
		error("*= frames differ");
		return *this;
	}
	
	// test dimensions give result with same dims
	if( !conformalDims(arg) || !arg.squareMat() ){ 
		error("*= dimension not conformal");
		return *this;
	}
	
	// join if needed (new object) arg_var now has no metadata
	DvObject_var arg_var;
	if(arg.seqSize() != 1) {
		arg_var = arg.Join(*this);
	}
	else arg_var = arg;
	
    // In the event that this->seqLen == 1 and arg->seqLen >1 join ensures only
    // one value is created for arg_var, so following is safe
    
	if (arg_var->totalSize()==1) {
		(*this) *= arg_var->Ddata[0]; // single value
	}
	else if(arg_var->arraySize() == 1){
		for(size_t i=0; i<seqSize(); i++) (*this)[i] *= arg_var->Ddata[i]; // single values at each record
	}
	else if (arg_var->seqSize() == 1 && arg_var->arraySize()==arraySize()) { // single matrix

		size_t common = arg_var->dims[0];
		size_t thisArrLen = this->arraySize();
		size_t argArrLen = arg_var->arraySize();
		size_t thisStride = thisArrLen/common;
		size_t argStride = argArrLen/common;
		
		size_t newArraySize = thisStride * argStride;
		
		slice sl_val(0, argArrLen, 1);
	
		for(size_t r=0; r<seqSize(); r++) {
			valarray <double> res(0.0, newArraySize);
						
			slice sl(r*thisArrLen, thisArrLen, 1);
		
			valarray <double> v1 = this->Ddata[sl];
			valarray <double> v2 = arg_var->Ddata[sl_val];

			// this should work for any dimension and rank 		
			for(size_t i=0; i<thisStride; i++){
				for(size_t j=0; j<argStride; j++){
					for(size_t p=0; p<common; p++){ // sum down inner dims
						res[i*argStride+j] += v1[i*common+p] * v2[p*argStride+j]; 
					}
				}
			}
			// copy result into valarray
			Ddata[sl] = res;
		}
	}
	else  { //  matrix sequences
		size_t common = arg_var->dims[0];
		size_t thisArrLen = this->arraySize();
		size_t argArrLen = arg_var->arraySize();
		size_t thisStride = thisArrLen/common;
		size_t argStride = argArrLen/common;
		
		size_t newArraySize = thisStride * argStride;
			
		for(size_t r=0; r<seqSize(); r++) {
			valarray <double> res(0.0, newArraySize);		// must zero at each record				
			slice sl(r*thisArrLen, thisArrLen, 1);
			slice sl_val(r*argArrLen, argArrLen, 1);
			slice sl_res(r*newArraySize, newArraySize, 1);
		
			valarray <double> v1 = this->Ddata[sl];
			valarray <double> v2 = arg_var->Ddata[sl_val];

			// this should work for any dimension and rank 		
			for(size_t i=0; i<thisStride; i++){
				for(size_t j=0; j<argStride; j++){
					for(size_t p=0; p<common; p++){ // sum down inner dims
						res[i*argStride+j] += v1[i*common+p] * v2[p*argStride+j]; 
					}
				}
			}
		
			// copy result into valarray
			Ddata[sl] = res;
		}
	}
		
	// convert units 
	
	DvString SIC = getSICProduct(arg);
	if(!SIC.empty()) change_xref(SI_CONVERSION, SIC);

	DvString unit = getUnitsProduct(arg);
	if(!unit.empty()) change_xref(UNITS, unit);

	// do frame
	setAttrsProduct(arg);
	
	return *this;
		
} 

// 
//  

DvObject &DvObject::operator/=(DvObject &arg)
{
    this->addErr(arg);
    
	if(seqLen == 0 || arg.seqLen == 0 ){
		error("/= empty object");
		return *this;
	}	
	
	if( !arg.is_dbl() || !this->is_dbl() ) {
		error("/= not defined for data type ");
		return *this;
	}
	
	
	// test dimensions give result with same dims
	if( !okDims(arg) ){ 
		error("/= dimension not same or scalar divisor");
		return *this;
	}
	
	// join if needed (new object) arg_var now has no metadata
	DvObject_var arg_var;
	if(arg.seqSize() != 1) {
		arg_var = arg.Join(*this);
	}
	else arg_var = arg;
	
	if (arg_var->totalSize() == this->totalSize()) {
		this->Ddata /= arg_var->Ddata; // single value
	}
	else if (arg_var->totalSize() == 1) {
		(*this) /= arg_var->Ddata[0]; // single value
	}
	else if(arg_var->arraySize() == 1){
		for(size_t i=0; i<seqSize(); i++) (*this)[i] /= arg_var->Ddata[i]; // single values at each record
	}
	else if (arg_var->seqSize() == 1 && arg_var->arraySize()==arraySize()) { // single matrix

		// NOTE: actual dimensions are not tested since array may have been unpacked
		size_t block = arraySize();
		for(size_t i=0; i<seqSize(); i++){
			Ddata[slice(i*block, block, 1)] /= arg_var->Ddata;  
		}
	}
		
	// convert units 
	
	DvString SIC = getSICRatio(arg);
	if(!SIC.empty()) change_xref(SI_CONVERSION, SIC);

	DvString unit = getUnitsRatio(arg);
	if(!unit.empty()) change_xref(UNITS, unit);

	return *this;
		
} 


//  + 
//  

DvObject *DvObject::operator+(DvObject &arg)
{
	if(seqLen == 0 || arg.seqLen == 0 ){
		DvObject *empty = new DvObject();
		empty->error("+ empty object");
		return empty;
	}

    if(this->seqSize() == 1 && arg.seqSize() > 1){
        
        // make 1st object into sequence
        DvObject_var thisSeq = new DvObject(this, arg.seqSize());
        DvObject_var tt = arg.getDep0();
        if(tt.is_ok()) thisSeq->change_xref(DEPEND_0, tt);
        DvObject *res = *thisSeq + arg;
        
        return res;
    }

    
	// needed for vectors, but benign for scalars or missing Frame info
	if( !sameFrame(arg) ){
		DvObject *empty = new DvObject();
		empty->error("+ frames differ");
		return empty;
	}
	
	// test dimensions (must match or arg must be single value
	if(!okDims(arg) ){ 
		DvObject *empty = new DvObject();
		empty->error("+ dimension mis-match");
		return empty;
	}
	
	// convert units if needed
	
	if( !sameBaseSI(arg)) {
		DvObject *empty = new DvObject();
		empty->error("+ base units differ");
		return empty;
	}

    DvObject *res = new DvObject(this);    // copies xrefs

    
	// join if needed (new object) arg_var now has no metadata
	DvObject_var arg_var;
	if(arg.seqSize() != 1 ) {
		arg_var = arg.Join(*this);
	}
	else arg_var = arg;
	
	if( sameUnits(arg) == false ) {
	 	// not sameconversion factor
		DvString argSIC = arg.getSIC();
		double argConv = convFactorFrom(argSIC);
		
		// copy if unchanged after join to avoid modifying input arg
		
		if(arg_var->get_id() == arg.get_id()) arg_var = new DvObject(arg);
		
		*arg_var *= argConv; // argConv is a double, no metadata testing
	}
	if(Ddata.size() > 0 ){

		if(arg_var->is_dbl()) {
			// special case of double with double for speed
			if (arg_var->totalSize()==totalSize() ) res->Ddata += arg_var->Ddata; // conformal sequences
			else if (arg_var->totalSize()==1) res->Ddata += arg_var->Ddata[0]; // single value
			else if (arg_var->seqSize() == 1 && arg_var->arraySize()==arraySize()) { // single matrix
				// NOTE: actual dimensions are not tested since array may have been unpacked
				size_t block = arraySize();
				for(size_t i=0; i<seqLen; i++){
					res->Ddata[slice(i*block, block, 1)] += arg_var->Ddata;  
				}
			}
		}
		else{
			if (arg_var->totalSize()==totalSize() ) res->Ddata += arg_var->toDouble(); // conformal sequences
			else if (arg_var->totalSize()==1) res->Ddata += (arg_var->toDouble())[0]; // single value
			else if (arg_var->seqSize() == 1 && arg_var->arraySize()==arraySize()) { // single matrix
				size_t block = arraySize();
				for(size_t i=0; i<seqLen; i++){
					res->Ddata[slice(i*block, block, 1)] += arg_var->toDouble();  
				}
			}
		}
	}
	
	else if(Idata.size() > 0 ){
		
		if (arg_var->totalSize()==totalSize() ) res->Idata += arg_var->toInt(); // conformal sequences
		else if (arg_var->totalSize()==1) res->Idata += (arg_var->toInt())[0]; // single value
		else if (arg_var->seqSize() == 1 && arg_var->arraySize()==arraySize()) { // single matrix
			size_t block = arraySize();
			for(size_t i=0; i<seqLen; i++){
				res->Idata[slice(i*block, block, 1)] += arg_var->toInt();  
			}
		}
	}
	
	else if(Sdata.size() > 0 ){
		
		if (arg_var->totalSize()==totalSize() ) res->Sdata += arg_var->toStr(); // conformal sequences
		else if (arg_var->totalSize()==1) res->Sdata += (arg_var->toStr())[0]; // single value
		else if (arg_var->seqSize() == 1 && arg_var->arraySize()==arraySize()) { // single matrix
			size_t block = arraySize();
			for(size_t i=0; i<seqLen; i++){
				res->Sdata[slice(i*block, block, 1)] += arg_var->toStr();  
			}
		}
	}
	
	else if(Tdata.size() > 0 && arg_var->is_dbl()){
		// only doubles can be added to time
		if (arg_var->totalSize()==totalSize() ) { // conformal sequences
			for(size_t i=0; i<Tdata.size(); i++) res->Tdata[i] += arg_var->Ddata[i];
		}
        else if (arg_var->totalSize()==1) { // single value
            for(size_t i=0; i<Tdata.size(); i++) res->Tdata[i] += arg_var->Ddata[0];
        } 
        else if (totalSize() == 1 && arg_var->totalSize() > 1 ) { // single value and seq
            
            for(size_t i=0; i<Tdata.size(); i++) res->Tdata[i] += arg_var->Ddata[0];
        } 
		else if (arg_var->seqSize() == 1 && arg_var->arraySize()==arraySize()) { // single matrix
			size_t block = arraySize();
			for(size_t i=0; i<seqLen; i++){
				for(size_t j=0; j<block; j++) res->Tdata[i*block + j] += arg_var->Ddata[j];  
			}
		}
	}
	else res->error("+"); // meaningless for events

	return res;
} 

//  -
//  

DvObject *DvObject::operator-(DvObject &arg)
{	
	if(seqLen == 0 || arg.seqLen == 0 ){
		DvObject *empty = new DvObject();
		empty->error("- empty object");
		return empty;
	}
    if(this->seqSize() == 1 && arg.seqSize() > 1){
        
        // make 1st object into sequence
        DvObject_var thisSeq = new DvObject(this, arg.seqSize());
        DvObject_var tt = arg.getDep0();
        if(tt.is_ok()) thisSeq->change_xref(DEPEND_0, tt);
        DvObject *res = *thisSeq - arg;
        
        return res;
    }
	
	// needed for vectors, but benign for scalars or missing Frame info
	if( !sameFrame(arg) ){
		DvObject *empty = new DvObject();
		empty->error("- frames differ");
		return empty;
	}
	
	// test dimensions (must match or arg must be single value
	if(!okDims(arg) ){ 
		DvObject *empty = new DvObject();
		empty->error("- dimension mis-match");
		return empty;
	}
	
	// convert units if needed
	
	if( !sameBaseSI(arg)) { 
		DvObject *empty = new DvObject();
		empty->error("- base units differ");
		return empty;
	}

	
	// join if needed (new object) arg_var now has no metadata
	DvObject_var arg_var;
	if(arg.seqSize() != 1) {
		arg_var = arg.Join(*this);
	}
	else arg_var = arg;
	
	// create new target object
	DvObject *res;
	if(Tdata.size() > 0 && arg_var->is_time()) res = new DvObject((double)0.0, this->Dims(), this->seqSize()); // special case
	else res = new DvObject(this);

	if( sameUnits(arg) == false ) {
	 	// not sameconversion factor
		DvString argSIC = arg.getSIC();
		double argConv = convFactorFrom(argSIC);
		// copy if unchanged after join to avoid modifying input arg
		
		if(arg_var->get_id() == arg.get_id()) arg_var = new DvObject(arg);
		
		*arg_var *= argConv; // argConv is a double, no metadata testing
	}
	
	if(Ddata.size() > 0 ){
		if(arg_var->is_dbl()) {
			// special case of double with double for speed
			if (arg_var->totalSize()==totalSize() ) res->Ddata -= arg_var->Ddata; // conformal sequences
			else if (arg_var->totalSize()==1) res->Ddata -= arg_var->Ddata[0]; // single value
			else if (arg_var->seqSize() == 1 && arg_var->arraySize()==arraySize()) { // single matrix
				// NOTE: actual dimensions are not tested since array may have been unpacked
				size_t block = arraySize();
				for(size_t i=0; i<seqLen; i++){
					res->Ddata[slice(i*block, block, 1)] -= arg_var->Ddata;  
				}
			}
		}
		else{
			if (arg_var->totalSize()==totalSize() ) res->Ddata -= arg_var->toDouble(); // conformal sequences
			else if (arg_var->totalSize()==1) res->Ddata -= (arg_var->toDouble())[0]; // single value
			else if (arg_var->seqSize() == 1 && arg_var->arraySize()==arraySize()) { // single matrix
				size_t block = arraySize();
				for(size_t i=0; i<seqLen; i++){
					res->Ddata[slice(i*block, block, 1)] -= arg_var->toDouble();  
				}
			}
		}
	}
	
	else if(Idata.size() > 0 ){
		
		if (arg_var->totalSize()==totalSize() ) res->Idata -= arg_var->toInt(); // conformal sequences
		else if (arg_var->totalSize()==1) res->Idata -= (arg_var->toInt())[0]; // single value
		else if (arg_var->seqSize() == 1 && arg_var->arraySize()==arraySize()) { // single matrix
			size_t block = arraySize();
			for(size_t i=0; i<seqLen; i++){
				res->Idata[slice(i*block, block, 1)] -= arg_var->toInt();  
			}
		}
	}
	
	else if(Tdata.size() > 0 && arg_var->is_dbl()){
		// only doubles can be added to time
		if (arg_var->totalSize()==totalSize() ) { // conformal sequences
			for(size_t i=0; i<Tdata.size(); i++) res->Tdata[i] -= arg_var->Ddata[i];
		}
		else if (arg_var->totalSize()==1) { // single value
			for(size_t i=0; i<Tdata.size(); i++) res->Tdata[i] -= arg_var->Ddata[0];
		} 
		else if (arg_var->seqSize() == 1 && arg_var->arraySize()==arraySize()) { // single matrix
			size_t block = arraySize();
			for(size_t i=0; i<seqLen; i++){
				for(size_t j=0; j<block; j++) res->Tdata[i*block + j] -= arg_var->Ddata[j];  
			}
		}
	}

	else if(Tdata.size() > 0 && arg_var->is_time()){
		// difference between times is double seconds
		if (arg_var->totalSize()==totalSize() ) { // conformal sequences
			for(size_t i=0; i<Tdata.size(); i++) res->Ddata[i] = Tdata[i] - arg_var->Tdata[i];
		}
		else if (arg_var->totalSize()==1) { // single value
			for(size_t i=0; i<Tdata.size(); i++) res->Ddata[i] = Tdata[i] - arg_var->Tdata[0];
		} 
		else if (arg_var->seqSize() == 1 && arg_var->arraySize()==arraySize()) { // single matrix
			size_t block = arraySize();
			for(size_t i=0; i<seqLen; i++){
				for(size_t j=0; j<block; j++) res->Ddata[i*block + j] = Tdata[i*block + j] - arg_var->Tdata[j];  
			}
		}
	}

	else res->error("- data type not supported"); // meaningless for events and strings

	return res;
} 

//  

DvObject *DvObject::operator*(DvObject &arg)
{	

	if( !arg.is_dbl() || !this->is_dbl() ) {
	
		//  handle non-double data types
	
		if( arg.is_int() && this->is_int() ) return this->multiplyInt(arg);

		DvObject_var argVar(arg);
		DvObject_var thisVar(*this);
	
		if( !arg.is_dbl() ){
			argVar = arg.convertToDbl();
		}	

		if( !this->is_dbl() ){
			thisVar = this->convertToDbl();
		}	
	
		return *thisVar * *argVar;

	}

	if(seqLen == 0 || arg.seqLen == 0 ){
		DvObject *empty = new DvObject();
		empty->error("* empty object");
		return empty;
	}	

    if(this->seqSize() == 1 && arg.seqSize() > 1){
        
        // make 1st object into sequence
        DvObject_var thisSeq = new DvObject(this, arg.seqSize());
        DvObject_var tt = arg.getDep0();
        if(tt.is_ok()) thisSeq->change_xref(DEPEND_0, tt);
        DvObject *res = *thisSeq * arg;
        
        return res;
    }
    

	// needed for vectors, but benign for scalars or missing Frame info
	if( !sameFrame(arg) ){
		DvObject *empty = new DvObject();
		empty->error("* frames differ");
		return empty;
	}
	
	// test dimensions give result with same dims
	if( !conformalDims(arg) ){ 
		DvObject *empty = new DvObject();
		empty->error("* dimension not commensurate");
		return empty;
	}


	if(this->arraySize() == 1 && arg.arraySize() > 1){
		// special case of multiply scalar sequence by array
		// invert order so that array gets converted to sequence of arrays
	
		DvObject *res = arg * *this;
		
		return res; 
	}
	
	// join if needed (new object) arg_var now has no metadata
	DvObject_var arg_var;
	if(arg.seqSize() != 1) {
		arg_var = arg.Join(*this);
	}
	else arg_var = arg;

	DvObject *res_ptr;
	
	if (arg_var->totalSize()==1) {
		res_ptr = new DvObject(this);		
		(*res_ptr) *= arg_var->Ddata[0]; // single value
	}
	else if(arg_var->arraySize() == 1){
		res_ptr = new DvObject(this);
		for(size_t i=0; i<seqSize(); i++) (*res_ptr)[i] *= arg_var->Ddata[i]; // single values at each record
	}
	else if(this->arraySize() == 1){
		res_ptr = new DvObject(arg_var);
		for(size_t i=0; i<arg_var->seqSize(); i++) (*res_ptr)[i] *= this->Ddata[i]; // single values at each record
	}
	else if (arg_var->seqSize() == 1 ) { // single matrix
	
		vector <size_t> dimObj; // dimensions of new array (may be empty for scalar result)
		for(size_t i=0; i<this->nDims()-1; i++) dimObj.push_back(this->dims[i]);
		for(size_t i=1; i<arg_var->nDims(); i++) dimObj.push_back(arg_var->dims[i]);
		
		// new sequence with correct dimensions
		res_ptr = new DvObject(0.0, dimObj, this->seqSize());

		size_t common = arg_var->dims[0];
		size_t thisArrLen = this->arraySize();
		size_t argArrLen = arg_var->arraySize();
		size_t thisStride = thisArrLen/common;
		size_t argStride = argArrLen/common;
		
		size_t newArraySize = thisStride * argStride;
		
		slice sl_val(0, argArrLen, 1);
	
		for(size_t r=0; r<seqSize(); r++) {
			valarray <double> res(0.0, newArraySize);
						
			slice sl(r*thisArrLen, thisArrLen, 1);
			slice sl_res(r*newArraySize, newArraySize, 1);
		
			valarray <double> v1 = this->Ddata[sl];
			valarray <double> v2 = arg_var->Ddata[sl_val];

			// this should work for any dimension and rank 		
			for(size_t i=0; i<thisStride; i++){
				for(size_t j=0; j<argStride; j++){
					for(size_t p=0; p<common; p++){ // sum down inner dims
						res[i*argStride+j] += v1[i*common+p] * v2[p*argStride+j]; 
					}
				}
			}
			// copy result into valarray
			res_ptr->Ddata[sl_res] = res;
		}
	}
	else  { //  matrix sequences
		vector <size_t> dimObj; // dimensions of new array (may be empty for scalar result)
		for(size_t i=0; i<this->nDims()-1; i++) dimObj.push_back(this->dims[i]);
		for(size_t i=1; i<arg_var->nDims(); i++) dimObj.push_back(arg_var->dims[i]);
		// new sequence with correct dimensions
		res_ptr = new DvObject(0.0, dimObj, this->seqSize());

		size_t common = arg_var->dims[0];
		size_t thisArrLen = this->arraySize();
		size_t argArrLen = arg_var->arraySize();
		size_t thisStride = thisArrLen/common;
		size_t argStride = argArrLen/common;
		
		size_t newArraySize = thisStride * argStride;
			
		for(size_t r=0; r<seqSize(); r++) {
			valarray <double> res(0.0, newArraySize);		// must zero at each record				
			slice sl(r*thisArrLen, thisArrLen, 1);
			slice sl_val(r*argArrLen, argArrLen, 1);
			slice sl_res(r*newArraySize, newArraySize, 1);
		
			valarray <double> v1 = this->Ddata[sl];
			valarray <double> v2 = arg_var->Ddata[sl_val];

			// this should work for any dimension and rank 		
			for(size_t i=0; i<thisStride; i++){
				for(size_t j=0; j<argStride; j++){
					for(size_t p=0; p<common; p++){ // sum down inner dims
						res[i*argStride+j] += v1[i*common+p] * v2[p*argStride+j]; 
					}
				}
			}
		
			// copy result into valarray
			res_ptr->Ddata[sl_res] = res;
		}
	}
	
	res_ptr->copy_xrefs_from(*this);
	
	// convert units 
	
	DvString SIC = getSICProduct(arg);
	if(!SIC.empty()) res_ptr->change_xref(SI_CONVERSION, SIC);

	DvString unit = getUnitsProduct(arg);
	if(!unit.empty()) res_ptr->change_xref(UNITS, unit);

	// set frame (and toFrame) and DEPEND_i
	
	res_ptr->setAttrsProduct(*arg_var); // use joined arg in case subset leaves only one record

	// LABLAXIS
	DvString labl_s = this->getXrefText("LABLAXIS");
	if(arg.xref_exists("LABLAXIS")){
		labl_s += " x ";
		labl_s += arg.getXrefText("LABLAXIS");
	}
	else if(arg.seqSize() == 1){
		labl_s += " x ";
		labl_s += arg.asStr(0);
	}
	res_ptr->change_xref("LABLAXIS", labl_s);

	return res_ptr;
} 

DvObject *DvObject::multiplyInt(DvObject &arg)
{	
	// used if both objects are int

	if(seqLen == 0 || arg.seqLen == 0 ){
		DvObject *empty = new DvObject();
		empty->error("* empty object");
		return empty;
	}	

    
    if(this->seqSize() == 1 && arg.seqSize() > 1){
        // make 1st object into sequence
        DvObject_var thisSeq = new DvObject(this, arg.seqSize());
        DvObject_var tt = arg.getDep0();
        if(tt.is_ok()) thisSeq->change_xref(DEPEND_0, tt);
        DvObject *res = thisSeq->multiplyInt(arg);
        
        return res;
    }
    
	// needed for vectors, but benign for scalars or missing Frame info
	if( !sameFrame(arg) ){
		DvObject *empty = new DvObject();
		empty->error("* frames differ");
		return empty;
	}
	
	// test dimensions give result with same dims
	if( !conformalDims(arg) ){ 
		DvObject *empty = new DvObject();
		empty->error("* dimension not commensurate");
		return empty;
	}
	

	// join if needed (new object) arg_var now has no metadata
	DvObject_var arg_var;
	if(arg.seqSize() != 1) {
		arg_var = arg.Join(*this); // drop arg xrefs, not used
	}
	else arg_var = arg;
	
	DvObject *res_ptr; 
	
	if (arg_var->totalSize()==1) {
		res_ptr = new DvObject(this);		
		(*res_ptr) *= arg_var->Idata[0]; // single value
	}
	else if(arg_var->arraySize() == 1){
		res_ptr = new DvObject(this);
		for(size_t i=0; i<seqSize(); i++) (*res_ptr)[i] *= arg_var->Idata[i]; // single values at each record
	}
	else if(this->arraySize() == 1){
		res_ptr = new DvObject(arg_var);
		for(size_t i=0; i<arg_var->seqSize(); i++) (*res_ptr)[i] *= this->Idata[i]; // single values at each record
	}
	else if (arg_var->seqSize() == 1 ) { // single matrix
	
		vector <size_t> dimObj; // dimensions of new array (may be empty for scalar result)
		for(size_t i=0; i<this->nDims()-1; i++) dimObj.push_back(this->dims[i]);
		for(size_t i=1; i<arg_var->nDims(); i++) dimObj.push_back(arg_var->dims[i]);
		
		// new sequence with correct dimensions
		res_ptr = new DvObject((int)0, dimObj, this->seqSize());

		size_t common = arg_var->dims[0];
		size_t thisArrLen = this->arraySize();
		size_t argArrLen = arg_var->arraySize();
		size_t thisStride = thisArrLen/common;
		size_t argStride = argArrLen/common;
		
		size_t newArraySize = thisStride * argStride;
		
		slice sl_val(0, argArrLen, 1);
	
		for(size_t r=0; r<seqSize(); r++) {
			valarray <int> res(0, newArraySize);
						
			slice sl(r*thisArrLen, thisArrLen, 1);
			slice sl_res(r*newArraySize, newArraySize, 1);
		
			valarray <int> v1 = this->Idata[sl];
			valarray <int> v2 = arg_var->Idata[sl_val];

			// this should work for any dimension and rank 		
			for(size_t i=0; i<thisStride; i++){
				for(size_t j=0; j<argStride; j++){
					for(size_t p=0; p<common; p++){ // sum down inner dims
						res[i*argStride+j] += v1[i*common+p] * v2[p*argStride+j]; 
					}
				}
			}
			// copy result into valarray
			res_ptr->Idata[sl_res] = res;
		}
	}
	else  { //  matrix sequences
		vector <size_t> dimObj; // dimensions of new array (may be empty for scalar result)
		for(size_t i=0; i<this->nDims()-1; i++) dimObj.push_back(this->dims[i]);
		for(size_t i=1; i<arg_var->nDims(); i++) dimObj.push_back(arg_var->dims[i]);
		// new sequence with correct dimensions
		res_ptr = new DvObject((int)0, dimObj, this->seqSize());

		size_t common = arg_var->dims[0];
		size_t thisArrLen = this->arraySize();
		size_t argArrLen = arg_var->arraySize();
		size_t thisStride = thisArrLen/common;
		size_t argStride = argArrLen/common;
		
		size_t newArraySize = thisStride * argStride;
			
		for(size_t r=0; r<seqSize(); r++) {
			valarray <int> res(0, newArraySize);		// must zero at each record				
			slice sl(r*thisArrLen, thisArrLen, 1);
			slice sl_val(r*argArrLen, argArrLen, 1);
			slice sl_res(r*newArraySize, newArraySize, 1);
		
			valarray <int> v1 = this->Idata[sl];
			valarray <int> v2 = arg_var->Idata[sl_val];

			// this should work for any dimension and rank 		
			for(size_t i=0; i<thisStride; i++){
				for(size_t j=0; j<argStride; j++){
					for(size_t p=0; p<common; p++){ // sum down inner dims
						res[i*argStride+j] += v1[i*common+p] * v2[p*argStride+j]; 
					}
				}
			}
		
			// copy result into valarray
			res_ptr->Idata[sl_res] = res;
		}
	}
	
	res_ptr->copy_xrefs_from(*this);
	
	// convert units 
	
	DvString SIC = getSICProduct(arg);
	if(!SIC.empty()) res_ptr->change_xref(SI_CONVERSION, SIC);

	DvString unit = getUnitsProduct(arg);
	if(!unit.empty()) res_ptr->change_xref(UNITS, unit);

	// set frame (and toFrame) and DEPEND_i
	
	res_ptr->setAttrsProduct(arg);

	// LABLAXIS
	DvString labl_s = this->getXrefText("LABLAXIS");
	if(arg.xref_exists("LABLAXIS")){
		labl_s += " x ";
		labl_s += arg.getXrefText("LABLAXIS");
	}
	else if(arg.seqSize() == 1){
		labl_s += " x ";
		labl_s += arg.asStr(0);
	}
	res_ptr->change_xref("LABLAXIS", labl_s);

	
	return res_ptr;


}


//  

DvObject *DvObject::operator/(DvObject &arg)
{
	if(seqLen == 0 || arg.seqLen == 0 ){
		DvObject *empty = new DvObject();
		empty->error("/ empty object");
		return empty;
	}	
    
    if(this->seqSize() == 1 && arg.seqSize() > 1){
        // make 1st object into sequence
        DvObject_var thisSeq = new DvObject(this, arg.seqSize());
        DvObject_var tt = arg.getDep0();
        if(tt.is_ok()) thisSeq->change_xref(DEPEND_0, tt);
        DvObject *res = *thisSeq / arg;
        
        return res;
    }
	
	if( !arg.is_dbl() || !this->is_dbl() ) {
		DvObject *empty = new DvObject();
		empty->error("/ not defined for data type");
		return empty;
	}
	
	
	// test dimensions give result with same dims
	if( !okDims(arg) ){ 
		DvObject *empty = new DvObject();
		empty->error("/ dimension not same or not scalar divisor");
		return empty;
	}


	DvObject *res = new DvObject(this);
	
	// join if needed (new object) arg_var now has no metadata
	DvObject_var arg_var;
	if(arg.seqSize() != 1) {
		arg_var = arg.Join(*this);
	}
	else arg_var = arg;
	
	if (arg_var->totalSize() == this->totalSize()) {
		res->Ddata /= arg_var->Ddata; // single value
	}
	else if (arg_var->totalSize() == 1) {
		(*res) /= arg_var->Ddata[0]; // single value
	}
	else if(arg_var->arraySize() == 1){
		for(size_t i=0; i<seqSize(); i++) (*res)[i] /= arg_var->Ddata[i]; // single values at each record
	}
	else if (arg_var->seqSize() == 1 && arg_var->arraySize()==arraySize()) { // single matrix

		// NOTE: actual dimensions are not tested since array may have been unpacked
		int block = arraySize();
		for(size_t i=0; i<seqSize(); i++){
			res->Ddata[slice(i*block, block, 1)] /= arg_var->Ddata;  
		}
	}
		
	// convert units 
	
	DvString SIC = getSICRatio(arg);
	if(!SIC.empty()) res->change_xref(SI_CONVERSION, SIC);

	DvString unit = getUnitsRatio(arg);
	if(!unit.empty()) res->change_xref(UNITS, unit);

	// LABLAXIS
	DvString labl_s = this->getXrefText( "LABLAXIS" );
	if(arg.xref_exists("LABLAXIS")){
		labl_s += "/";
		labl_s += arg.getXrefText( "LABLAXIS" );
	}
	else if(arg.seqSize() == 1){
		
		labl_s += " / ";
		labl_s += arg.asStr(0);
	}
	res->change_xref( "LABLAXIS", labl_s);

	return res;
		
} 

// Conversions

valarray <double>  DvObject::toDouble(double safeDefault){
	size_t valarraySize = totalSize();
	
	if(valarraySize == 0) {
		valarray <double> res(safeDefault, 1);
		return res;
	}
	valarray <double> res(safeDefault, valarraySize);
	
	if(Ddata.size() > 0 ) {
		return Ddata;
	}
	else if(Idata.size() > 0 ) {
		for(size_t i=0; i<Idata.size(); i++) res[i] = (double) Idata[i];
	}
	else if(Sdata.size() > 0 ) {
		for(size_t i=0; i<Sdata.size(); i++) res[i] = Sdata[i].toDouble();
	}
	else if(Tdata.size() > 0 ) {
		for(size_t i=0; i<Tdata.size(); i++) res[i] = Tdata[i].getEpoch2000();
	}
	else if(Edata.size() > 0 ) {
		for(size_t i=0; i<Edata.size(); i++) res[i] = Edata[i].duration();
	}
	return res;
}

valarray <int>  DvObject::toInt(int safeDefault){
	size_t valarraySize = totalSize();
	
	if(valarraySize == 0) {
		valarray <int> res(safeDefault, 1);
		return res;
	}
	valarray <int> res(safeDefault, valarraySize);
	
	if(Ddata.size() > 0 ) {
		for(size_t i=0; i<Ddata.size(); i++) res[i] = (int) Ddata[i];
	}
	else if(Idata.size() > 0 ) {
		return Idata;
	}
	else if(Sdata.size() > 0 ) {
		for(size_t i=0; i<Sdata.size(); i++) res[i] = Sdata[i].toInt();
	}
	else if(Tdata.size() > 0 ) {
		for(size_t i=0; i<Tdata.size(); i++) res[i] = (int) Tdata[i].getEpoch2000();
	}
	else if(Edata.size() > 0 ) {
		for(size_t i=0; i<Edata.size(); i++) res[i] = (int) Edata[i].duration();
	}
	return res;
}

valarray <DvString>  DvObject::toStr(DvString safeDefault){
	size_t valarraySize = totalSize();
	
	if(valarraySize == 0) {
		valarray <DvString> res(safeDefault, 1);
		return res;
	}
	
	valarray <DvString> res(safeDefault, valarraySize);
	if(Ddata.size() > 0 ) {
		for(size_t i=0; i<Ddata.size(); i++) res[i].set(Ddata[i]);
	}
	else if(Idata.size() > 0 ) {
		for(size_t i=0; i<Idata.size(); i++) res[i].set(Idata[i]);
	}
	else if(Sdata.size() > 0 ) {
		return Sdata;
	}
	else if(Tdata.size() > 0 ) {
		for(size_t i=0; i<Tdata.size(); i++) res[i] = Tdata[i].getISOstring().c_str();
	}
	else if(Edata.size() > 0 ) {
		for(size_t i=0; i<Edata.size(); i++) res[i] = Edata[i].getISOstring().c_str();
	}
	return res;
}

valarray <DvTime>  DvObject::toTime(DvTime safeDefault){
	size_t valarraySize = totalSize();
	
	if(valarraySize == 0) {
		valarray <DvTime> res(safeDefault, 1);
		return res;
	}
	valarray <DvTime> res(safeDefault, valarraySize);
	
	if(Ddata.size() > 0 ) {
		for(size_t i=0; i<Ddata.size(); i++) res[i].setFromEpoch2000(Ddata[i]);
	}
	if(Idata.size() > 0 ) {
		for(size_t i=0; i<Idata.size(); i++) res[i].setFromEpoch2000((double)Idata[i]);
	}
	else if(Sdata.size() > 0 ) {
		for(size_t i=0; i<Sdata.size(); i++) res[i].setFromISOstring(Sdata[i]);
	}
	else if(Tdata.size() > 0 ) {
		return Tdata;
	}
	else if(Edata.size() > 0 ) {
		for(size_t i=0; i<Edata.size(); i++) res[i] = Edata[i].start();
	}
	return res;
}


double DvObject::asDouble(size_t posn, double safeDefault){
	double res(safeDefault);
	
	if(Ddata.size() > posn ) {
		return Ddata[posn];
	}
	else if(Idata.size() > posn ) {
		res = (double) Idata[posn];
	}
	else if(Sdata.size() > posn ) {
		res = Sdata[posn].toDouble();
	}
	else if(Tdata.size() > posn ) {
		res = Tdata[posn].getEpoch2000();
	}
	else if(Edata.size() > posn ) {
		res = Edata[posn].duration();
	}
	return res;
}

int DvObject::asInt(size_t posn, int safeDefault){
	int res(safeDefault);
	
	if(Ddata.size() > posn ) {
		res = (int) Ddata[posn];
	}
	else if(Idata.size() > posn ) {
		return Idata[posn];
	}
	else if(Sdata.size() > posn ) {
		res = Sdata[posn].toInt();
	}
	else if(Tdata.size() > posn ) {
		res = (int) Tdata[posn].getEpoch2000();
	}
	else if(Edata.size() > posn ) {
		res = (int) Edata[posn].duration();
	}
	return res;
}

DvString  DvObject::asStr(size_t posn, DvString safeDefault){
	DvString res(safeDefault);
	if(Ddata.size() > posn ) {
		res.set(Ddata[posn]);
	}
	else if(Idata.size() > posn ) {
		res.set(Idata[posn]);
	}
	else if(Sdata.size() > posn ) {
		return Sdata[posn];
	}
	else if(Tdata.size() > posn ) {
		res = Tdata[posn].getISOstring().c_str();
	}
	else if(Edata.size() > posn ) {
		res = Edata[posn].getISOstring().c_str();
	}
	return res;
}

DvString  DvObject::asText(size_t posn, DvString safeDefault){
	// same as asStr() but strings are quoted
	DvString res(safeDefault);
	if(Ddata.size() > posn ) {
		res.set(Ddata[posn]);
	}
	else if(Idata.size() > posn ) {
		res.set(Idata[posn]);
	}
	else if(Sdata.size() > posn ) {
		res="\"";
		res += Sdata[posn];
		res += "\"";
	}
	else if(Tdata.size() > posn ) {
		res = Tdata[posn].getISOstring().c_str();
	}
	else if(Edata.size() > posn ) {
		res = Edata[posn].getISOstring().c_str();
	}
	return res;
}


DvTime  DvObject::asTime(size_t posn, DvTime safeDefault){
	DvTime res(safeDefault);
	
	if(Ddata.size() > posn ) {
		res.setFromEpoch2000(Ddata[posn]);
	}
	if(Idata.size() > posn ) {
		res.setFromEpoch2000((double)Idata[posn]);
	}
	else if(Sdata.size() > posn ) {
		res.setFromISOstring(Sdata[posn]);
	}
	else if(Tdata.size() > posn ) {
		return Tdata[posn];
	}
	else if(Edata.size() > posn ) {
		res = Edata[posn].start();
	}
	return res;
}



DvString DvObject::getTypeAsText(){
	if(is_str()) return DvString("CHAR");
	if(is_dbl()) return DvString("DOUBLE");
	if(is_int()) return DvString("INT");
	if(is_time()) return DvString("ISO_TIME");
	if(is_event()) return DvString("ISO_TIME_RANGE");
	
	return DvString("CHAR"); // ensure something is returned
}

// equivalent object as a double

DvObject_var DvObject::convertToDbl()
{
	DvObject_var res = new DvObject(0.0, this->Dims(), this->seqSize());
	for(size_t i=0; i<res->totalSize(); i++) res->dbl(i) = this->asDouble(i);
	
	// xrefs
	DvNode *xref = this->first_xref();
	while(xref){	
	
  		// take xref for these records	
		if( (is_int() && xref->obj()->is_int()) ||
		    (is_time() && xref->obj()->is_str()) ||
			(is_time() && xref->obj()->is_time()) ||
		    (is_time() && xref->obj()->is_event()) ){
				DvObject_var newXref = xref->obj()->convertToDbl();
				res->change_xref(xref->name(), newXref);
		}
		else res->change_xref(xref->name(), xref->obj());
	
		xref = xref->next;	
	}
	return res;
}


void DvObject::ensureSIC(){
	DvObject_var SIC = get_xref(SI_CONVERSION);
	if(SIC.is_nil()) return; // can do nothing
	
	// Now do xrefs
	DvNode * next = this->first_xref();
	while(next){
		if(next->name() == DELTA_PLUS) {
			if(next->obj()->xref_exists(SI_CONVERSION) == false) next->obj()->change_xref(SI_CONVERSION, SIC);
		}
		if(next->name() == DELTA_MINUS) {
			if(next->obj()->xref_exists(SI_CONVERSION) == false) next->obj()->change_xref(SI_CONVERSION, SIC);
		}
		if(next->name() == SCALEMIN) {
			if(next->obj()->xref_exists(SI_CONVERSION) == false) next->obj()->change_xref(SI_CONVERSION, SIC);
		}
		if(next->name() == SCALEMAX) {
			if(next->obj()->xref_exists(SI_CONVERSION) == false) next->obj()->change_xref(SI_CONVERSION, SIC);
		}
		if(next->name() == VALIDMIN) {
			if(next->obj()->xref_exists(SI_CONVERSION) == false) next->obj()->change_xref(SI_CONVERSION, SIC);
		}
		if(next->name() == VALIDMAX) {
			if(next->obj()->xref_exists(SI_CONVERSION) == false) next->obj()->change_xref(SI_CONVERSION, SIC);
		}
		next->obj()->ensureSIC();
		next=next->next;
	}
	

}


DvObject_var DvObject::changeUnitsTo(DvString SIC, DvString Units)
{
    if( !this->sameBaseSI(SIC.c_str()) ) {
        DvObject_var nilObj = new DvObject();
        nilObj->error("base SI units incompatible");
        return nilObj;
    }
    
	DvString thisConvStr = getBaseSI().before('>');
	DvString argConvStr = SIC.before('>');
	if(thisConvStr.empty() || argConvStr.empty() ) {
		return new DvObject(this);
	}
	
	// must be double as unit change may make int fractional
	DvObject_var res = this->convertToDbl();
	
	double argConv = thisConvStr.toDouble() / argConvStr.toDouble();	
	
	// Convert
  	*res *= argConv; // operates on the valarray

  	res->copy_xrefs_from(*this);
	
	DvObject_var newSIC = new DvObject(SIC);
	res->change_xref(SI_CONVERSION, newSIC);
  
  // Now convert the parameters and xrefs that have to be converted

	DvObject_var xref = res->get_xref(DELTA_PLUS);
	if(xref->is_ok()){
		xref = xref->changeUnitsTo(SIC, Units);
		res->change_xref(DELTA_PLUS, xref);
 	}
	xref = res->get_xref(DELTA_MINUS);
	if(xref->is_ok()){
		xref = xref->changeUnitsTo(SIC, Units);
		res->change_xref(DELTA_MINUS, xref);
 	} 
	xref = res->get_xref(SCALEMIN);
	if(xref->is_ok()){
		xref = xref->changeUnitsTo(SIC, Units);
		res->change_xref(SCALEMIN, xref);
 	} 
	xref = res->get_xref(SCALEMAX);
	if(xref->is_ok()){
		xref = xref->changeUnitsTo(SIC, Units);
		res->change_xref(SCALEMAX, xref);
 	} 
	xref = res->get_xref(VALIDMIN);
	if(xref->is_ok()){
		xref = xref->changeUnitsTo(SIC, Units);
		res->change_xref(VALIDMIN, xref);
 	} 
	xref = res->get_xref(VALIDMAX);
	if(xref->is_ok()){
		xref = xref->changeUnitsTo(SIC, Units);
		res->change_xref(VALIDMAX, xref);
 	} 
	  
  
	if( !Units.empty() )
	{
		xref = new DvObject(Units);
		res->change_xref(UNITS, xref);
	}
	else{
		// create Units from SI_conversion
		// not pretty, but safe
		DvString newUnits = SIC.after(">");
		if(atof(SIC.c_str()) != 1){
			newUnits += " / ";
			newUnits += argConvStr;
		}
		res->change_xref(UNITS, newUnits);

	}
  
	return res;
}

void DvObject::apply_mask(DvMask &msk){
	// new sequence length after masking
	size_t newLen = msk.resultSize();
	size_t oldLen = this->seqSize();
	size_t arrLen = this->arraySize();
	seqLen = newLen;

	if(msk.size() != oldLen) {
		error("[Mask] wrong length");
	 	// mask wrong length
		return;
	}
	if(is_dbl()){
		// copy data out of the way
		valarray <double> temp(0., this->Ddata.size() );
		temp = Ddata;
	
		// fix valarray size
		Ddata.resize(arrLen*newLen, 0.);
		size_t j=0;
		for(size_t i=0; i<oldLen; i++){
			if(msk[i]) {
				// use slices to access old and new records
				slice nsl(j*arrLen, arrLen, 1); 
				slice osl(i*arrLen, arrLen, 1); 
				Ddata[nsl] = temp[osl];
				j++;
			}
		}
	}
	else if(is_int()){
		// copy data out of the way
		valarray <int> temp(0, this->Idata.size() );
		temp = Idata;
		// fix valarray size
		Idata.resize(arrLen*newLen, 0);
		size_t j=0;
		for(size_t i=0; i<oldLen; i++){
			if(msk[i]) {
 				// use slices to access old and new records
				slice nsl(j*arrLen, arrLen, 1); 
				slice osl(i*arrLen, arrLen, 1); 
				Idata[nsl] = temp[osl];
				j++;
			}
		}
	}
	else if(is_str()){
		// copy data out of the way
		valarray <DvString> temp(DvString(), this->Sdata.size() );
		temp = Sdata;
	
		// fix valarray size
		Sdata.resize(arrLen*newLen, DvString());
		size_t j=0;
		for(size_t i=0; i<oldLen; i++){
			if(msk[i]) {
				// use slices to access old and new records
				slice nsl(j*arrLen, arrLen, 1); 
				slice osl(i*arrLen, arrLen, 1); 
				Sdata[nsl] = temp[osl];
				j++;
			}
		}
	}
	else if(is_time()){
		// copy data out of the way
		valarray <DvTime> temp(DvTime(), this->Tdata.size() );
		temp = Tdata;
	
		// fix valarray size
		Tdata.resize(arrLen*newLen, DvTime());
		size_t j=0;
		for(size_t i=0; i<oldLen; i++){
			if(msk[i]) {
				// use slices to access old and new records
				slice nsl(j*arrLen, arrLen, 1); 
				slice osl(i*arrLen, arrLen, 1); 
				Tdata[nsl] = temp[osl];
				j++;
			}
		}
	}
	else if(is_event()){
		// copy data out of the way
		valarray <DvEvent> temp(DvEvent(), this->Edata.size() );
		temp = Edata;
	
		// fix valarray size
		Edata.resize(arrLen*newLen, DvEvent());
		size_t j=0;
		for(size_t i=0; i<oldLen; i++){
			if(msk[i]) {
				// use slices to access old and new records
				slice nsl(j*arrLen, arrLen, 1); 
				slice osl(i*arrLen, arrLen, 1); 
				Edata[nsl] = temp[osl];
				j++;
			}
		}
	}
	
    if(newLen == 0) this->error("Empty object after mask");
        
	// Now do xrefs
	DvNode * next = this->first_xref();
	while(next){
		if(next->obj()->seqSize() == oldLen) {
			// take a copy in case this object is used elsewhere
			DvObject_var xref = new DvObject(next->obj());
			xref->apply_mask(msk);
			next->setObj(xref);
		}
		next=next->next;
	}
	
	
}



DvObject_var DvObject::sqrt(){

	if(is_dbl()){
		DvObject_var res = new DvObject(this); // constructor copies xrefs
		res->sqrtThis();
		return res;
	}
	else return new DvObject();
}

void DvObject::sqrtThis(){	
	
	if(is_dbl()){
		for(size_t i=0; i<totalSize(); i++) Ddata[i] = std::sqrt(Ddata[i]);
		
	}
		
	// Now do xrefs
	
	DvString SIC =	getSICPower(0.5);
	if( !SIC.empty()) this->change_xref(SI_CONVERSION, SIC);

	DvString units =	getUnitsPower(0.5);
	if( !units.empty()) this->change_xref(UNITS, units);
	
	DvString labl_s = "sqrt(";
	labl_s += this->getXrefText("LABLAXIS");
	labl_s += ")";
	this->change_xref( "LABLAXIS", labl_s);

	this->delete_xref(DELTA_PLUS);
	this->delete_xref(DELTA_MINUS);
}

DvObject_var DvObject::abs(){

	DvObject_var res = new DvObject(this); // constructor copies xrefs
	res->absThis();
	return res;
	
}

void DvObject::absThis(){	
	
	if( isThreeVector() && is_dbl()){
		size_t len = seqSize();
		valarray <double> res(len);
		for(size_t i=0; i<len; i++){
			double x = Ddata[i*3];
			double y = Ddata[i*3+1];
			double z = Ddata[i*3+2];
			res[i] = std::sqrt(x*x + y*y + z*z);
		}
		Ddata.resize(len);
		Ddata = res;
		dims.clear();
		this->setFrameAttr("scalar>na");
	}
	else if(is_dbl()){
		for(size_t i=0; i<totalSize(); i++) Ddata[i] = std::fabs(Ddata[i]);
		
	}
	else if(is_int()){
		for(size_t i=0; i<totalSize(); i++) Idata[i] = std::abs(Idata[i]);
		
	}
	
	// LABLAXIS
	DvString lab = this->getXrefText( "LABLAXIS" );
	if(lab.empty()) lab = this->getXrefText( "ObjectName" );
	
	DvString labl_s = "|";
	labl_s += lab;
	labl_s += "|";
	this->change_xref( "LABLAXIS", labl_s);

}

DvObject_var DvObject::log(){

	if(is_dbl()){
		DvObject_var res = new DvObject(this); // constructor copies xrefs
		if( res->logThis() ) return res;
	}
	return new DvObject();
}

bool DvObject::logThis(){	
	
	if( hasUnits() ){
		error("[log] object has units");
		return false;
	}

	if(is_dbl()){
		for(size_t i=0; i<totalSize(); i++) Ddata[i] = std::log(Ddata[i]);
	}
		
	// Now modify xrefs
	
	// LABLAXIS
	DvString lab = this->getXrefText( "LABLAXIS" );
	if(lab.empty()) lab = this->getXrefText( "ObjectName" );
	
	DvString labl_s = "ln(";
	labl_s += lab;
	labl_s += ")";
	this->change_xref( "LABLAXIS", labl_s);

	this->delete_xref(DELTA_PLUS);
	this->delete_xref(DELTA_MINUS);

	return true;
	
}

DvObject_var DvObject::log10(){

	if(is_dbl()){
		DvObject_var res = new DvObject(this); // constructor copies xrefs
		if( res->log10This() ) return res;
	}
	return new DvObject();
}

bool DvObject::log10This(){	
	
	if( hasUnits() ){
		error("[log] object has units");
		return false;
	}
	
	if(is_dbl()){
		for(size_t i=0; i<totalSize(); i++) Ddata[i] = std::log10(Ddata[i]);
	}
		
	// Now modify xrefs
	
	// LABLAXIS
	DvString lab = this->getXrefText( "LABLAXIS" );
	if(lab.empty()) lab = this->getXrefText( "ObjectName" );
	
	DvString labl_s = "log(";
	labl_s += lab;
	labl_s += ")";
	this->change_xref( "LABLAXIS", labl_s);

	this->delete_xref(DELTA_PLUS);
	this->delete_xref(DELTA_MINUS);

	return true;
	

}

DvObject_var DvObject::inverse(){

	DvObject_var res;
	
	if(is_dbl() && arraySize() == 1 ){
		
		res = new DvObject(this); // copies xrefs
		for(size_t i=0; i<totalSize(); i++) res->dbl(i) = 1.0 / Ddata[i];		

	}
	else if(is_dbl() && nDims() == 2 && dims[0] == dims[1]){
		res = matrixInverse(); // constructor copies xrefs	

	}
	else {
		res = new DvObject();

		res->error("Inverse not yet supported on this object type");
		return res;
	}
	
	// Now do xrefs
	DvString SIC =	getSICInverse();
	if( !SIC.empty()) res->change_xref(SI_CONVERSION, SIC);

	DvString units =	getUnitsInverse();
	if( !units.empty()) res->change_xref(UNITS, units);

	// LABLAXIS
	DvString lab = this->getXrefText( "LABLAXIS" );
	if(lab.empty()) lab = this->getXrefText( "ObjectName" );
	
	DvString labl_s = "1 / ";
	labl_s += lab;
	res->change_xref( LABLAXIS, labl_s );

	this->delete_xref(DELTA_PLUS);
	this->delete_xref(DELTA_MINUS);
			
	return res;
}



DvObject_var DvObject::chgSign(){

	if(is_dbl() || is_int()){
		DvObject_var res = new DvObject(this); // constructor copies xrefs
		res->chgSignThis();
		return res;
	}
	else return new DvObject();
}

void DvObject::chgSignThis(){	
	
	if(is_dbl()){
		for(size_t i=0; i<totalSize(); i++) Ddata[i] = -1.0 * Ddata[i];
		
	}
	if(is_int()){
		for(size_t i=0; i<totalSize(); i++) Idata[i] = -1 * Idata[i];
		
	}
}

DvObject_var DvObject::power(double p){

	if(isThreeVector()){
		DvObject_var thisSquare = new DvObject(this);
		thisSquare = this->dot(thisSquare);
		double phalf = 0.5 * p;
		
		return thisSquare->power(phalf);
	}

	if(is_dbl() ){
		DvObject_var res = new DvObject(this); // constructor copies xrefs
		res->powerThis(p);
		return res;
	}
	else return new DvObject();
}

void DvObject::powerThis(double p){	
	
	if(is_dbl()){
		for(size_t i=0; i<totalSize(); i++) Ddata[i] = std::pow(Ddata[i], p);
		
	}
	
	// Now do xrefs
	DvString SIC = getSICPower(p);
	if( !SIC.empty()) this->change_xref(SI_CONVERSION, SIC);

	DvString units = getUnitsPower(p);
	if( !units.empty()) this->change_xref(UNITS, units);

	this->setFrameAttr("scalar>na");
	
	// LABLAXIS
	DvString lab = this->getXrefText( "LABLAXIS" );
	if(lab.empty()) lab = this->getXrefText( "ObjectName" );
	
	DvString labl_s = lab;
	labl_s += " ^ ";
	labl_s += p;
	this->change_xref( "LABLAXIS", labl_s );

	this->delete_xref(DELTA_PLUS);
	this->delete_xref(DELTA_MINUS);

}

DvObject_var DvObject::power(int p){

	if(isThreeVector()){
		DvObject_var thisSquare = new DvObject(this);
		thisSquare = this->dot(thisSquare);
		double phalf = 0.5 * p;
		
		return thisSquare->power(phalf);
	}

	if(is_dbl() || is_int()){
		DvObject_var res = new DvObject(this); // constructor copies xrefs
		res->powerThis(p);
		return res;
	}
	else return new DvObject();
}

void DvObject::powerThis(int p){	
		
	if(is_dbl()){
		for(size_t i=0; i<totalSize(); i++) Ddata[i] = std::pow(Ddata[i], p);
		
	}
	else if(is_int()){
		for(size_t i=0; i<totalSize(); i++) {
			int val = Idata[i];
			for(int j=1; j<p; j++) val *= Idata[i];
			Idata[i] = val;
		}
	}
		// Now do xrefs
	DvString SIC = getSICPower((double)p);
	if( !SIC.empty()) this->change_xref(SI_CONVERSION, SIC);

	DvString units = getUnitsPower((double)p);
	if( !units.empty()) this->change_xref(UNITS, units);

	this->setFrameAttr("scalar>na");
	
	// LABLAXIS
	DvString lab = this->getXrefText( "LABLAXIS" );
	if(lab.empty()) lab = this->getXrefText( "ObjectName" );
	
	DvString labl_s = lab;
	labl_s += " ^ ";
	labl_s += p;
	this->change_xref( "LABLAXIS", labl_s );

}

DvObject_var DvObject::remainder(double d){

	if(is_dbl() ){
		DvObject_var res = new DvObject(this); // constructor copies xrefs
		res->remainderThis(d);
		return res;
	}
	else return new DvObject();
}

void DvObject::remainderThis(double d){	
	
	if(is_dbl()){
		for(size_t i=0; i<totalSize(); i++) {
			if(d < 1.e-30) Ddata[i] = 0.; // trap divide by zero
			else Ddata[i] = std::fmod(Ddata[i], d);
		}
	}
}

DvObject_var DvObject::remainder(DvObject_var &obj){

	if(is_dbl() && obj->is_dbl() && (seqSize() == obj->seqSize() || obj->seqSize() == 1) ){
		DvObject_var res = new DvObject(this); // constructor copies xrefs
		if( res->remainderThis(obj) ) return res;
	}
	return new DvObject();
}

bool DvObject::remainderThis(DvObject_var &obj){	
	
	if(is_dbl() && obj->is_dbl()){
		if( seqSize() == obj->seqSize()){
			if(arraySize() == obj->arraySize()){
				for(size_t i=0; i<totalSize(); i++) {
					if(obj->dbl(i) < 1.e-30 ) Ddata[i] = 0.; // trap divide by zero
					else Ddata[i] = std::fmod(Ddata[i], obj->dbl(i));
				}
			}
		}
		else if( obj->seqSize() == 1){
			if(arraySize() == obj->arraySize()){
				for(size_t i=0; i<seqSize(); i++) {
					for(size_t j=0; j<arraySize(); j++){
						if(obj->dbl(j) < 1.e-30 ) this->dbl(i, j) = 0.; // trap divide by zero
						else this->dbl(i, j) = std::fmod(this->dbl(i, j), obj->dbl(j));
					}
				}
			}
		}
		else {
			this->error("[remainder] size mismatch");	
			return false;
		}
	}
	else {
		this->error("[remainder] only handled for double data");	
		return false;
	}
	return true;
}

DvObject_var DvObject::exp(){

	if(is_dbl() ){
		DvObject_var res = new DvObject(this); // constructor copies xrefs
		if( res->expThis() ) return res;
	}
	return new DvObject();
}

bool DvObject::expThis(){	
	
	if( hasUnits() ){
		error("[exp] object has units");
		return false;
	}
	
	if(is_dbl()){
		for(size_t i=0; i<totalSize(); i++) Ddata[i] = std::exp(Ddata[i]);
	}
	
	// LABLAXIS
	DvString lab = this->getXrefText( "LABLAXIS" );
	if(lab.empty()) lab = this->getXrefText( "ObjectName" );
	
	DvString labl_s("exp(");
	labl_s += lab;
	labl_s += ")";
	this->change_xref( "LABLAXIS", labl_s);

 	// Frame
	this->setFrameAttr("scalar>na");

	return true;
}

DvObject_var DvObject::dot(DvObject_var &obj){

	if( !sameFrame(*obj) ) {
		DvObject_var res = new DvObject();
		res->error("[dot product] frames differ");
		return res;
	}

	if(this->seqSize() == 1 && obj->seqSize() > 1){
		// make 1st object into sequence
		DvObject_var thisSeq = new DvObject(this, obj->seqSize());
		DvObject_var tt = obj->getDep0();
		if(tt.is_ok()) thisSeq->change_xref(DEPEND_0, tt);
		DvObject_var res = thisSeq->dot(obj);
		
		return res;
	}

	if( is_dbl() && obj->is_dbl() && isThreeVector() && obj->isThreeVector() ){
	
		size_t len = seqSize();
		valarray <double> data(len);
		
		if(obj->seqSize() == 1){
			double x = obj->Ddata[0];
			double y = obj->Ddata[1];
			double z = obj->Ddata[2];
			for(size_t i=0; i<len; i++){
				data[i] = Ddata[i*3]*x + Ddata[i*3+1]*y + Ddata[i*3+2]*z;
			}
		}
		else{
	
			DvObject_var arg = obj->Join(*this);  // may be unchanged
	
			for(size_t i=0; i<len; i++){
				data[i] = Ddata[i*3]*arg->Ddata[i*3] + Ddata[i*3+1]*arg->Ddata[i*3+1] + Ddata[i*3+2]*arg->Ddata[i*3+2];
			}			
		}
		
		DvObject_var res = new DvObject(data);
		res->copy_xrefs_from(*this);
		res->setFrameAttr("scalar>na");

		// convert units 	
		DvString SIC = getSICProduct(*obj);
		if(!SIC.empty()) res->change_xref(SI_CONVERSION, SIC);

		DvString unit = getUnitsProduct(*obj);
		if(!unit.empty()) res->change_xref(UNITS, unit);

		return res;
	}
	else  {
		DvObject_var res = new DvObject();
		res->error("[dot product] not three vector");
		return res;
	}
}

DvObject_var DvObject::vec(DvObject_var &obj){
    
    if( !sameFrame(*obj) ) {
        DvObject_var res = new DvObject();
        res->error("[cross product] frames differ");
        return res;
    }
    if(this->seqSize() == 1 && obj->seqSize() > 1){
        // make 1st object into sequence
        DvObject_var thisSeq = new DvObject(this, obj->seqSize());
        DvObject_var tt = obj->getDep0();
        if(tt.is_ok()) thisSeq->change_xref(DEPEND_0, tt);
        DvObject_var res = thisSeq->vec(obj);
        
        return res;
    }
    
    if( is_dbl() && obj->is_dbl() && isThreeVector() && obj->isThreeVector() ){
        size_t len = seqSize();
        DvObject_var res = new DvObject(*this); // sets up dimensions correctly and copies xrefs
        
        if(obj->seqSize() == 1){
            
            double x = obj->Ddata[0];
            double y = obj->Ddata[1];
            double z = obj->Ddata[2];
            for(size_t i=0; i<len; i++){
                res->Ddata[i*3]   = Ddata[i*3+1]*z - Ddata[i*3+2]*y;
                res->Ddata[i*3+1] = Ddata[i*3+2]*x - Ddata[i*3]*z;
                res->Ddata[i*3+2] = Ddata[i*3]*y   - Ddata[i*3+1]*x;
            }
        }
        else{
            
            DvObject_var arg = obj->Join(*this); // may be unchanged
            
            for(size_t i=0; i<len; i++){
                res->Ddata[i*3]   = Ddata[i*3+1]*arg->Ddata[i*3+2] - Ddata[i*3+2]*arg->Ddata[i*3+1];
                res->Ddata[i*3+1] = Ddata[i*3+2]*arg->Ddata[i*3]   - Ddata[i*3]*arg->Ddata[i*3+2];
                res->Ddata[i*3+2] = Ddata[i*3]*arg->Ddata[i*3+1]   - Ddata[i*3+1]*arg->Ddata[i*3];
            }
        }
        
        // convert units 	
        DvString SIC = getSICProduct(*obj);
        if(!SIC.empty()) res->change_xref(SI_CONVERSION, SIC);
        
        DvString unit = getUnitsProduct(*obj);
        if(!unit.empty()) res->change_xref(UNITS, unit);
        
        // catch dummy frame
        if(res->getFrame() == "any"){
            DvString otherFrame = obj->getFrameAttr();
            res->setFrameAttr(otherFrame);
        }
        
        return res;
    }
    else  {
        DvObject_var res = new DvObject();
        res->error("[cross product] not three vector");
        return res;
    }
}

DvObject_var DvObject::normalize(){
    
    DvObject_var res;
    if( is_dbl() && isThreeVector() ){
        size_t len = seqSize();
        res = *this / *(this->abs());
    }
    else  {
        res = new DvObject();
        res->error("[Normalize] not three vector");
    }
    
    return res;
}

DvObject_var DvObject::multElements(DvObject_var &obj){
	
	DvObject_var res = new DvObject();

	size_t arrLen = arraySize(); 
	if ( arrLen != obj->arraySize() )  {
		res->error("[multiply elements] array dimensions do not match");
		return res;
	}

	if(seqSize() == 1 && obj->seqSize() != 1) {
		// do other way round to get sequence
		DvObject_var thisObj = this; // get var ptr
		return obj->multElements(thisObj);
	}
		
	if( is_dbl() && obj->is_dbl()  ){
		size_t len = seqSize();
		res = new DvObject(*this); // sets up dimensions correctly and copies xrefs

		if(obj->seqSize() == 1){
			for(size_t i=0; i<len; i++){
				for(size_t j=0; j<arrLen; j++){
					res->dbl(i, j) *= obj->dbl(j);
				}
			}
		}
		else{
			DvObject_var arg = obj->Join(*this); // may be unchanged
		
			for(size_t i=0; i<len; i++){
				for(size_t j=0; j<arrLen; j++){
					res->dbl(i, j) *= arg->dbl(i, j);
				}
			}
		}
	}
	else if( is_dbl() && obj->is_int()  ){
		size_t len = seqSize();
		res = new DvObject(*this); // sets up dimensions correctly and copies xrefs

		if(obj->seqSize() == 1){
			for(size_t i=0; i<len; i++){
				for(size_t j=0; j<arrLen; j++){
					res->dbl(i, j) *= obj->itg(j);
				}
			}
		}
		else{
	
			DvObject_var arg = obj->Join(*this); // may be unchanged
		
			for(size_t i=0; i<len; i++){
				for(size_t j=0; j<arrLen; j++){
					res->dbl(i, j) *= arg->itg(i, j);
				}
			}
		}
	}
	else if( is_int() && obj->is_int()  ){
		size_t len = seqSize();
		DvObject *res = new DvObject(*this); // sets up dimensions correctly and copies xrefs

		if(obj->seqSize() == 1){
			for(size_t i=0; i<len; i++){
				for(size_t j=0; j<arrLen; j++){
					res->itg(i, j) *= obj->itg(j);
				}
			}
		}
		else{
	
			DvObject_var arg = obj->Join(*this); // may be unchanged
		
			for(size_t i=0; i<len; i++){
				for(size_t j=0; j<arrLen; j++){
					res->itg(i, j) *= arg->itg(i, j);
				}
			}
		}
	}
	else  {
		res->error("[piecewise product] not possible");
		return res;
	}
	
	// convert units 	
	DvString SIC = getSICProduct(*obj);
	if(!SIC.empty()) res->change_xref(SI_CONVERSION, SIC);

	DvString unit = getUnitsProduct(*obj);
	if(!unit.empty()) res->change_xref(UNITS, unit);
		
	return res;

}
DvObject_var DvObject::divideElements(DvObject_var &obj){

	DvObject_var res = new DvObject();
	
	size_t arrLen = arraySize(); 
	if ( arrLen != obj->arraySize() )  {
		res->error("[multiply elements] array dimensions do not match");
		return res;
	}

	if(seqSize() == 1 && obj->seqSize() != 1) {
		// do other way round to get sequence
		DvObject_var thisObj = this; // get var ptr
		DvObject_var invObj = obj->inverse();
		return invObj->multElements(thisObj); // now multiply
	}
		
	if( is_dbl() && obj->is_dbl()  ){
		size_t len = seqSize();
		res = new DvObject(*this); // sets up dimensions correctly and copies xrefs

		if(obj->seqSize() == 1){
			for(size_t i=0; i<len; i++){
				for(size_t j=0; j<arrLen; j++){
					res->dbl(i, j) /= obj->dbl(j);
				}
			}
		}
		else{
	
			DvObject_var jObj = obj->Join(*this); // may be unchanged
		
			for(size_t i=0; i<len; i++){
				for(size_t j=0; j<arrLen; j++){
					res->dbl(i, j) /= jObj->dbl(i, j);
				}
			}
		}
		
		// convert units 	
		DvString SIC = getSICProduct(*obj);
		if(!SIC.empty()) res->change_xref(SI_CONVERSION, SIC);

		DvString unit = getUnitsProduct(*obj);
		if(!unit.empty()) res->change_xref(UNITS, unit);
		
	}
	else if( is_dbl() && obj->is_int()  ){
		size_t len = seqSize();
		
		res = new DvObject(*this); // sets up dimensions correctly and copies xrefs

		if(obj->seqSize() == 1){
			for(size_t i=0; i<len; i++){
				for(size_t j=0; j<arrLen; j++){
					res->dbl(i, j) /= obj->itg(j);
				}
			}
		}
		else{
	
			DvObject_var jObj = obj->Join(*this); // may be unchanged
		
			for(size_t i=0; i<len; i++){
				for(size_t j=0; j<arrLen; j++){
					res->dbl(i, j) /= jObj->itg(i, j);
				}
			}
		}
		
		// convert units 	
		DvString SIC = getSICProduct(*obj);
		if(!SIC.empty()) res->change_xref(SI_CONVERSION, SIC);

		DvString unit = getUnitsProduct(*obj);
		if(!unit.empty()) res->change_xref(UNITS, unit);
		
	}
	else if( is_int() && obj->is_int()  ){
		size_t len = seqSize();
		res = new DvObject(*this); // sets up dimensions correctly and copies xrefs

		if(obj->seqSize() == 1){
			for(size_t i=0; i<len; i++){
				for(size_t j=0; j<arrLen; j++){
					res->itg(i, j) /= obj->itg(j);
				}
			}
		}
		else{
	
			DvObject_var jObj = obj->Join(*this); // may be unchanged
		
			for(size_t i=0; i<len; i++){
				for(size_t j=0; j<arrLen; j++){
					res->itg(i, j) /= jObj->itg(i, j);
				}
			}
		}
		
	}
	else  {
		res->error("[piecewise product] not possible");
		return res;
	}
	
	// convert units 
	
	DvString SIC = getSICRatio(*obj);
	if(!SIC.empty()) res->change_xref(SI_CONVERSION, SIC);

	DvString unit = getUnitsRatio(*obj);
	if(!unit.empty()) res->change_xref(UNITS, unit);

	// LABLAXIS
	DvString labl_s = this->getXrefText( "LABLAXIS" );
	if(obj->xref_exists("LABLAXIS")){
		labl_s += "/";
		labl_s += obj->getXrefText( "LABLAXIS" );
	}
	else if(obj->seqSize() == 1){
		
		labl_s += " / ";
		labl_s += obj->asStr(0);
	}
	res->change_xref( "LABLAXIS", labl_s);

	return res;

}


double DvObject::max(){

	if(is_dbl() ){
		double val = Ddata.max();
		if( isnan(val) ){
		
			// handle case with NaN
			val = -1.e30;
			for(size_t i=0; i<Ddata.size(); i++) if( !isnan(Ddata[i]) && Ddata[i] > val) val = Ddata[i];
		}
		return val;
	}
	if(is_int() ){
		return (double) Idata.max();
	}
	if(is_time() ){
		double max = Tdata[Tdata.size()-1].getEpoch2000();
		for(int i=Tdata.size()-2; i>=0; i--){
			if(Tdata[i].getEpoch2000() > max) max = Tdata[i].getEpoch2000();
		}
		return max;
	}
	else return 1.e30;
}

double DvObject::min(){

	if(is_dbl() ){
		double val = Ddata.min();
		if( isnan(val) ){
		
			// handle case with NaN
			val = 1.e30;
			for(size_t i=0; i<Ddata.size(); i++) if( !isnan(Ddata[i]) && Ddata[i] < val) val = Ddata[i];
		}
		return val;
	}
	if(is_int() ){
		return (double) Idata.min();
	}
	if(is_time() ){
		double min = Tdata[0].getEpoch2000();
		for(size_t i=1; i<Tdata.size(); i++){
			if(Tdata[i].getEpoch2000() < min) min = Tdata[i].getEpoch2000();
		}
		return min;
	}
	else return -1.e30;
}

// trig functions
// -------------------------------------


DvObject_var DvObject::toRad(){

	if(is_dbl() ){
		DvObject_var res = new DvObject(this); // constructor copies xrefs
		res->toRadThis();
		return res;
	}
	else return new DvObject();
}

void DvObject::toRadThis(){	
	
	double thisConv = DV_PI/180.;

	if(is_dbl()){
		if( isRad() ) return;
		else if(isDeg()){
			for(size_t i=0; i<totalSize(); i++) Ddata[i] = Ddata[i]*thisConv;
		}
		else if( sameBaseSI( "1>rad" ) ){
			thisConv = getBaseSI().before('>').toDouble();
			for(size_t i=0; i<totalSize(); i++) Ddata[i] = Ddata[i]*thisConv;
		}
		else if( sameBaseSI("1>deg") ){
			thisConv *= getBaseSI().before('>').toDouble(); // to rad (already has deg2rad)
			for(size_t i=0; i<totalSize(); i++) Ddata[i] = Ddata[i]*thisConv;
		}
		else {
			error("[toRad] unknown units, assuming radian");
			return;
		}
	}
	else{
		error("[toRad] not possible on data type");
		return;
	}
	
	// fix attributes
	DvString SIC("1.0>rad");
	change_xref(SI_CONVERSION, SIC);

	DvString units("rad");
	change_xref(UNITS, units);

	if(xref_exists(SCALEMIN) ){
		DvObject_var minVal = get_xref(SCALEMIN);
		*minVal *= thisConv;
		change_xref(SCALEMIN, minVal);
	}
	
	if(xref_exists(SCALEMAX) ){
		DvObject_var maxVal = get_xref(SCALEMAX);
		*maxVal *= thisConv;
		change_xref(SCALEMAX, maxVal);
	}
	
}

DvObject_var DvObject::toDeg(){

	if(is_dbl() ){
		DvObject_var res = new DvObject(this); // constructor copies xrefs
		res->toDegThis();
		return res;
	}
	else return new DvObject();
}

void DvObject::toDegThis(){	
	
	double thisConv = 180./DV_PI;

	if(is_dbl()){
		if( isDeg() ) return;
		else if(isRad()){
			for(size_t i=0; i<totalSize(); i++) Ddata[i] = Ddata[i]*thisConv;
		}
		else if( sameBaseSI("1>deg") ) {
			thisConv = getBaseSI().before('>').toDouble();
			for(size_t i=0; i<totalSize(); i++) Ddata[i] = Ddata[i]*thisConv;
		}
		else if ( sameBaseSI("1>rad") ){
			thisConv *= getBaseSI().before('>').toDouble(); // to rad (already has rad2deg)
			for(size_t i=0; i<totalSize(); i++) Ddata[i] = Ddata[i]*thisConv;
		}
		else {
			error("[toDeg] unknown units, assuming degree");
			return;
		}
	}
	else{
		error("[toDeg] not possible on data type");
		return;
	}
	// fix attributes
	DvString SIC("1.0>deg");
	change_xref(SI_CONVERSION, SIC);

	DvString units("deg");
	change_xref(UNITS, units);

	if(xref_exists(SCALEMIN) ){
		DvObject_var minVal = get_xref(SCALEMIN);
		*minVal *= thisConv;
		change_xref(SCALEMIN, minVal);
	}
	
	if(xref_exists(SCALEMAX) ){
		DvObject_var maxVal = get_xref(SCALEMAX);
		*maxVal *= thisConv;
		change_xref(SCALEMAX, maxVal);
	}
		
}

	// Cartesian to Polar coords (radians)

DvObject_var DvObject::XYtoRP(){

	DvObject_var res = new DvObject(this); // constructor copies xrefs
	res->XYtoRPThis();
	return res;
}

void DvObject::XYtoRPThis(){	
// Transform two-dimensional Cartesian, or xy, coordinates to polar coordinates
// ((r),(p)hi)
// phi is a counterclockwise angular displacement in radians from the
// positive x-axis.
// r is the magnitude of the vector
// rp_out is an indexable class where rp_out[0]=(r), rp_out[1]=(p)hi
// Similarly, xy_in[0]=x, xy_in[1]=y

	if(is_dbl() ){
		double r=0.,p=0.;
		for(size_t n=0; n<seqSize(); n++){
			r = std::sqrt(Ddata[2*n]*Ddata[2*n] + Ddata[2*n+1]*Ddata[2*n+1]);
			p = atan_yx(Ddata[2*n+1], Ddata[2*n]);
			Ddata[2*n] = r;
			Ddata[2*n+1] = p;
		}
	}
	else{
		error("[XYtoRP] not possible on data type");
		return;
	}
	
	// fix attributes
	DvString SIC = getXrefText(SI_CONVERSION);
	vector <size_t> dims;
	dims.push_back(2);
	DvObject_var newSIC = new DvObject(SIC, dims);
	newSIC->str(1) = "1>rad";
	change_xref(SI_CONVERSION, newSIC);

	DvString unit = getXrefText(UNITS);
	DvObject_var newUnits = new DvObject(unit, dims);
	newUnits->str(1) = "rad";
	change_xref(UNITS, newUnits);;

}

DvObject_var DvObject::XYtoRP_deg(){

	DvObject_var res = new DvObject(this); // constructor copies xrefs
	res->XYtoRP_degThis();
	return res;
}

void DvObject::XYtoRP_degThis(){	

	double thisConv = 180./DV_PI;

	if(is_dbl() ){
		double r=0.,p=0.;
		for(size_t n=0; n<seqSize(); n++){
			r = std::sqrt(Ddata[2*n]*Ddata[2*n] + Ddata[2*n+1]*Ddata[2*n+1]);
			p = atan_yx(Ddata[2*n+1], Ddata[2*n]);
			Ddata[2*n] = r;
			Ddata[2*n+1] = p * thisConv;
		}
	}
	else{
		error("[XYtoRP] not possible on data type");
		return;
	}
	
	// fix attributes
	DvString SIC = getXrefText(SI_CONVERSION);
	vector <size_t> dims;
	dims.push_back(2);
	DvObject_var newSIC = new DvObject(SIC, dims);
	newSIC->str(1) = "1>deg";
	change_xref(SI_CONVERSION, newSIC);

	DvString unit = getXrefText(UNITS);
	DvObject_var newUnits = new DvObject(unit, dims);
	newUnits->str(1) = "deg";
	change_xref(UNITS, newUnits);;

}

	// Polar to Cartesian  coords (radians)

DvObject_var DvObject::RPtoXY(){

	DvObject_var res = new DvObject(this); // constructor copies xrefs
	res->RPtoXYThis();
	return res;
}

void DvObject::RPtoXYThis(){	
// Transform  polar coordinates ((r),(p)hi) to two-dimensional Cartesian,
// or xy, coordinates
// phi is a counterclockwise angular displacements in radians from the
// positive x-axis.

	double thisConv = 1.;

	if(is_dbl() ){
		if(isVectDeg()){
			thisConv = DV_PI/180.;
		}
		double r=0.,p=0.;
		for(size_t n=0; n<seqSize(); n++){
			r = Ddata[2*n];
			p = Ddata[2*n+1] * thisConv;
			Ddata[2*n]   = r * std::cos(p); // x
			Ddata[2*n+1] = r * std::sin(p); // y
		}
	}
	else{
		error("[RPtoXY] not possible on data type");
		return;
	}
	
	// fix attributes
	DvObject_var SIC = get_xref(SI_CONVERSION);
	change_xref(SI_CONVERSION, SIC->str()); // uses 1st (r) component

	DvObject_var units = get_xref(UNITS);
	change_xref(UNITS, units->str());

		
}

	// Cylindrical Polar to Cartesian  coords (radians)

DvObject_var DvObject::RPZtoXYZ(){

	DvObject_var res = new DvObject(this); // constructor copies xrefs
	res->RPZtoXYZThis();
	return res;
}

void DvObject::RPZtoXYZThis(){	
// Transform cylindrical coordinates ((R),(p)hi,(z)) to Cartesian, or
// xyz, coordinates
// R is the distance from the origin to a point in the x-y plane, phi is a
// counterclockwise angular displacement in radians from the 
// positive x-axis, and z is the height above the x-y plane.
// Rpz_in is an indexable class where Rpz_in[0]=(R), Rpz_in[1]=(p)hi
// and Rpz_in[2]=(z).
// Similarly, xyz_out[0]=x, xyz_out[1]=y, xyz_out[2]=z

	double thisConv = 1.;

	if(is_dbl() ){
		if(isVectDeg()){
			thisConv = DV_PI/180.;
		}
		double r=0.,p=0.;
		for(size_t n=0; n<seqSize(); n++){
			r = Ddata[3*n];
			p = Ddata[3*n+1] * thisConv;
			Ddata[3*n]   = r * std::cos(p); // x
			Ddata[3*n+1] = r * std::sin(p); // y
			// z unchanged
		}
	}
	else{
		error("[RPtoXY] not possible on data type");
		return;
	}
	
	// fix attributes
	DvObject_var SIC = get_xref(SI_CONVERSION);
	change_xref(SI_CONVERSION, SIC->str()); // uses 1st (r) component

	DvObject_var units = get_xref(UNITS);
	change_xref(UNITS, units->str()); // uses 1st (r) component

		
}

	// Cartesian to Cylindrical Polar coords 

DvObject_var DvObject::XYZtoRPZ(){

	DvObject_var res = new DvObject(this); // constructor copies xrefs
	res->XYZtoRPZThis();
	return res;
}

void DvObject::XYZtoRPZThis(){	
// Transform Cartesian, or xyz, coordinates to cylindrical coordinates
// ((R),(p)hi,(z))
// R is the distance from the origin to a point in the x-y plane, phi is a
// counterclockwise angular displacement in degrees from the 
// positive x-axis, and z is the height above the x-y plane.
// Rpz_out is an indexable class where Rpz_out[0]=(R), Rpz_out[1]=(p)hi
// and Rpz_out[2]=(z).
// Similarly, xyz_in[0]=x, xyz_in[1]=y, xyz_in[2]=z

	if(is_dbl() ){
		double r=0.,p=0.;
		for(size_t n=0; n<seqSize(); n++){
			r = std::sqrt(Ddata[3*n]*Ddata[3*n] + Ddata[3*n+1]*Ddata[3*n+1]);
			p = atan_yx(Ddata[3*n+1], Ddata[3*n]);
			Ddata[3*n] = r;
			Ddata[3*n+1] = p;
			// z unchanged
		}
	}
	else{
		error("[XYZtoRPZ] not possible on data type");
		return;
	}
	
	// fix attributes
	DvString SIC = getXrefText(SI_CONVERSION);
	vector <size_t> dims;
	dims.push_back(3);
	DvObject_var newSIC = new DvObject(SIC, dims);
	newSIC->str(1) = "1>rad";
	change_xref(SI_CONVERSION, newSIC);

	DvString unit = getXrefText(UNITS);
	DvObject_var newUnits = new DvObject(unit, dims);
	newUnits->str(1) = "rad";
	change_xref(UNITS, newUnits);;

}

DvObject_var DvObject::XYZtoRPZ_deg(){

	DvObject_var res = new DvObject(this); // constructor copies xrefs
	res->XYZtoRPZ_degThis();
	return res;
}

void DvObject::XYZtoRPZ_degThis(){	

	double thisConv = 180./DV_PI;

	if(is_dbl() ){
		double r=0.,p=0.;
		for(size_t n=0; n<seqSize(); n++){
			r = std::sqrt(Ddata[3*n]*Ddata[3*n] + Ddata[3*n+1]*Ddata[3*n+1]);
			p = atan_yx(Ddata[3*n+1], Ddata[3*n]);
			Ddata[3*n] = r;
			Ddata[3*n+1] = p * thisConv;
			// z unchanged
		}
	}
	else{
		error("[XYZtoRPZ deg] not possible on data type");
		return;
	}
	
	// fix attributes
	DvString SIC = getXrefText(SI_CONVERSION);
	vector <size_t> dims;
	dims.push_back(3);
	DvObject_var newSIC = new DvObject(SIC, dims);
	newSIC->str(1) = "1>deg";
	change_xref(SI_CONVERSION, newSIC);

	DvString unit = getXrefText(UNITS);
	DvObject_var newUnits = new DvObject(unit, dims);
	newUnits->str(1) = "deg";
	change_xref(UNITS, newUnits);;

}


	//  Polar RTP to Cartesian  coords 

DvObject_var DvObject::RTPtoXYZ(){

	DvObject_var res = new DvObject(this); // constructor copies xrefs
	res->RTPtoXYZThis();
	return res;
}

void DvObject::RTPtoXYZThis(){	
// Transform spherical coordinates ((r),(t)heta,(p)hi) to Cartesian, or
// xyz, coordinates
// Theta and phi are angular displacements in radians from the positive z-axis
// and positive x-axis, respectively
// theta in range [0,pi), phi in range [0,2pi)
// rtp_in is an indexable class where rtp_in[0]=(r), rtp_in[1]=(t)heta
// and rtp_in[2]=(p)hi.
// Similarly, xyz_out[0]=x, xyz_out[1]=y, xyz_out[2]=z

	double thisConv = 1.;

	if(is_dbl() ){
		if(isVectDeg()){
			thisConv = DV_PI/180.;
		}
		double r=0., t=0., p=0., xy=0.;
		for(size_t n=0; n<seqSize(); n++){
			r = Ddata[3*n];
			t = Ddata[3*n+1] * thisConv;
			p = Ddata[3*n+2] * thisConv;
			xy = r * std::sin(t);
			Ddata[3*n]   = xy * std::cos(p); // x
			Ddata[3*n+1] = xy * std::sin(p); // y
			Ddata[3*n+2] = r  * std::cos(t); // z
		}
	}
	else{
		error("[RTPtoXYZ] not possible on data type");
		return;
	}
	
	// fix attributes
	DvObject_var SIC = get_xref(SI_CONVERSION);
	change_xref(SI_CONVERSION, SIC->str()); // uses 1st (r) component

	DvObject_var units = get_xref(UNITS);
	change_xref(UNITS, units->str()); // uses 1st (r) component

		
}

	// Cartesian to Polar RTP coords 

DvObject_var DvObject::XYZtoRTP(){

	DvObject_var res = new DvObject(this); // constructor copies xrefs
	res->XYZtoRTPThis();
	return res;
}

void DvObject::XYZtoRTPThis(){	
// Transform Cartesian, or xyz, coordinates to spherical coordinates
// ((r),(t)heta,(p)hi)
// Theta and phi are angular displacements in degrees from the positive z-axis
// and positive x-axis, respectively
// theta in range [0,180), phi in range [0,360)
// rtp_out is an indexable class where rtp_out[0]=(r), rtp_out[1]=(t)heta
// and rtp_out[2]=(p)hi.
// Similarly, xyz_in[0]=x, xyz_in[1]=y, xyz_in[2]=z

	if(is_dbl() ){
		double r=0., t=0., p=0., xy_sq=0.;
		for(size_t n=0; n<seqSize(); n++){
			xy_sq = Ddata[3*n]*Ddata[3*n] + Ddata[3*n+1]*Ddata[3*n+1];
			r = std::sqrt(xy_sq + Ddata[3*n+2]*Ddata[3*n+2]);
			p = atan_yx(Ddata[3*n+1], Ddata[3*n]);
			t = atan_yx(std::sqrt(xy_sq), Ddata[3*n+2]);
			Ddata[3*n] = r;
			Ddata[3*n+1] = t;
			Ddata[3*n+2] = p;
		}
	}
	else{
		error("[XYZtoRTP] not possible on data type");
		return;
	}
	
	// fix attributes
	DvString SIC = getXrefText(SI_CONVERSION);
	vector <size_t> dims;
	dims.push_back(3);
	DvObject_var newSIC = new DvObject(SIC, dims);
	newSIC->str(1) = "1>rad";
	newSIC->str(2) = "1>rad";
	change_xref(SI_CONVERSION, newSIC);

	DvString unit = getXrefText(UNITS);
	DvObject_var newUnits = new DvObject(unit, dims);
	newUnits->str(1) = "rad";
	newUnits->str(2) = "rad";
	change_xref(UNITS, newUnits);;

}

DvObject_var DvObject::XYZtoRTP_deg(){
    
	DvObject_var res = new DvObject(this); // constructor copies xrefs
	res->XYZtoRTP_degThis();
	return res;
}

void DvObject::XYZtoRTP_degThis(){	

	double thisConv = 180./DV_PI;

	if(is_dbl() ){
		double r=0., t=0., p=0., xy_sq=0.;
		for(size_t n=0; n<seqSize(); n++){
			xy_sq = Ddata[3*n]*Ddata[3*n] + Ddata[3*n+1]*Ddata[3*n+1];
			r = std::sqrt(xy_sq + Ddata[3*n+2]*Ddata[3*n+2]);
			p = atan_yx(Ddata[3*n+1], Ddata[3*n]);
			t = atan_yx(std::sqrt(xy_sq), Ddata[3*+2]);
			Ddata[3*n] = r;
			Ddata[3*n+1] = t * thisConv;
			Ddata[3*n+2] = p * thisConv;
		}
	}
	else{
		error("[XYZtoRTP deg] not possible on data type");
		return;
	}
	
	// fix attributes
	DvString SIC = getXrefText(SI_CONVERSION);
	vector <size_t> dims;
	dims.push_back(3);
	DvObject_var newSIC = new DvObject(SIC, dims);
	newSIC->str(1) = "1>deg";
	newSIC->str(2) = "1>deg";
	change_xref(SI_CONVERSION, newSIC);

	DvString unit = getXrefText(UNITS);
	DvObject_var newUnits = new DvObject(unit, dims);
	newUnits->str(1) = "deg";
	newUnits->str(2) = "deg";
	change_xref(UNITS, newUnits);;

}



	//  Polar RLP to Cartesian  coords 

DvObject_var DvObject::RLPtoXYZ(){

	DvObject_var res = new DvObject(this); // constructor copies xrefs
	res->RLPtoXYZThis();
	return res;
}

void DvObject::RLPtoXYZThis(){	
// Transform coordinates expressed in terms of magnitude, latitude, and
// longitude (r,(l)atitude,(p)hi) to Cartesian, or xyz, coordinates
// latitude = pi/2 - theta, where theta is the angular distance in radians 
// from the positive z axis, phi is the angular distance in radians from the
// positive x axis and is equal to the longitude, and r is the magnitude of 
// the radial distance.
// latitude in range [pi/2,-pi/2), phi in range [0,2pi)
// rlp_in is an indexable class where rlp_in[0]=(r), rlp_in[1]=(l)atitude,
// and rlp_in[2]=(p)hi
// Similarly, xyz_out[0]=x, xyz_out[1]=y, xyz_out[2]=z

	double thisConv = 1.;

	if(is_dbl() ){
		if(isVectDeg()){
			thisConv = DV_PI/180.;
		}
		double r=0., l=0., p=0., xy=0.;
		for(size_t n=0; n<seqSize(); n++){
			r = Ddata[3*n];
			l = Ddata[3*n+1] * thisConv;
			p = Ddata[3*n+2] * thisConv;
			xy = r * std::cos(l);
			Ddata[3*n]   = xy * std::cos(p); // x
			Ddata[3*n+1] = xy * std::sin(p); // y
			Ddata[3*n+2] = r  * std::sin(l); // z
		}
	}
	else{
		error("[RLPtoXYZ] not possible on data type");
		return;
	}
	
	// fix attributes
	DvObject_var SIC = get_xref(SI_CONVERSION);
	change_xref(SI_CONVERSION, SIC->str()); // uses 1st (r) component

	DvObject_var units = get_xref(UNITS);
	change_xref(UNITS, units->str()); // uses 1st (r) component

		
}


	// Cartesian to Polar RLP coords 

DvObject_var DvObject::XYZtoRLP(){

	DvObject_var res = new DvObject(this); // constructor copies xrefs
	res->XYZtoRLPThis();
	return res;
}

void DvObject::XYZtoRLPThis(){	
// Transform Cartesian, or xyz, coordinates to coordinates expressed in terms
// of magnitude, latitude, and longitude (r,(l)atitude,(p)hi)
// latitude = pi/2 - theta, where theta is the angular distance in radians 
// from the positive z axis, phi is the angular distance in radians from the
// positive x axis and is equal to the longitude, and r is the magnitude of 
// the radial distance.
// latitude in range [pi/2,-pi/2), phi in range [0,2pi)

	if(is_dbl() ){
		double r=0., l=0., p=0., xy_sq=0.;
		for(size_t n=0; n<seqSize(); n++){
			xy_sq = Ddata[3*n]*Ddata[3*n] + Ddata[3*n+1]*Ddata[3*n+1];
			r = std::sqrt(xy_sq + Ddata[3*n+2]*Ddata[3*n+2]);
			p = atan_yx(Ddata[3*n+1], Ddata[3*n]);
			l = std::atan2(Ddata[3*n+2], std::sqrt(xy_sq));
			Ddata[3*n] = r;
			Ddata[3*n+1] = l;
			Ddata[3*n+2] = p;
		}
	}
	else{
		error("[XYZtoRLP] not possible on data type");
		return;
	}
	
	// fix attributes
	DvString SIC = getXrefText(SI_CONVERSION);
	vector <size_t> dims;
	dims.push_back(3);
	DvObject_var newSIC = new DvObject(SIC, dims);
	newSIC->str(1) = "1>rad";
	newSIC->str(2) = "1>rad";
	change_xref(SI_CONVERSION, newSIC);

	DvString unit = getXrefText(UNITS);
	DvObject_var newUnits = new DvObject(unit, dims);
	newUnits->str(1) = "rad";
	newUnits->str(2) = "rad";
	change_xref(UNITS, newUnits);;

}


DvObject_var DvObject::XYZtoRLP_deg(){

	DvObject_var res = new DvObject(this); // constructor copies xrefs
	res->XYZtoRLPThis();
	return res;
}

void DvObject::XYZtoRLP_degThis(){	

	double thisConv = 180./DV_PI;

	if(is_dbl() ){
		double r=0., l=0., p=0., xy_sq=0.;
		for(size_t n=0; n<seqSize(); n++){
			xy_sq = Ddata[3*n]*Ddata[3*n] + Ddata[3*n+1]*Ddata[3*n+1];
			r = std::sqrt(xy_sq + Ddata[3*n+2]*Ddata[3*n+2]);
			p = atan_yx(Ddata[3*n+1], Ddata[3*n]);
			l = std::atan2(Ddata[3*n+2], std::sqrt(xy_sq));
			Ddata[3*n] = r;
			Ddata[3*n+1] = l * thisConv;
			Ddata[3*n+2] = p * thisConv;
		}
	}
	else{
		error("[XYZtoRLP deg] not possible on data type");
		return;
	}
	
	// fix attributes
	DvString SIC = getXrefText(SI_CONVERSION);
	vector <size_t> dims;
	dims.push_back(3);
	DvObject_var newSIC = new DvObject(SIC, dims);
	newSIC->str(1) = "1>deg";
	newSIC->str(2) = "1>deg";
	change_xref(SI_CONVERSION, newSIC);

	DvString unit = getXrefText(UNITS);
	DvObject_var newUnits = new DvObject(unit, dims);
	newUnits->str(1) = "deg";
	newUnits->str(2) = "deg";
	change_xref(UNITS, newUnits);;

}



// Trig

DvObject_var DvObject::sin(){

	if(is_dbl() ){
		DvObject_var res = new DvObject(this); // constructor copies xrefs
		res->sinThis();
		return res;
	}
	else return new DvObject();
}
void DvObject::sinThis(){	
	
	if(is_dbl()){
		// convert this to radians
		toRadThis();
		
		for(size_t i=0; i<totalSize(); i++) Ddata[i] = std::sin(Ddata[i]);
	}
	else{
		error("[sin] not possible on data type");
		return;
	}
	
	// fix attributes
	DvString SIC("1.0>(sin)");
	change_xref(SI_CONVERSION, SIC);

	DvString units("");
	change_xref(UNITS, units);

	DvObject_var minVal = new DvObject(-1.0);
	change_xref(SCALEMIN, minVal);

	DvObject_var maxVal = new DvObject(1.0);
	change_xref(SCALEMAX, maxVal);

	DvString lab = this->getXrefText( "LABLAXIS" );
	if(lab.empty()) lab = this->getXrefText( "ObjectName" );
	
	DvString labAxis("sin( ");
	labAxis += lab;
	labAxis += " )";
	change_xref(LABLAXIS, labAxis);

	DvString fieldname("sin( ");
	labAxis += getXrefText(FIELDNAM);
	labAxis += " )";
	change_xref(FIELDNAM, labAxis);
}

//-----
DvObject_var DvObject::cos(){

	if(is_dbl() ){
		DvObject_var res = new DvObject(this); // constructor copies xrefs
		res->cosThis();
		return res;
	}
	else return new DvObject();
}
void DvObject::cosThis(){	
	
	if(is_dbl()){
		// convert this to radians
		toRadThis();
		
		for(size_t i=0; i<totalSize(); i++) Ddata[i] = std::cos(Ddata[i]);
	}
	else{
		error("[cos] not possible on data type");
		return;
	}
	
	// fix attributes
	DvString SIC("1.0>(cos)");
	change_xref(SI_CONVERSION, SIC);

	DvString units("");
	change_xref(UNITS, units);

	DvObject_var minVal = new DvObject(-1.0);
	change_xref(SCALEMIN, minVal);

	DvObject_var maxVal = new DvObject(1.0);
	change_xref(SCALEMAX, maxVal);

	DvString lab = this->getXrefText( "LABLAXIS" );
	if(lab.empty()) lab = this->getXrefText( "ObjectName" );
	
	DvString labAxis("cos( ");
	labAxis += lab;
	labAxis += " )";
	change_xref(LABLAXIS, labAxis);

	DvString fieldname("cos( ");
	labAxis += getXrefText(FIELDNAM);
	labAxis += " )";
	change_xref(FIELDNAM, labAxis);
}

//-----
DvObject_var DvObject::tan(){

	if(is_dbl() ){
		DvObject_var res = new DvObject(this); // constructor copies xrefs
		res->tanThis();
		return res;
	}
	else return new DvObject();
}
void DvObject::tanThis(){	
	
	if(is_dbl()){
		// convert this to radians
		toRadThis();
		
		for(size_t i=0; i<totalSize(); i++) Ddata[i] = std::tan(Ddata[i]);
	}
	else{
		error("[tan] not possible on data type");
		return;
	}
	
	// fix attributes
	DvString SIC("1.0>(tan)");
	change_xref(SI_CONVERSION, SIC);

	DvString units("");
	change_xref(UNITS, units);

	DvString lab = this->getXrefText( "LABLAXIS" );
	if(lab.empty()) lab = this->getXrefText( "ObjectName" );
	
	DvString labAxis("tan( ");
	labAxis += lab;
	labAxis += " )";
	change_xref(LABLAXIS, labAxis);

	DvString fieldname("tan( ");
	labAxis += getXrefText(FIELDNAM);
	labAxis += " )";
	change_xref(FIELDNAM, labAxis);
}

//-----
DvObject_var DvObject::cot(){

	if(is_dbl() ){
		DvObject_var res = new DvObject(this); // constructor copies xrefs
		res->cotThis();
		return res;
	}
	else return new DvObject();
}
void DvObject::cotThis(){	
	
	if(is_dbl()){
		// convert this to radians
		toRadThis();
		double tanVal;
		for(size_t i=0; i<totalSize(); i++) {
			tanVal = std::tan(Ddata[i]);
			if(tanVal == 0.0) Ddata[i] = DvNaN;
			else Ddata[i] = 1.0/tanVal;
		}
	}
	else{
		error("[cot] not possible on data type");
		return;
	}
	
	// fix attributes
	DvString SIC("1.0>(cot)");
	change_xref(SI_CONVERSION, SIC);

	DvString units("");
	change_xref(UNITS, units);

	DvString lab = this->getXrefText( "LABLAXIS" );
	if(lab.empty()) lab = this->getXrefText( "ObjectName" );
	
	DvString labAxis("cot( ");
	labAxis += lab;
	labAxis += " )";
	change_xref(LABLAXIS, labAxis);

	DvString fieldname("cot( ");
	labAxis += getXrefText(FIELDNAM);
	labAxis += " )";
	change_xref(FIELDNAM, labAxis);
}

//-----
DvObject_var DvObject::sinh(){

	if(is_dbl() ){
		DvObject_var res = new DvObject(this); // constructor copies xrefs
		res->sinhThis();
		return res;
	}
	else return new DvObject();
}
void DvObject::sinhThis(){	
	
	if(is_dbl()){

		for(size_t i=0; i<totalSize(); i++) Ddata[i] = std::sinh(Ddata[i]);
	}
	else{
		error("[sinh] not possible on data type");
		return;
	}
	if( hasUnits() ){
		error("[sinh] object has units");
	}
		
	// fix attributes
	DvString SIC("1.0>(sinh)");
	change_xref(SI_CONVERSION, SIC);

	DvString units("");
	change_xref(UNITS, units);

	DvString lab = this->getXrefText( "LABLAXIS" );
	if(lab.empty()) lab = this->getXrefText( "ObjectName" );
	
	DvString labAxis("sinh( ");
	labAxis += lab;
	labAxis += " )";
	change_xref(LABLAXIS, labAxis);

	DvString fieldname("sinh( ");
	labAxis += getXrefText(FIELDNAM);
	labAxis += " )";
	change_xref(FIELDNAM, labAxis);

}

//-----
DvObject_var DvObject::cosh(){

	if(is_dbl() ){
		DvObject_var res = new DvObject(this); // constructor copies xrefs
		res->coshThis();
		return res;
	}
	else return new DvObject();
}
void DvObject::coshThis(){	
	
	if(is_dbl()){

		for(size_t i=0; i<totalSize(); i++) Ddata[i] = std::cosh(Ddata[i]);
	}
	else{
		error("[cosh] not possible on data type");
		return;
	}
	
	// fix attributes
	DvString SIC("1.0>(cosh)");
	change_xref(SI_CONVERSION, SIC);

	DvString units("");
	change_xref(UNITS, units);

	DvString lab = this->getXrefText( "LABLAXIS" );
	if(lab.empty()) lab = this->getXrefText( "ObjectName" );
	
	DvString labAxis("cosh( ");
	labAxis += lab;
	labAxis += " )";
	change_xref(LABLAXIS, labAxis);

	DvString fieldname("cosh( ");
	labAxis += getXrefText(FIELDNAM);
	labAxis += " )";
	change_xref(FIELDNAM, labAxis);
	
	if( hasUnits() ){
		error("[cosh] object has units");
	}
}

//-----
DvObject_var DvObject::tanh(){

	if(is_dbl() ){
		DvObject_var res = new DvObject(this); // constructor copies xrefs
		res->tanhThis();
		return res;
	}
	else return new DvObject();
}
void DvObject::tanhThis(){	
	
	if(is_dbl()){

		for(size_t i=0; i<totalSize(); i++) Ddata[i] = std::tanh(Ddata[i]);
	}
	else{
		error("[tanh] not possible on data type");
		return;
	}
	
	// fix attributes
	DvString SIC("1.0>(tanh)");
	change_xref(SI_CONVERSION, SIC);

	DvString units("");
	change_xref(UNITS, units);

	DvObject_var minVal = new DvObject(-1.0);
	change_xref(SCALEMIN, minVal);

	DvObject_var maxVal = new DvObject(1.0);
	change_xref(SCALEMAX, maxVal);

	DvString lab = this->getXrefText( "LABLAXIS" );
	if(lab.empty()) lab = this->getXrefText( "ObjectName" );
	
	DvString labAxis("tanh( ");
	labAxis += lab;
	labAxis += " )";
	change_xref(LABLAXIS, labAxis);

	DvString fieldname("tanh( ");
	labAxis += getXrefText(FIELDNAM);
	labAxis += " )";
	change_xref(FIELDNAM, labAxis);
	
	if( hasUnits() ){
		error("[tanh] object has units");
	}
}

//-----
DvObject_var DvObject::coth(){

	if(is_dbl() ){
		DvObject_var res = new DvObject(this); // constructor copies xrefs
		res->cothThis();
		return res;
	}
	else return new DvObject();
}
void DvObject::cothThis(){	
	
	if(is_dbl()){
		double tanhVal;
		for(size_t i=0; i<totalSize(); i++) {
			tanhVal = std::tanh(Ddata[i]);
			if(tanhVal == 0.0) Ddata[i] = DvNaN;
			else Ddata[i] = 1.0/tanhVal;
		}
	}
	else{
		error("[coth] not possible on data type");
		return;
	}
	
	// fix attributes
	DvString SIC("1.0>(coth)");
	change_xref(SI_CONVERSION, SIC);

	DvString units("");
	change_xref(UNITS, units);

	DvObject_var minVal = new DvObject(-1.0);
	change_xref(SCALEMIN, minVal);

	DvObject_var maxVal = new DvObject(1.0);
	change_xref(SCALEMAX, maxVal);

	DvString lab = this->getXrefText( "LABLAXIS" );
	if(lab.empty()) lab = this->getXrefText( "ObjectName" );
	
	DvString labAxis("coth( ");
	labAxis += lab;
	labAxis += " )";
	change_xref(LABLAXIS, labAxis);

	DvString fieldname("coth( ");
	labAxis += getXrefText(FIELDNAM);
	labAxis += " )";
	change_xref(FIELDNAM, labAxis);
	
	if( hasUnits() ){
		error("[coth] object has units");
	}
}

// -----
DvObject_var DvObject::asin(bool indeg){

	if(is_dbl() ){
		DvObject_var res = new DvObject(this); // constructor copies xrefs
		res->asinThis(indeg);
		return res;
	}
	else return new DvObject();
}
void DvObject::asinThis(bool indeg){	
	
	if(is_dbl()){
		double rad2deg = 180./DV_PI;
		
		// get angle (radian)
		if(indeg) {
			for(size_t i=0; i<totalSize(); i++) {
				double val = Ddata[i];
				if( val < -1.0 || val > 1.0) Ddata[i] = DvNaN;
				else Ddata[i] = std::asin(val) * rad2deg;
			}
		}
		else {
			for(size_t i=0; i<totalSize(); i++) {
				double val = Ddata[i];
				if( val < -1.0 || val > 1.0) Ddata[i] = DvNaN;
				else Ddata[i] = std::asin(val);
			}
		}
	}
	else{
		error("[asin] not possible on data type");
		return;
	}
	
	// fix attributes
	if(indeg){
		DvString SIC("1.0>deg");
		change_xref(SI_CONVERSION, SIC);

		DvString units("deg");
		change_xref(UNITS, units);

		DvObject_var minVal = new DvObject(-90.0);
		change_xref(SCALEMIN, minVal);

		DvObject_var maxVal = new DvObject(90.0);
		change_xref(SCALEMAX, maxVal);
	}
	else{
		DvString SIC("1.0>rad");
		change_xref(SI_CONVERSION, SIC);

		DvString units("rad");
		change_xref(UNITS, units);

		DvObject_var minVal = new DvObject(-DV_PI/2.);
		change_xref(SCALEMIN, minVal);

		DvObject_var maxVal = new DvObject(DV_PI/2.);
		change_xref(SCALEMAX, maxVal);
	}
	

	DvString lab = this->getXrefText( "LABLAXIS" );
	if(lab.empty()) lab = this->getXrefText( "ObjectName" );
	
	DvString labAxis("asin( ");
	labAxis += lab;
	labAxis += " )";
	change_xref(LABLAXIS, labAxis);

	DvString fieldname("asin( ");
	labAxis += getXrefText(FIELDNAM);
	labAxis += " )";
	change_xref(FIELDNAM, labAxis);

	if( hasUnits() ){
		error("[asin] object has units");
	}
}

// -----
DvObject_var DvObject::acos(bool indeg){

	if(is_dbl() ){
		DvObject_var res = new DvObject(this); // constructor copies xrefs
		res->acosThis(indeg);
		return res;
	}
	else return new DvObject();
}
void DvObject::acosThis(bool indeg){	
	
	if(is_dbl()){
		double rad2deg = 180./DV_PI;
		
		// get angle (radian)
		if(indeg) {
			for(size_t i=0; i<totalSize(); i++) {
				double val = Ddata[i];
				if( val < -1.0 || val > 1.0) Ddata[i] = DvNaN;
				else Ddata[i] = std::acos(val) * rad2deg;
			}
		}
		else {
			for(size_t i=0; i<totalSize(); i++) {
				double val = Ddata[i];
				if( val < -1.0 || val > 1.0) Ddata[i] = DvNaN;
				else Ddata[i] = std::acos(val);
			}
		}
	}
	else{
		error("[acos] not possible on data type");
		return;
	}
	
	// fix attributes
	if(indeg){
		DvString SIC("1.0>deg");
		change_xref(SI_CONVERSION, SIC);

		DvString units("deg");
		change_xref(UNITS, units);

		DvObject_var minVal = new DvObject(-90.0);
		change_xref(SCALEMIN, minVal);

		DvObject_var maxVal = new DvObject(90.0);
		change_xref(SCALEMAX, maxVal);
	}
	else{
		DvString SIC("1.0>rad");
		change_xref(SI_CONVERSION, SIC);

		DvString units("rad");
		change_xref(UNITS, units);

		DvObject_var minVal = new DvObject(-DV_PI/2.);
		change_xref(SCALEMIN, minVal);

		DvObject_var maxVal = new DvObject(DV_PI/2.);
		change_xref(SCALEMAX, maxVal);
	}
	

	DvString lab = this->getXrefText( "LABLAXIS" );
	if(lab.empty()) lab = this->getXrefText( "ObjectName" );
	
	DvString labAxis("acos( ");
	labAxis += lab;
	labAxis += " )";
	change_xref(LABLAXIS, labAxis);

	DvString fieldname("acos( ");
	labAxis += getXrefText(FIELDNAM);
	labAxis += " )";
	change_xref(FIELDNAM, labAxis);
	
	if( hasUnits() ){
		error("[acos] object has units");
	}
}

// -----
DvObject_var DvObject::atan(bool indeg){

	if(is_dbl() ){
		DvObject_var res = new DvObject(this); // constructor copies xrefs
		res->atanThis(indeg);
		return res;
	}
	else return new DvObject();
}
void DvObject::atanThis(bool indeg){	
	
	if(is_dbl()){
		double rad2deg = 180./DV_PI;
		
		// get angle (radian)
		if(indeg) {
			for(size_t i=0; i<totalSize(); i++) {
				Ddata[i] = std::atan( Ddata[i]) * rad2deg;
			}
		}
		else {
			for(size_t i=0; i<totalSize(); i++) {
				Ddata[i] = std::atan(Ddata[i]);
			}
		}
	}
	else{
		error("[atan] not possible on data type");
		return;
	}
	
	// fix attributes
	if(indeg){
		DvString SIC("1.0>deg");
		change_xref(SI_CONVERSION, SIC);

		DvString units("deg");
		change_xref(UNITS, units);

		DvObject_var minVal = new DvObject(-90.0);
		change_xref(SCALEMIN, minVal);

		DvObject_var maxVal = new DvObject(90.0);
		change_xref(SCALEMAX, maxVal);
	}
	else{
		DvString SIC("1.0>rad");
		change_xref(SI_CONVERSION, SIC);

		DvString units("rad");
		change_xref(UNITS, units);

		DvObject_var minVal = new DvObject(-DV_PI/2.);
		change_xref(SCALEMIN, minVal);

		DvObject_var maxVal = new DvObject(DV_PI/2.);
		change_xref(SCALEMAX, maxVal);
	}
	

	DvString lab = this->getXrefText( "LABLAXIS" );
	if(lab.empty()) lab = this->getXrefText( "ObjectName" );
	
	DvString labAxis("atan( ");
	labAxis += lab;
	labAxis += " )";
	change_xref(LABLAXIS, labAxis);

	DvString fieldname("atan( ");
	labAxis += getXrefText(FIELDNAM);
	labAxis += " )";
	change_xref(FIELDNAM, labAxis);
	
	if( hasUnits() ){
		error("[atan] object has units");
	}
}

// -----
DvObject_var DvObject::atan2(DvObject_var &x, bool indeg){
	// atan(*this / x) in range [-pi,pi]
	
	if(totalSize() != x->totalSize() && x->totalSize() != 1){
		error("[atan2] y and x must have same number of elements, or x a single value");
		return new DvObject();
	}
	if(is_dbl() ){
		DvObject_var res = new DvObject(this); // constructor copies xrefs
		res->atan2This(x, indeg);
		return res;
	}
	else return new DvObject();
}
void DvObject::atan2This(DvObject_var &x, bool indeg){	
	// atan(*this / x) in range [-pi,pi]
	
	if(is_dbl() && x->is_dbl()){
		double rad2deg = 180./DV_PI;
		if(x->totalSize() == 1){
			// get angle (radian)
			double xVal = x->Ddata[0];
			if(indeg) {
				for(size_t i=0; i<totalSize(); i++) {
					Ddata[i] = std::atan2( Ddata[i], xVal) * rad2deg;
				}
			}
			else {
				for(size_t i=0; i<totalSize(); i++) {
					Ddata[i] = std::atan2(Ddata[i], xVal);
				}
			}
		}
		else{
			// get angle (radian)
			if(indeg) {
				for(size_t i=0; i<totalSize(); i++) {
					Ddata[i] = std::atan2( Ddata[i], x->Ddata[i]) * rad2deg;
				}
			}
			else {
				for(size_t i=0; i<totalSize(); i++) {
					Ddata[i] = std::atan2(Ddata[i], x->Ddata[i]);
				}
			}
		}
	}
	else{
		error("[atan2] not possible on data type");
		return;
	}
	
	// fix attributes
	if(indeg){
		DvString SIC("1.0>deg");
		change_xref(SI_CONVERSION, SIC);

		DvString units("deg");
		change_xref(UNITS, units);

		DvObject_var minVal = new DvObject(-180.0);
		change_xref(SCALEMIN, minVal);

		DvObject_var maxVal = new DvObject(180.0);
		change_xref(SCALEMAX, maxVal);
	}
	else{
		DvString SIC("1.0>rad");
		change_xref(SI_CONVERSION, SIC);

		DvString units("rad");
		change_xref(UNITS, units);

		DvObject_var minVal = new DvObject(-DV_PI);
		change_xref(SCALEMIN, minVal);

		DvObject_var maxVal = new DvObject(DV_PI);
		change_xref(SCALEMAX, maxVal);
	}
	

	DvString lab = this->getXrefText( "LABLAXIS" );
	if(lab.empty()) lab = this->getXrefText( "ObjectName" );
	
	DvString labAxis("atan2( ");
	labAxis += lab;
	labAxis += " )";
	change_xref(LABLAXIS, labAxis);

	DvString fieldname("atan2( ");
	labAxis += getXrefText(FIELDNAM);
	labAxis += " )";
	change_xref(FIELDNAM, labAxis);
	
}

// -----
DvObject_var DvObject::acot(bool indeg){

	if(is_dbl() ){
		DvObject_var res = new DvObject(this); // constructor copies xrefs
		res->acotThis(indeg);
		return res;
	}
	else return new DvObject();
}
void DvObject::acotThis(bool indeg){	
	
	if(is_dbl()){
		double rad2deg = 180./DV_PI;
		
		// get angle (radian)
		if(indeg) {
			for(size_t i=0; i<totalSize(); i++) {
				if(Ddata[i] == 0) Ddata[i] = 90.0;
				else Ddata[i] = std::atan( 1.0 / Ddata[i]) * rad2deg;
			}
		}
		else {
			for(size_t i=0; i<totalSize(); i++) {
				if(Ddata[i] == 0) Ddata[i] = DV_PI * 0.5;
				else Ddata[i] = std::atan( 1.0 / Ddata[i]);
			}
		}
	}
	else{
		error("[acot] not possible on data type");
		return;
	}
	
	// fix attributes
	if(indeg){
		DvString SIC("1.0>deg");
		change_xref(SI_CONVERSION, SIC);

		DvString units("deg");
		change_xref(UNITS, units);

		DvObject_var minVal = new DvObject(-90.0);
		change_xref(SCALEMIN, minVal);

		DvObject_var maxVal = new DvObject(90.0);
		change_xref(SCALEMAX, maxVal);
	}
	else{
		DvString SIC("1.0>rad");
		change_xref(SI_CONVERSION, SIC);

		DvString units("rad");
		change_xref(UNITS, units);

		DvObject_var minVal = new DvObject(-DV_PI/2.);
		change_xref(SCALEMIN, minVal);

		DvObject_var maxVal = new DvObject(DV_PI/2.);
		change_xref(SCALEMAX, maxVal);
	}
	

	DvString lab = this->getXrefText( "LABLAXIS" );
	if(lab.empty()) lab = this->getXrefText( "ObjectName" );
	
	DvString labAxis("acot( ");
	labAxis += lab;
	labAxis += " )";
	change_xref(LABLAXIS, labAxis);

	DvString fieldname("acot( ");
	labAxis += getXrefText(FIELDNAM);
	labAxis += " )";
	change_xref(FIELDNAM, labAxis);
	
	if( hasUnits() ){
		error("[acot] object has units");
	}
}

// -----
DvObject_var DvObject::asinh(){

	if(is_dbl() ){
		DvObject_var res = new DvObject(this); // constructor copies xrefs
		res->asinhThis();
		return res;
	}
	else return new DvObject();
}
void DvObject::asinhThis(){	
	
	if(is_dbl()){
		
		for(size_t i=0; i<totalSize(); i++) {
			double val = Ddata[i];
			Ddata[i] = std::log(val + std::sqrt(val*val + 1));
		}
	}
	else{
		error("[asinh] not possible on data type");
		return;
	}
	
	// fix attributes

	DvString lab = this->getXrefText( "LABLAXIS" );
	if(lab.empty()) lab = this->getXrefText( "ObjectName" );
	
	DvString labAxis("asinh( ");
	labAxis += lab;
	labAxis += " )";
	change_xref(LABLAXIS, labAxis);

	DvString fieldname("asinh( ");
	labAxis += getXrefText(FIELDNAM);
	labAxis += " )";
	change_xref(FIELDNAM, labAxis);
}

// -----
DvObject_var DvObject::acosh(){

	if(is_dbl() ){
		DvObject_var res = new DvObject(this); // constructor copies xrefs
		res->acoshThis();
		return res;
	}
	else return new DvObject();
}
void DvObject::acoshThis(){	
	
	if(is_dbl()){
		
		for(size_t i=0; i<totalSize(); i++) {
			double val = Ddata[i];
			if(val <= 1.0) Ddata[i] = DvNaN;
			else {
				Ddata[i] = std::log(val + std::sqrt(val*val - 1));
			}
		}
	}
	else{
		error("[acosh] not possible on data type");
		return;
	}
	
	// fix attributes

	DvString lab = this->getXrefText( "LABLAXIS" );
	if(lab.empty()) lab = this->getXrefText( "ObjectName" );
	
	DvString labAxis("acosh( ");
	labAxis += lab;
	labAxis += " )";
	change_xref(LABLAXIS, labAxis);

	DvString fieldname("acosh( ");
	labAxis += getXrefText(FIELDNAM);
	labAxis += " )";
	change_xref(FIELDNAM, labAxis);
}

// -----
DvObject_var DvObject::atanh(){

	if(is_dbl() ){
		DvObject_var res = new DvObject(this); // constructor copies xrefs
		res->atanhThis();
		return res;
	}
	else return new DvObject();
}
void DvObject::atanhThis(){	
	
	if(is_dbl()){
		
		for(size_t i=0; i<totalSize(); i++) {
			double val = Ddata[i];
			if(val <= -1.0 || val >= 1.0 ) Ddata[i] = DvNaN;
			else {
				Ddata[i] = 0.5 * std::log( (1+val)/(1-val));
			}
		}
	}
	else{
		error("[atanh] not possible on data type");
		return;
	}
	
	// fix attributes

	DvString lab = this->getXrefText( "LABLAXIS" );
	if(lab.empty()) lab = this->getXrefText( "ObjectName" );
	
	DvString labAxis("atanh( ");
	labAxis += lab;
	labAxis += " )";
	change_xref(LABLAXIS, labAxis);

	DvString fieldname("atanh( ");
	labAxis += getXrefText(FIELDNAM);
	labAxis += " )";
	change_xref(FIELDNAM, labAxis);
}

// -----
DvObject_var DvObject::acoth(){

	if(is_dbl() ){
		DvObject_var res = new DvObject(this); // constructor copies xrefs
		res->acothThis();
		return res;
	}
	else return new DvObject();
}
void DvObject::acothThis(){	
	
	if(is_dbl()){
		
		for(size_t i=0; i<totalSize(); i++) {
			double val = Ddata[i];
			if(val < -1.0 || val > 1.0 ){
				Ddata[i] = 0.5 * std::log( (val+1)/(val-1));
			}
			else Ddata[i] = DvNaN;
		}
	}
	else{
		error("[acoth] not possible on data type");
		return;
	}
	
	// fix attributes

	DvString lab = this->getXrefText( "LABLAXIS" );
	if(lab.empty()) lab = this->getXrefText( "ObjectName" );
	
	DvString labAxis("acoth( ");
	labAxis += lab;
	labAxis += " )";
	change_xref(LABLAXIS, labAxis);

	DvString fieldname("acoth( ");
	labAxis += getXrefText(FIELDNAM);
	labAxis += " )";
	change_xref(FIELDNAM, labAxis);
}


