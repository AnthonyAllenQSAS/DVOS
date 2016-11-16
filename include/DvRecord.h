//
//  DvRecord.h
//
//  Dvos -  Data variable object classes
//
//  Tony Allen
//  Apr 2013
//

#ifndef DVREC_H
#define DVREC_H


#include "DvObj.h"


//-----------------------------------------------------------------------------
// class DvRecord
//-----------------------------------------------------------------------------

class DvRecord  {
	DvObject *obj;
	size_t rec;
	bool isDouble;

  friend class DvObject;
  	
  public:	
	DvRecord(){ 
		obj = NULL; 
		rec = 0;
		isDouble = true;
	}
	
	DvRecord(DvObject *dobj, size_t drec=0){ 
		obj = dobj; 
		// single value (or array) will behave like sequence of any length with first val repeated
		if( drec < dobj->seqSize() ) rec = drec;
		else rec = 0; 
		isDouble = (obj->Ddata.size() > 0);
	}
	
	~DvRecord(){;};
	
	//  operators are implemented where applicable
	// =
	DvRecord &operator=(double val){ 
		size_t arrLen = obj->arraySize(); 
		slice sl(rec*arrLen, arrLen, 1); 
		if(isDouble) obj->Ddata[sl] = val; 
		else if(obj->is_time()) obj->Tdata[sl] = val; 
		else if(obj->is_int()) obj->Idata[sl] = (int) val; 
		else if(obj->is_str()) obj->Sdata[sl] = DvString(val); 
		return *this;
	}
	DvRecord &operator=(valarray<double> &val){ 
		size_t arrLen = obj->arraySize();
		slice sl(rec*arrLen, arrLen, 1);
		if(val.size() == arrLen && isDouble)  obj->Ddata[sl] = val; 
		return *this;
	}
	DvRecord &operator=(int val){ 
		size_t arrLen = obj->arraySize(); 
		slice sl(rec*arrLen, arrLen, 1); 
		if(obj->is_int()) obj->Idata[sl] = val; 
		else if(isDouble) obj->Ddata[sl] = val; 
		else if(obj->is_time()) obj->Tdata[sl] = val; 
		else if(obj->is_str()) obj->Sdata[sl] = DvString(val); 
		return *this;
	}
	DvRecord &operator=(valarray<int> &val){ 
		size_t arrLen = obj->arraySize();
		slice sl(rec*arrLen, arrLen, 1);
		if(val.size() == arrLen && obj->is_int()) obj->Idata[sl] = val; 
		return *this;
	}
	DvRecord &operator=(DvString val){ 
		size_t arrLen = obj->arraySize(); 
		slice sl(rec*arrLen, arrLen, 1); 
		if(obj->is_str()) obj->Sdata[sl] = val;  
		else if(isDouble) obj->Ddata[sl] = val.toDouble();  
		else if(obj->is_int()) obj->Idata[sl] = val.toInt();  
		else if(obj->is_time()) obj->Tdata[sl] = val.toTime();  
		else if(obj->is_event()) obj->Edata[sl] = val.toEvent();  
		return *this;
	}
	DvRecord &operator=(valarray<DvString> &val){ 
		size_t arrLen = obj->arraySize();
		slice sl(rec*arrLen, arrLen, 1);
		if(val.size() == arrLen && obj->is_str()) obj->Sdata[sl] = val; 
		return *this;
	}
	DvRecord &operator=(DvTime val){ 
		size_t arrLen = obj->arraySize(); 
		slice sl(rec*arrLen, arrLen, 1); 
		if(obj->is_time()) obj->Tdata[sl] = val;  
		else if(isDouble) obj->Ddata[sl] = val.getEpoch2000();  
		else if(obj->is_str()) obj->Sdata[sl] = val;  
		return *this;
	}
	DvRecord &operator=(valarray<DvTime> &val){ 
		size_t arrLen = obj->arraySize();
		slice sl(rec*arrLen, arrLen, 1);
		if(val.size() == arrLen && obj->is_time()) obj->Tdata[sl] = val; 
		return *this;
	}
	DvRecord &operator=(DvEvent val){ 
		size_t arrLen = obj->arraySize(); 
		slice sl(rec*arrLen, arrLen, 1); 
		if(obj->is_event()) obj->Edata[sl] = val;  
		else if(obj->is_str()) obj->Sdata[sl] = val;  
		return *this;
	}
	DvRecord &operator=(valarray<DvEvent> &val){ 
		size_t arrLen = obj->arraySize();
		slice sl(rec*arrLen, arrLen, 1);
		if(val.size() == arrLen && obj->is_event()) obj->Edata[sl] = val; 
		return *this;
	}

	
	// +=
	DvRecord &operator+=(double val){ 
		size_t arrLen = obj->arraySize(); 
		size_t ist = rec*arrLen;
		if(isDouble) for(size_t i=0; i<arrLen; i++) obj->Ddata[ist+i] += val; 
		else if(obj->is_int()) for(size_t i=0; i<arrLen; i++) obj->Idata[ist+i] += (int) val; 
		else if(obj->is_str()) for(size_t i=0; i<arrLen; i++) obj->Sdata[ist+i] += val; 
		else if(obj->is_time()) for(size_t i=0; i<arrLen; i++) obj->Tdata[ist+i] += val; 
		return *this;
	}
	DvRecord &operator+=(valarray<double> &val){ 
		size_t arrLen = obj->arraySize();
		slice sl(rec*arrLen, arrLen, 1);
		if(val.size() == arrLen && isDouble) obj->Ddata[sl] += val; 
		return *this;
	}
	DvRecord &operator+=(int val){ 
		size_t arrLen = obj->arraySize();  
		size_t ist = rec*arrLen;
		if(isDouble) for(size_t i=0; i<arrLen; i++) obj->Ddata[ist+i] += (double) val; 
		else if(obj->is_int()) for(size_t i=0; i<arrLen; i++) obj->Idata[ist+i] += val; 
		else if(obj->is_str()) for(size_t i=0; i<arrLen; i++) obj->Sdata[ist+i] += val; 
		else if(obj->is_time()) for(size_t i=0; i<arrLen; i++) obj->Tdata[ist+i] += val; 
		return *this;
	}
	DvRecord &operator+=(valarray<int> &val){ 
		size_t arrLen = obj->arraySize();
		slice sl(rec*arrLen, arrLen, 1);
		if(val.size() == arrLen && obj->is_int()) obj->Idata[sl] += val; 
		return *this;
	}
	DvRecord &operator+=(DvString &val){ 
		size_t arrLen = obj->arraySize();  
		size_t ist = rec*arrLen;
		if(obj->is_str()) for(size_t i=0; i<arrLen; i++) obj->Sdata[ist+i] += val; 
		else if(isDouble) for(size_t i=0; i<arrLen; i++) obj->Ddata[ist+i] += val.toDouble(); 
		else if(obj->is_int()) for(size_t i=0; i<arrLen; i++) obj->Idata[ist+i] += val.toInt(); 
		else if(obj->is_time()) for(size_t i=0; i<arrLen; i++) obj->Tdata[ist+i] += val.toDouble(); // add seconds
		return *this;
	}
	DvRecord &operator+=(valarray<DvString> &val){ 
		size_t arrLen = obj->arraySize();
		slice sl(rec*arrLen, arrLen, 1);
		if(val.size() == arrLen && obj->is_str()) obj->Sdata[sl] += val; 
		return *this;
	}
	// +=
	DvRecord &operator+=(DvTime val){ 
		size_t arrLen = obj->arraySize(); 
		size_t ist = rec*arrLen;
		if(obj->is_time()) for(size_t i=0; i<arrLen; i++) obj->Tdata[ist+i] += val; 
		return *this;
	}
	
	// -=
	DvRecord &operator-=(double val){ 
		size_t arrLen = obj->arraySize();  
		size_t ist = rec*arrLen;
		if(isDouble) for(size_t i=0; i<arrLen; i++) obj->Ddata[ist+i] -= val;  
		else if(obj->is_int()) for(size_t i=0; i<arrLen; i++) obj->Idata[ist+i] -= (int) val;  
		else if(obj->is_time()) for(size_t i=0; i<arrLen; i++) obj->Tdata[ist+i] -= val;  
		return *this;
	}
	DvRecord &operator-=(valarray<double> &val){ 
		size_t arrLen = obj->arraySize();
		slice sl(rec*arrLen, arrLen, 1);
		if(val.size() == arrLen && isDouble) obj->Ddata[sl] -= val; 
		return *this;
	}
	DvRecord &operator-=(int val){ 
		size_t arrLen = obj->arraySize();  
		size_t ist = rec*arrLen;
		if(isDouble) for(size_t i=0; i<arrLen; i++) obj->Ddata[ist+i] -= val;  
		if(obj->is_int()) for(size_t i=0; i<arrLen; i++) obj->Idata[ist+i] -= val;  
		if(obj->is_time()) for(size_t i=0; i<arrLen; i++) obj->Tdata[ist+i] -= val;  
		return *this;
	}
	DvRecord &operator-=(DvTime val){ 
		size_t arrLen = obj->arraySize();  
		size_t ist = rec*arrLen;
		if(obj->is_time()) for(size_t i=0; i<arrLen; i++) obj->Tdata[ist+i] -= val;  
		return *this;
	}
	DvRecord &operator-=(valarray<int> &val){ 
		size_t arrLen = obj->arraySize();
		slice sl(rec*arrLen, arrLen, 1);
		if(val.size() == arrLen && obj->is_int()) obj->Idata[sl] -= val; 
		return *this;
	}
	
	// *=
	DvRecord &operator*=(double val){ 
		size_t arrLen = obj->arraySize(); 
		size_t ist = rec*arrLen;
		if(isDouble) for(size_t i=0; i<arrLen; i++) obj->Ddata[ist+i] *= val;  
		else if(obj->is_int()) for(size_t i=0; i<arrLen; i++) obj->Idata[ist+i] *= (int) val;  
		else if(obj->is_time()) for(size_t i=0; i<arrLen; i++) obj->Tdata[ist+i] *= val;  
		return *this;
	}
	DvRecord &operator*=(valarray<double> &val){ // this is NOT matrix multiplication
		size_t arrLen = obj->arraySize();
		slice sl(rec*arrLen, arrLen, 1);
		if(val.size() == arrLen && isDouble) obj->Ddata[sl] *= val; 
		return *this;
	}
	DvRecord &operator*=(int val){ 
		size_t arrLen = obj->arraySize(); 
		size_t ist = rec*arrLen;
		if(isDouble) for(size_t i=0; i<arrLen; i++) obj->Ddata[ist+i] *= val;  
		if(obj->is_int()) for(size_t i=0; i<arrLen; i++) obj->Idata[ist+i] *= val;  
		return *this;
	}
	DvRecord &operator*=(valarray<int> &val){ // this is NOT matrix multiplication
		size_t arrLen = obj->arraySize();
		slice sl(rec*arrLen, arrLen, 1);
		if(val.size() == arrLen && obj->is_int()) obj->Idata[sl] *= val; 
		return *this;
	}
	
	// /=
	DvRecord &operator/=(double val){ 
		size_t arrLen = obj->arraySize(); 
		size_t ist = rec*arrLen;
		if(isDouble) for(size_t i=0; i<arrLen; i++) obj->Ddata[ist+i] /= val; 
		if(obj->is_int()) for(size_t i=0; i<arrLen; i++) obj->Idata[ist+i] /= (int) val; 
		if(obj->is_time()) for(size_t i=0; i<arrLen; i++) obj->Tdata[ist+i] /= val; 
		return *this;
	}
	DvRecord &operator/=(valarray<double> &val){ 
		size_t arrLen = obj->arraySize();
		slice sl(rec*arrLen, arrLen, 1);
		if(val.size() == arrLen && isDouble) obj->Ddata[sl] /= val; 
		return *this;
	}
	DvRecord &operator/=(int val){  // This will round down to int
		size_t arrLen = obj->arraySize(); 
		size_t ist = rec*arrLen;
		if(isDouble) for(size_t i=0; i<arrLen; i++) obj->Ddata[ist+i] /= val;
		if(obj->is_int()) for(size_t i=0; i<arrLen; i++) obj->Idata[ist+i] /= val; 
		return *this;
	}
	DvRecord &operator/=(valarray<int> &val){  // This will round down to int
		size_t arrLen = obj->arraySize();
		slice sl(rec*arrLen, arrLen, 1);
		if(val.size() == arrLen && obj->is_int()) obj->Idata[sl] /= val; 
		return *this;
	}

	
	// Generic =
	DvRecord &operator=(DvRecord val){ 
		int arrLen = obj->arraySize();
		slice sl(rec*arrLen, arrLen, 1);
		arrLen = val.obj->arraySize();
		slice sl_val(val.rec*arrLen, arrLen, 1);
		if(obj->is_dbl() && val.obj->is_dbl()) {
        obj->Ddata[sl] = val.obj->Ddata[sl_val];
        } 
		else if(obj->is_int() && val.obj->is_int()) obj->Idata[sl] = val.obj->Idata[sl_val]; 
		else if(obj->is_str() && val.obj->is_str()) obj->Sdata[sl] = val.obj->Sdata[sl_val]; 
		else if(obj->is_time() && val.obj->is_time()) obj->Tdata[sl] = val.obj->Tdata[sl_val]; 
		else if(obj->is_event() && val.obj->is_event()) obj->Edata[sl] = val.obj->Edata[sl_val]; 
		
		return *this;
	}


	// Generic +=
	DvRecord &operator+=(DvRecord val){ 
		size_t arrLen = obj->arraySize();
		if(arrLen != val.obj->arraySize()) {
			obj->error("DvRecord += operator array length mismatch");
			return *this;
		}
		// += not working under Mac		
		// slice sl(rec*arrLen, arrLen, 1);
		// slice sl_val(val.rec*arrLen, arrLen, 1);

		// if(obj->is_dbl() && val.obj->is_dbl()) obj->Ddata[sl] +=  val.obj->Ddata[sl_val]; 
		// else if(obj->is_int() && val.obj->is_int()) obj->Idata[sl] += val.obj->Idata[sl_val]; 
		// else if(obj->is_str() && val.obj->is_str()) obj->Sdata[sl] += val.obj->Sdata[sl_val]; 
		// else if(obj->is_time() && val.obj->is_time()) obj->Tdata[sl] += val.obj->Tdata[sl_val]; 
	
		if(obj->is_dbl() && val.obj->is_dbl()) {
			for(size_t k=0; k<arrLen; k++) obj->Ddata[rec*arrLen+k] +=  val.obj->Ddata[val.rec*arrLen+k]; 
		}
		else if(obj->is_int() && val.obj->is_int()) {
			for(size_t k=0; k<arrLen; k++) obj->Idata[rec*arrLen+k] += val.obj->Idata[val.rec*arrLen+k]; 
		}
		else if(obj->is_str() && val.obj->is_str()){
			for(size_t k=0; k<arrLen; k++) obj->Sdata[rec*arrLen+k] += val.obj->Sdata[val.rec*arrLen+k]; 
		}
		else if(obj->is_time() && val.obj->is_time()){
			for(size_t k=0; k<arrLen; k++) obj->Tdata[rec*arrLen+k] += val.obj->Tdata[val.rec*arrLen+k]; 
		}
		
		return *this;
	}

	// Generic -=
	DvRecord &operator-=(DvRecord val){ 
		size_t arrLen = obj->arraySize();
		if(arrLen != val.obj->arraySize()) {
			obj->error("DvRecord -= operator array length mismatch");
			return *this;
		}
		
		// -= not working under Mac		
		// slice sl(rec*arrLen, arrLen, 1);
		// slice sl_val(val.rec*arrLen, arrLen, 1);
		// if(obj->is_dbl() && val.obj->is_dbl()) obj->Ddata[sl] -= val.obj->Ddata[sl_val]; 
		// else if(obj->is_int() && val.obj->is_int()) obj->Idata[sl] -= val.obj->Idata[sl_val]; 
		// else if(obj->is_time() && val.obj->is_time()) obj->Tdata[sl] -= val.obj->Tdata[sl_val]; 
		
		if(obj->is_dbl() && val.obj->is_dbl()) {
			for(size_t k=0; k<arrLen; k++) obj->Ddata[rec*arrLen+k] -=  val.obj->Ddata[val.rec*arrLen+k]; 
		}
		else if(obj->is_int() && val.obj->is_int()) {
			for(size_t k=0; k<arrLen; k++) obj->Idata[rec*arrLen+k] -= val.obj->Idata[val.rec*arrLen+k]; 
		}
		else if(obj->is_time() && val.obj->is_time()){
			for(size_t k=0; k<arrLen; k++) obj->Tdata[rec*arrLen+k] -= val.obj->Tdata[val.rec*arrLen+k]; 
		}

		return *this;
	}

	// Generic *=
	DvRecord &operator*=(DvRecord val){ 
		
		// should only get here for valid matrix multiplication
		// result will stay the same size 
		
		size_t common = val.obj->dims[0];
		size_t thisStride = this->obj->arraySize()/common;
		size_t argStride = val.obj->arraySize()/common;
		
		size_t newArraySize = thisStride * argStride;
		
		valarray <double> res(0.0, newArraySize);
		
		int arrLen = this->obj->arraySize();
		slice sl(rec*arrLen, arrLen, 1);
		
		arrLen = val.obj->arraySize();
		slice sl_val(val.rec*arrLen, arrLen, 1);
		
		valarray <double> v1 = obj->Ddata[sl];
		valarray <double> v2 = val.obj->Ddata[sl_val];
		// this should work for any dimension and rank 
		// so long as val is square (to keep result same dimensions) (for *= only)
		
		for(size_t i=0; i<thisStride; i++){
			for(size_t j=0; j<argStride; j++){
				for(size_t p=0; p<common; p++){ // sum down inner dims
					res[i*argStride+j] += v1[i*common+p] * v2[p*argStride+j]; 
				}
			}
		}
		
		// copy result into valarray
		this->obj->Ddata[sl] = res;
		

		return *this;
	}

	// Generic /=
	DvRecord &operator/=(DvRecord val){ 
		size_t arrLen = obj->arraySize();
		if(arrLen != val.obj->arraySize()) {
			obj->error("DvRecord /= operator array length mismatch");
			return *this;
		}
		
		// /= not working under Mac		
		// slice sl(rec*arrLen, arrLen, 1);
		// slice sl_val(val.rec*arrLen, arrLen, 1);
		// if(isDouble && val.obj->is_dbl()) obj->Ddata[sl] /= val.obj->Ddata[sl_val]; 
		
		if(isDouble && val.obj->is_dbl()) {
			for(size_t k=0; k<arrLen; k++) obj->Ddata[rec*arrLen+k] /=  val.obj->Ddata[val.rec*arrLen+k];  
		}
		return *this;
	}
	


};




#endif // #ifndef DVREC_H

