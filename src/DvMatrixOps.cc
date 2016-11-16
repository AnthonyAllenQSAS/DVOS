//
//  DvAlgorithms.cc
//
//  Dvos -  Data variable object classes
//
//  Algorithm calls 
//
//  Tony Allen
//  Apr 2013
//

#include <iostream>

#include "DvObj.h"
#include "DvRecord.h"


// These algorithms are not the fastest, but are chosen to be robust
// Faster algorithms should be added at some point
// However, note that they allow slices to operate on all
// records at once, whereas other methods for inverse etc. may not.



DvObject_var DvObject::matrixInverse(){
	
	// operates only on Ddata 
	// fast for long sequences of 3x3 matrices (e.g. rotations)
	// not so fast for large matrices, e.g. solving large sets of linear equations
	// Robust against the occasional singular matrix
	
	if(this->nDims() != 2) {
		DvObject_var res = new DvObject();
		res->error("Operation only available for 2D array");
		return res;
	}
	if(dims[0] != dims[1]){
		DvObject_var res = new DvObject();
		res->error("Operation only available for square array");
		return res;
	}
	if(not_dbl()){
		DvObject_var res = new DvObject();
		res->error("Operation only available for double type");
		return res;
	}


	DvObject_var adj = this->adjoint();
	
	DvObject_var det = this->det();
	
	
	adj->copy_xrefs_from(*this);
	det->copy_xref_from(DEPEND_0, *this);	
	
	DvMask msk(seqSize());
	bool hasSingular = false;
	for(size_t i=0; i<seqSize(); i++) {
		if( fabs(det->dbl(i)) < 1e-30 ) {
			msk[i] = false;
			hasSingular = true;
		}
	}

	if(hasSingular){	
		adj->apply_mask(msk);
		det->apply_mask(msk);
	}

	*adj /= *det;

	return adj;
	
}

valarray <double> DvObject::detValarray(){

	// internal use only, on Ddata 
	size_t arrayDim = dims[0];
	size_t len = seqSize();	
	
	valarray <double> mres(0.0, len); // scalar sequence of determinants
	
	int arrLen = arrayDim*arrayDim;
	
	if(arrayDim == 1){
			mres = Ddata[slice(0, len, 1)];
	}
	else if(arrayDim == 2){
		// hard wire for speed
		for(size_t n=0; n<len; n++){
			// step through sequence 
			int k = n*arrLen;
			mres[n] = Ddata[k]*Ddata[k+3] - Ddata[k+1]*Ddata[k+2];
		}
	}
	else if(arrayDim == 3){
		// hard wire for speed
		
		for(size_t n=0; n<len; n++){
			// step through sequence 
			int k = n*arrLen;
			mres[n] = Ddata[k]*(Ddata[k+4]*Ddata[k+8] - Ddata[k+5]*Ddata[k+7])
					- Ddata[k+3]*(Ddata[k+1]*Ddata[k+8] - Ddata[k+2]*Ddata[k+7])
					+ Ddata[k+6]*(Ddata[k+1]*Ddata[k+5] - Ddata[k+2]*Ddata[k+4]);
		}
	}
	else{
		// step along dim
		
		// does all records together
		for (size_t i=0; i<arrayDim; i+=2){

			DvObject_var m = subMatrix(i,0,false);
			valarray <double> d(m->detValarray());
			d *= elementValarray(i*arrayDim);
			mres += d;
		}
		for (size_t i=1; i<arrayDim; i+=2){
			
			DvObject_var m = subMatrix(i,0,false);
			valarray <double> d(m->detValarray());
			d *= elementValarray(i*arrayDim);
			mres -= d;
		}
	}
	
	return mres;
}


DvObject_var DvObject::det(){

	// operates only on Ddata 

	if( not_dbl() ){
		DvObject_var res = new DvObject();
		res->error("cannot form determinant of this data type");
		return res;
	}
	if(this->nDims() != 2) {
		DvObject_var res = new DvObject();
		res->error("Operation only available for 2D array");
		return res;
	}
	if(dims[0] != dims[1]){
		DvObject_var res = new DvObject();
		res->error("Operation only available for square array");
		return res;
	}
	
	
	size_t arrayDim = dims[0];
	
	valarray <double> d(this->detValarray());
	DvObject_var res = new DvObject(d);

	// SI_conversion
	DvString newSIC = getSICPower(arrayDim);
	if(!newSIC.empty()) res->change_xref(SI_CONVERSION, newSIC);
	
	return res;
}

DvObject_var DvObject::adjoint(){

	// operates only on Ddata 

	if( not_dbl() ){
		DvObject_var res = new DvObject();
		res->error("cannot form adjoint of this data type");
		return res;
	}
	if(this->nDims() != 2) {
		DvObject_var res = new DvObject();
		res->error("Adjoint only available for 2D array");
		return res;
	}
	if(dims[0] != dims[1]){
		DvObject_var res = new DvObject();
		res->error("Adjoint only available for square array");
		return res;
	}
	
	int arrayDim = dims[0];
	int arrLen = arraySize();
	int len = seqSize();
	
	// step along dim

	valarray <double> mres(0.0, len*arrLen);
	int isign = 1;
	int jsign = 1;
	for(size_t i=0; i<dims[0]; i++) {
		for(size_t j=0; j<dims[1]; j++) {
			// apply to all records at once
			DvObject_var sub = subMatrix(i,j,false);
			valarray <double> cofact(sub->detValarray());

			int sign = isign*jsign;
			cofact *= sign;
			mres[slice(j*dims[0]+i, len, arrLen)] = cofact[slice(0,len,1)] ; // this transposes
			jsign *= -1;
			
		}
		jsign = 1; // reset
		isign *= -1;
	}
	DvObject_var res = new DvObject(mres, dims, seqSize());
	// SI_conversion
	DvString newSIC = getSICPower(arrayDim-1);
	if(!newSIC.empty()) res->change_xref(SI_CONVERSION, newSIC);

    DvObject_var dep1 = this->get_xref(DEPEND_1);
    res->change_xref(DEPEND_1, dep1);
    
    DvObject_var dep2 = this->get_xref(DEPEND_2);
    res->change_xref(DEPEND_2, dep2);
    
    DvObject_var dep0 = this->get_xref(DEPEND_0);
    res->change_xref(DEPEND_0, dep0);
    
	return res;

}


DvObject_var DvObject::element(size_t i){

	DvObject_var res;
	// return scalar sequence of element at i
	if(is_dbl()){
		valarray <double> val(0.0, seqSize());
		val = Ddata[slice(i, seqSize(), arraySize() )];
		res = new DvObject(val);
		
	}
	else if(is_int()){
		valarray <int> val(0, seqSize());
		val = Idata[slice(i, seqSize(), arraySize() )];
		res = new DvObject(val);
		
	}
	else if(is_str()){
		valarray <DvString> val("", seqSize());
		val = Sdata[slice(i, seqSize(), arraySize() )];
		res = new DvObject(val);
		
	}
	else if(is_time()){
		valarray <DvTime> val(DvTime(), seqSize());
		val = Tdata[slice(i, seqSize(), arraySize() )];
		res = new DvObject(val);
		
	}
	else if(is_event()){
		valarray <DvEvent> val(DvEvent(), seqSize());
		val = Edata[slice(i, seqSize(), arraySize() )];
		res = new DvObject(val);
		
	}
	else return new DvObject();
	
	DvString SIC = getXrefText(SI_CONVERSION);
	if( !SIC.empty() ) res->change_xref(SI_CONVERSION, SIC);
	return res;
}



std::slice_array <double> DvObject::elementValarray(size_t i){

	return Ddata[slice(i, seqSize(), arraySize() )];
}


DvObject_var DvObject::transpose(bool with){

	if(this->nDims() != 2) {
		DvObject_var res = new DvObject();
		res->error("Transpose only available for 2D array");
		return res;
	}

	// swap dimensions
	vector <size_t> newDim;
	newDim.push_back(dims[1]);
	newDim.push_back(dims[0]);
	
	DvObject_var res;
	// return scalar sequence of element at i
	if(is_dbl()){
		res = new DvObject(0.0, newDim, seqSize());
		for(size_t i=0; i<dims[0]; i++) {
			for(size_t j=0; j<dims[1]; j++) {
				// apply to all records at once
				res->setSlice(slice(j*newDim[1]+i, seqSize(), arraySize()), Ddata[slice(i*dims[1]+j, seqSize(), arraySize() )]);
			}
		}
	}
	else if(is_int()){
		res = new DvObject((int)0, newDim, seqSize());
		for(size_t i=0; i<dims[0]; i++) {
			for(size_t j=0; j<dims[1]; j++) {
				// apply to all records at once
				res->setSlice(slice(j*newDim[1]+i, seqSize(), arraySize()), Idata[slice(i*dims[1]+j, seqSize(), arraySize() )]);
			}
		}
	}
	else if(is_str()){
		res = new DvObject(DvString(), newDim, seqSize());
		for(size_t i=0; i<dims[0]; i++) {
			for(size_t j=0; j<dims[1]; j++) {
				// apply to all records at once
				res->setSlice(slice(j*newDim[1]+i, seqSize(), arraySize()), Sdata[slice(i*dims[1]+j, seqSize(), arraySize() )]);
			}
		}
	}
	else if(is_time()){
		res = new DvObject(DvTime(), newDim, seqSize());
		for(size_t i=0; i<dims[0]; i++) {
			for(size_t j=0; j<dims[1]; j++) {
				// apply to all records at once
				res->setSlice(slice(j*newDim[1]+i, seqSize(), arraySize()), Tdata[slice(i*dims[1]+j, seqSize(), arraySize() )]);
			}
		}
	}
	else if(is_event()){
		res = new DvObject(DvEvent(), newDim, seqSize());
		for(size_t i=0; i<dims[0]; i++) {
			for(size_t j=0; j<dims[1]; j++) {
				// apply to all records at once
				res->setSlice(slice(j*newDim[1]+i, seqSize(), arraySize()), Edata[slice(i*dims[1]+j, seqSize(), arraySize() )]);
			}
		}
	}
	else res = new DvObject();
    if(with){
        res->copy_xrefs_from(*this);
        // swap depends
        DvObject_var dep1 = this->get_xref(DEPEND_1);
        res->change_xref(DEPEND_2, dep1);
        DvObject_var dep2 = this->get_xref(DEPEND_2);
        res->change_xref(DEPEND_1, dep2);
    }
    else{
        DvString SIC = getXrefText(SI_CONVERSION);
        if( !SIC.empty() ) res->change_xref(SI_CONVERSION, SIC);
    }
	return res;
}


DvObject_var DvObject::subMatrix(size_t row, size_t col, bool with){
	
	// subMatrix
	DvObject_var res;
	vector <size_t> newDim;
	newDim.push_back(dims[0]-1);
	newDim.push_back(dims[1]-1);
	size_t newArrSize = newDim[0]*newDim[1];
	
	// return matrix sequence of subMatrix at i,j (omitting row i and col j)
	
	if(is_dbl()){
		res = new DvObject(0.0, newDim, seqSize());
		size_t n=0;
		for(size_t i=0; i<dims[0]; i++) {
			for(size_t j=0; j<dims[1]; j++) {
				if(i != row && j != col){
					if(n >= newArrSize){
						res->error("array bound error in subMatrix()");
						return res;
					}
					// apply to all records at once
					res->setSlice(slice(n, seqSize(), newArrSize), Ddata[slice(i*dims[1]+j, seqSize(), arraySize() )]);
					n++;
				}
			}
		}
	}
	else if(is_int()){
		res = new DvObject((int)0, newDim, seqSize());
		size_t n=0;
		for(size_t i=0; i<dims[0]; i++) {
			for(size_t j=0; j<dims[1]; j++) {
				if(i != row && j != col){
					if(n >= newArrSize){
						res->error("array bound error in subMatrix()");
						return res;
					}
					// apply to all records at once
					res->setSlice(slice(n, seqSize(), newArrSize), Idata[slice(i*dims[1]+j, seqSize(), arraySize() )]);
					n++;
				}
			}
		}
	}
	else if(is_str()){
		res = new DvObject(DvString(), newDim, seqSize());
		size_t n=0;
		for(size_t i=0; i<dims[0]; i++) {
			for(size_t j=0; j<dims[1]; j++) {
				if(i != row && j != col){
					if(n >= newArrSize){
						res->error("array bound error in subMatrix()");
						return res;
					}
					// apply to all records at once
					res->setSlice(slice(n, seqSize(), newArrSize), Sdata[slice(i*dims[1]+j, seqSize(), arraySize() )]);
					n++;
				}
			}
		}
	}
	else if(is_time()){
		res = new DvObject(DvTime(), newDim, seqSize());
		size_t n=0;
		for(size_t i=0; i<dims[0]; i++) {
			for(size_t j=0; j<dims[1]; j++) {
				if(i != row && j != col){
					if(n >= newArrSize){
						res->error("array bound error in subMatrix()");
						return res;
					}
					// apply to all records at once
					res->setSlice(slice(n, seqSize(), newArrSize), Tdata[slice(i*dims[1]+j, seqSize(), arraySize() )]);
					n++;
				}
			}
		}
	}
	else if(is_event()){
		res = new DvObject(DvEvent(), newDim, seqSize());
		size_t n=0;
		for(size_t i=0; i<dims[0]; i++) {
			for(size_t j=0; j<dims[1]; j++) {
				if(i != row && j != col){
					if(n >= newArrSize){
						res->error("array bound error in subMatrix()");
						return res;
					}
					// apply to all records at once
					res->setSlice(slice(n, seqSize(), newArrSize), Edata[slice(i*dims[1]+j, seqSize(), arraySize() )]);
					n++;
				}
			}
		}
	}
	else res = new DvObject();
	
    if(with){

        res->copy_xrefs_from(*this);
        
        DvObject_var dep1 = get_xref(DEPEND_1);
        if(dep1.is_ok()){
            dep1 = dep1->subDimension(row);
            res->change_xref(DEPEND_1, dep1);
        }
        DvObject_var dep2 = get_xref(DEPEND_2);
        if(dep2.is_ok()){
            dep2 = dep2->subDimension(col);
            res->change_xref(DEPEND_2, dep2);
        }
        DvObject_var lab1 = get_xref(LABL_PTR_1);
        if(lab1.is_ok()){
            lab1 = lab1->subDimension(row);
            res->change_xref(LABL_PTR_1, lab1);
        }
        DvObject_var lab2 = get_xref(LABL_PTR_2);
        if(lab2.is_ok()){
            lab2 = lab2->subDimension(col);
            res->change_xref(LABL_PTR_2, lab2);
        }
    }
    else{
        DvString SIC = getXrefText(SI_CONVERSION);
        if( !SIC.empty() ) res->change_xref(SI_CONVERSION, SIC);
    }
	return res;
}


DvObject_var DvObject::subDimension(size_t s, bool with){
    DvObject_var ret;
    if(this->nDims() != 1) return ret;
    
    vector <size_t> newDim;
    newDim.push_back(this->dims[0] - 1);
    
    ret = this->create(newDim, this->seqSize());

    if(is_dbl()){
        for(size_t n=0; n<seqSize(); n++){
            size_t k=0;
            for(size_t i=0; i<dims[0]; i++){
                if(i != s) {
                    ret->dbl(n,k) = this->dbl(n,i);
                    k++;
                }
            }
        }
    }
    else if(is_int()){
        for(size_t n=0; n<seqSize(); n++){
            size_t k=0;
            for(size_t i=0; i<dims[0]; i++){
                if(i != s) {
                    ret->itg(n,k) = this->itg(n,i);
                    k++;
                }
            }
        }
    }
    else if(is_str()){
        for(size_t n=0; n<seqSize(); n++){
            size_t k=0;
            for(size_t i=0; i<dims[0]; i++){
                if(i != s) {
                    ret->str(n,k) = this->str(n,i);
                    k++;
                }
            }
        }
    }
    else if(is_time()){
        for(size_t n=0; n<seqSize(); n++){
            size_t k=0;
            for(size_t i=0; i<dims[0]; i++){
                if(i != s) {
                    ret->time(n,k) = this->time(n,i);
                    k++;
                }
            }
        }
    }
    else if(is_event()){
        for(size_t n=0; n<seqSize(); n++){
            size_t k=0;
            for(size_t i=0; i<dims[0]; i++){
                if(i != s) {
                    ret->event(n,k) = this->event(n,i);
                    k++;
                }
            }
        }
    }
    
    if(with) ret->copy_xrefs_from(*this);
    
    // recursively fix xrefs
    DvNode *xref = ret->first_xref();
    while(xref){	
        
        if(xref->obj()->nDims() == 1 && xref->obj()->Dims()[0] == this->dims[0]) {
            xref->obj() = xref->obj()->subDimension(s,with);
        }
        
        xref = xref->next;
    }
    
    return ret;
}

DvObject_var DvObject::trace(){

    DvObject_var ret;
	
	if(this->nDims() != 2) {
		ret = new DvObject();
		ret->error("[trace] requires 2D array");
		return ret;
	}

	if(dims[0] != dims[1]) {
		ret = new DvObject();
		ret->error("[trace] requires square array");
		return ret;
	}
	int arrayDim = dims[0];

	int len = this->seqSize();
	if(this->is_dbl() ){
		valarray <double> val(0.0, len);
		for(int n=0; n<len; n++) {
			for( int i=0; i<arrayDim; i++) val[n] += this->dbl(n, i, i);
		}
		
		ret = new DvObject(val);
	}
	else if(this->is_int() ){
		valarray <int> val(0, len);
		
		for(int n=0; n<len; n++) {
			for( int i=0; i<arrayDim; i++) val[n] += this->itg(n, i, i);
		}
		ret = new DvObject(val);
	}
	
	if(ret.is_nil()) {
		ret->error("[trace] not possible for data type");
		return ret;
	}

    // Add timeline (or DEPEND_0) 
	
	ret->copy_xref_from(DEPEND_0, *this);

	// useful xrefs
	
	ret->copy_xref_from(UNITS, *this);
	ret->copy_xref_from(SI_CONVERSION, *this);
	
	DvString xref = getXrefText(FIELDNAM);
	if(!xref.empty()){
		xref += " trace";
		ret->change_xref(FIELDNAM, xref);
	}
	
	DvString labl = "trace(";
	labl += this->getXrefText( LABLAXIS );
	labl += ")";
	ret->change_xref(LABLAXIS, labl);


 	// Frame
	ret->setFrameAttr("scalar>na" );
	
	return ret;
}

// ---------------------

DvObject_var DvObject::subsetDim(size_t d, size_t from, size_t to)
{
	if(from == to) return sliceDim(d, from);
	
	if(d >= dims.size()){
		DvObject *ret = new DvObject();
		ret->error("[subset dimension] requested dim > rank");
		return ret;
	}
	if( to >= dims[d] || from >= dims[d] ){
		DvObject *ret = new DvObject();
		ret->error("[subset dimension] requested range outside dimension");
		return ret;
	}
	
	DvObject_var result;
	if(to > from){ 
		vector <size_t> newDims = dims;
		size_t newDimLen = to - from + 1;
		newDims[d] = newDimLen;
	
		result = this->create(newDims, seqLen); 
	
	
		valarray <size_t> lengths(dims.size()+1); // +1 to include record steps
		valarray <size_t> strides(dims.size()+1);
		size_t arrLen = arraySize();
	
		size_t start = from;
		if(dims.size() > 0){
			lengths[dims.size()] = dims[dims.size()-1];
			strides[dims.size()] = 1;
		}
		for (size_t i=dims.size(); i>1; i--){
			lengths[i-1] = dims[i-2];
			strides[i-1] = dims[i-1]*strides[i];
			if(d < i-1) start *= dims[i-1];
		}
		lengths[d+1] = newDimLen; // now correct for subset dimension, d+1 now we include records
		lengths[0] = seqLen;
		strides[0] = arrLen;


		gslice gsl(start, lengths, strides);
 		
		if(is_dbl() )result->Ddata = Ddata[gsl];
		else if(is_int() )result->Idata = Idata[gsl];
		else if(is_str() )result->Sdata = Sdata[gsl];
		else if(is_time() )result->Tdata = Tdata[gsl];
		else if(is_event() )result->Edata = Edata[gsl];
	
	}
	else { 
		// Handles cyclic array slice
		// if from > to, get (from, N)+(0, to)
		vector <size_t> newDims = dims;
		size_t newDimLen = dims[d] - from + to + 1;
		newDims[d] = newDimLen;
	
		result = this->create(newDims, seqLen); 
	
		valarray <size_t> lengths(dims.size()+1); // +1 to include record steps
		valarray <size_t> strides(dims.size()+1);
		size_t arrLen = arraySize();

		// for gslice in result
		valarray <size_t> strides_new(dims.size()+1);
		size_t arrLen_new = result->arraySize();

		// start with (from,N)	
		size_t start = from;
		if(dims.size() > 0){
			lengths[dims.size()] = dims[dims.size()-1];
			strides[dims.size()] = 1;
		}
		for (size_t i=dims.size(); i>1; i--){
			lengths[i-1] = dims[i-2];
			strides[i-1] = dims[i-1]*strides[i];
			if(d < i-1) start *= dims[i-1];
		}
		lengths[d+1] = dims[d] - from; // now correct for subset dimension, d+1 now we include records
		lengths[0] = seqLen;
		strides[0] = arrLen;

	
		size_t start_new = 0;
		strides_new[newDims.size()] = 1;
		for (size_t i=newDims.size(); i>1; i--){
			strides_new[i-1] = newDims[i-1]*strides_new[i];
		}
		strides_new[0] = arrLen_new;

		gslice gsl(start, lengths, strides);		
		gslice gsl_new(start_new, lengths, strides_new);

		if(is_dbl() )result->Ddata[gsl_new] = Ddata[gsl];
		else if(is_int() )result->Idata[gsl_new] = Idata[gsl];
		else if(is_str() )result->Sdata[gsl_new] = Sdata[gsl];
		else if(is_time() )result->Tdata[gsl_new] = Tdata[gsl];
		else if(is_event() )result->Edata[gsl_new] = Edata[gsl];
	
		// Now fill in (0, to)
		start = 0;
		lengths[d+1] = to+1; // now correct for subset dimension, d+1 now we include records

		start_new = dims[d] - from;
		for (size_t i=newDims.size(); i>1; i--){
			if(d < i-1) start_new *= newDims[i-1];
		}
		gsl = gslice(start, lengths, strides);		
		gsl_new = gslice(start_new, lengths, strides_new);

		if(is_dbl())result->Ddata[gsl_new] = Ddata[gsl];
		else if(is_int())result->Idata[gsl_new] = Idata[gsl];
		else if(is_str())result->Sdata[gsl_new] = Sdata[gsl];
		else if(is_time())result->Tdata[gsl_new] = Tdata[gsl];
		else if(is_event())result->Edata[gsl_new] = Edata[gsl];

	}
	
	result->copy_xrefs_from(*this);
	
	DvString dep_i = "DEPEND_";
	dep_i.append(d+1);
	if(result->xref_exists(dep_i) ) {
		DvObject_var dxref = result->get_xref(dep_i);
		if(dxref->arraySize() > 1){
			dxref = dxref->subsetDim(0, from, to);
			result->change_xref(dep_i, dxref);
		}
	}
	
	DvString label_i = "LABEL_";
	label_i.append(d+1);
	if(result->xref_exists(label_i)) {
		DvObject_var lxref = result->get_xref(label_i);
		if(lxref->arraySize() > 1){
			lxref = lxref->subsetDim(0, from, to);
			result->change_xref(label_i, lxref);
		}
	}
	
	DvString rep_i = "REPRESENTATION_";
	rep_i.append(d+1);
	if(result->xref_exists(rep_i)) {
		DvObject_var rxref = result->get_xref(rep_i);
		if(rxref->arraySize() > 1){
			rxref = rxref->subsetDim(0, from, to);
			result->change_xref(rep_i, rxref);
		}
	}
	
	DvString delta_p = "DELTA_PLUS";
	if(result->xref_exists(delta_p)) {
		DvObject_var dxref = result->get_xref(delta_p);
		if(dxref->arraySize() > 1){
			dxref = dxref->subsetDim(0, from, to);
			result->change_xref(delta_p, dxref);
		}
	}
	
	DvString delta_m = "DELTA_MINUS";
	if(result->xref_exists(delta_m)) {
		DvObject_var dxref = result->get_xref(delta_m);
		if(dxref->arraySize() > 1){
			dxref = dxref->subsetDim(0, from, to);
			result->change_xref(delta_m, dxref);
		}
	}
		
	return result;
}

// ---------------------

DvObject_var DvObject::sliceDim(size_t d, size_t at)
{
	if(d >= dims.size()){
		DvObject_var ret = new DvObject();
		ret->error("[slice dimension] requested dim > rank");
		return ret;
	}
	if(at >= dims[d]){
		DvObject_var ret = new DvObject();
		ret->error("[slice dimension] requested row outside dimension");
		return ret;
	}

	vector <size_t> newDims = dims;
	newDims.erase(newDims.begin()+d);

	DvObject_var result = this->create(newDims, seqLen); 
	
	valarray <size_t> lengths(dims.size()+1); // +1 to include record steps
	valarray <size_t> strides(dims.size()+1);
	size_t arrLen = arraySize();
	
	size_t start = at;
	if(dims.size() > 0){
		lengths[dims.size()] = dims[dims.size()-1];
		strides[dims.size()] = 1;
	}
	for (size_t i=dims.size(); i>1; i--){
		lengths[i-1] = dims[i-2];
		strides[i-1] = dims[i-1]*strides[i];
		if(d < i-1) start *= dims[i-1];
	}
	lengths[d+1] = 1; // now correct for slice dimension, d+1 now we include records
	lengths[0] = seqLen;
	strides[0] = arrLen;


	gslice gsl(start, lengths, strides);
 		
	if(is_dbl() )result->Ddata = Ddata[gsl];
	else if(is_int() )result->Idata = Idata[gsl];
	else if(is_str() )result->Sdata = Sdata[gsl];
	else if(is_time() )result->Tdata = Tdata[gsl];
	else if(is_event() )result->Edata = Edata[gsl];


	result->copy_xrefs_from(*this);

	
	for(size_t i=d+1; i<=dims.size(); i++){
		DvString dep_i = "DEPEND_";
		dep_i.append(i);
		result->delete_xref(dep_i); // safe if xref does not exist
		// get next and move down in i
		DvString dep_ip1 = "DEPEND_";
		dep_ip1.append(i+1);
		if(result->xref_exists(dep_ip1)) {
			DvObject_var dxref = result->get_xref(dep_ip1);
			result->change_xref(dep_i, dxref);
		}
	}
	
	for(size_t i=d+1; i<=dims.size(); i++){
		DvString lab_i = "LABEL_";
		lab_i.append(i);
		result->delete_xref(lab_i); // safe if xref does not exist
		// get next and move down in i
		DvString lab_ip1 = "LABEL_";
		lab_ip1.append(i+1);
		if(result->xref_exists(lab_ip1)) {
			DvObject_var lxref = result->get_xref(lab_ip1);
			result->change_xref(lab_i, lxref);
		}
	}

	for(size_t i=d+1; i<=dims.size(); i++){
		DvString rep_i = "REPRESENTATION_";
		rep_i.append(i);
		result->delete_xref(rep_i); // safe if xref does not exist
		// get next and move down in i
		DvString rep_ip1 = "REPRESENTATION_";
		rep_ip1.append(i+1);
		if(result->xref_exists(rep_ip1)) {
			DvObject_var rxref = result->get_xref(rep_ip1);
			result->change_xref(rep_i, rxref);
		}
	}

    DvObject_var rank = result->get_xref("rank");
    if(rank.is_ok()){
        int rankVal = rank->asInt();
        rankVal -= 1;
        rank = new DvObject(rankVal);
        result->change_xref("rank", rank);
    }
    result->delete_xref(FRAME);

	return result;

}

void DvObject::setSliceDim(size_t d, size_t at, std::slice_array <double> val)
{
	// 
	
	if(d >= dims.size()){
		this->error("[slice dimension] requested dim > rank");
		return;
	}
	if(at >= dims[d]){
		this->error("[slice dimension] requested row outside dimension");
		return ;
	}
	
	valarray <size_t> lengths(dims.size()+1); // +1 to include record steps
	valarray <size_t> strides(dims.size()+1);
	size_t arrLen = arraySize();
	
	size_t start = at;
	if(dims.size() > 0){
		lengths[dims.size()] = dims[dims.size()-1];
		strides[dims.size()] = 1;
	}
	for (size_t i=dims.size(); i>1; i--){
		lengths[i-1] = dims[i-2];
		strides[i-1] = dims[i-1]*strides[i];
		if(d < i-1) start *= dims[i-1];
	}
	lengths[d+1] = 1; // now correct for slice dimension, d+1 now we include records
	lengths[0] = seqLen;
	strides[0] = arrLen;


	gslice gsl(start, lengths, strides);
 		
	this->Ddata[gsl] = val;

	return;
}

void DvObject::setSliceDim(size_t d, size_t at, std::slice_array <int> val)
{
	// 
	
	if(d >= dims.size()){
		this->error("[slice dimension] requested dim > rank");
		return;
	}
	if(at >= dims[d]){
		this->error("[slice dimension] requested row outside dimension");
		return ;
	}
	
	valarray <size_t> lengths(dims.size()+1); // +1 to include record steps
	valarray <size_t> strides(dims.size()+1);
	size_t arrLen = arraySize();
	
	size_t start = at;
	if(dims.size() > 0){
		lengths[dims.size()] = dims[dims.size()-1];
		strides[dims.size()] = 1;
	}
	for (size_t i=dims.size(); i>1; i--){
		lengths[i-1] = dims[i-2];
		strides[i-1] = dims[i-1]*strides[i];
		if(d < i-1) start *= dims[i-1];
	}
	lengths[d+1] = 1; // now correct for slice dimension, d+1 now we include records
	lengths[0] = seqLen;
	strides[0] = arrLen;


	gslice gsl(start, lengths, strides);
 		
	this->Idata[gsl] = val;

	return;
}

void DvObject::setSliceDim(size_t d, size_t at, std::slice_array <DvString> val)
{
	// 
	
	if(d >= dims.size()){
		this->error("[slice dimension] requested dim > rank");
		return;
	}
	if(at >= dims[d]){
		this->error("[slice dimension] requested row outside dimension");
		return ;
	}
	
	valarray <size_t> lengths(dims.size()+1); // +1 to include record steps
	valarray <size_t> strides(dims.size()+1);
	size_t arrLen = arraySize();
	
	size_t start = at;
	if(dims.size() > 0){
		lengths[dims.size()] = dims[dims.size()-1];
		strides[dims.size()] = 1;
	}
	for (size_t i=dims.size(); i>1; i--){
		lengths[i-1] = dims[i-2];
		strides[i-1] = dims[i-1]*strides[i];
		if(d < i-1) start *= dims[i-1];
	}
	lengths[d+1] = 1; // now correct for slice dimension, d+1 now we include records
	lengths[0] = seqLen;
	strides[0] = arrLen;


	gslice gsl(start, lengths, strides);
 		
	this->Sdata[gsl] = val;

	return;
}

void DvObject::setSliceDim(size_t d, size_t at, std::slice_array <DvTime> val)
{
	// 
	
	if(d >= dims.size()){
		this->error("[slice dimension] requested dim > rank");
		return;
	}
	if(at >= dims[d]){
		this->error("[slice dimension] requested row outside dimension");
		return ;
	}
	
	valarray <size_t> lengths(dims.size()+1); // +1 to include record steps
	valarray <size_t> strides(dims.size()+1);
	size_t arrLen = arraySize();
	
	size_t start = at;
	if(dims.size() > 0){
		lengths[dims.size()] = dims[dims.size()-1];
		strides[dims.size()] = 1;
	}
	for (size_t i=dims.size(); i>1; i--){
		lengths[i-1] = dims[i-2];
		strides[i-1] = dims[i-1]*strides[i];
		if(d < i-1) start *= dims[i-1];
	}
	lengths[d+1] = 1; // now correct for slice dimension, d+1 now we include records
	lengths[0] = seqLen;
	strides[0] = arrLen;


	gslice gsl(start, lengths, strides);
 		
	this->Tdata[gsl] = val;

	return;
}

void DvObject::setSliceDim(size_t d, size_t at, std::slice_array <DvEvent> val)
{
	// 
	
	if(d >= dims.size()){
		this->error("[slice dimension] requested dim > rank");
		return;
	}
	if(at >= dims[d]){
		this->error("[slice dimension] requested row outside dimension");
		return ;
	}
	
	valarray <size_t> lengths(dims.size()+1); // +1 to include record steps
	valarray <size_t> strides(dims.size()+1);
	size_t arrLen = arraySize();
	
	size_t start = at;
	if(dims.size() > 0){
		lengths[dims.size()] = dims[dims.size()-1];
		strides[dims.size()] = 1;
	}
	for (size_t i=dims.size(); i>1; i--){
		lengths[i-1] = dims[i-2];
		strides[i-1] = dims[i-1]*strides[i];
		if(d < i-1) start *= dims[i-1];
	}
	lengths[d+1] = 1; // now correct for slice dimension, d+1 now we include records
	lengths[0] = seqLen;
	strides[0] = arrLen;


	gslice gsl(start, lengths, strides);
 		
	this->Edata[gsl] = val;

	return;
}

// ---------------------

DvObject_var DvObject::sumDim(size_t d, size_t from, size_t to, bool stripBad)
{
    if(is_str() || is_time() || is_event() ){
		DvObject_var ret = new DvObject();
		ret->error("[sum dimension] Not valid for this data type");
		return ret;
	}
	
	if(from == to) return sliceDim(d, from);
	
	if(d >= dims.size()){
		DvObject_var ret = new DvObject();
		ret->error("[sum dimension] requested dim > rank");
		return ret;
	}

	if( to >= dims[d] || from >= dims[d] ){
		DvObject_var ret = new DvObject();
		ret->error("[sum dimension] requested range outside dimension");
		return ret;
	}

	vector <size_t> newDims = dims;
	newDims.erase(newDims.begin()+d);
	
	DvObject_var result = this->create(newDims, seqLen); 
	size_t newArrSize = result->arraySize();
	
	valarray <size_t> lengths(dims.size()+1); // +1 to include record steps
	valarray <size_t> strides(dims.size()+1);
	size_t arrLen = arraySize();
	
	lengths[dims.size()] = dims[dims.size()-1];
	strides[dims.size()] = 1;
	for (size_t i=dims.size()-1; i>0; i--){
		lengths[i] = dims[i-1];
		strides[i] = dims[i]*strides[i+1];
	}
	lengths[d+1] = 1; // now correct for subset dimension, d+1 now we include records
	lengths[0] = seqLen;
	strides[0] = arrLen;

	size_t start_new = 0;
	valarray <size_t> lengths_new(newDims.size()+1);
	valarray <size_t> strides_new(newDims.size()+1);
	if(newDims.size() > 0){
		lengths_new[newDims.size()] = newDims[newDims.size()-1];
		strides_new[newDims.size()] = 1;
	}
	
	for (size_t i=newDims.size(); i>1; i--){
		lengths_new[i-1] = newDims[i-2];
		strides_new[i-1] = newDims[i-1]*strides_new[i];
	}
	lengths_new[0] = seqLen;
	strides_new[0] = newArrSize;

	gslice gsl_new(start_new, lengths_new, strides_new);
	// code here is duplicated for speed under different options
	if(stripBad){
		if(is_dbl()){
			double dFILL = get_dFILL(); 
			valarray <double> val_new(0.0, newArrSize*seqLen) ;
            
            if(to >= from){
                // normal case
                for(size_t k=from; k<=to; k++){
                    size_t start = k;
                    for (size_t i=dims.size()-1; i>d; i--){
                        start *= dims[i];
                    }
                    gslice gsl(start, lengths, strides);
                    valarray <double> val = Ddata[gsl];
                    for(size_t j=0; j<val.size(); j++) {
                        if( !isnan(val[j]) && val[j] != dFILL ) val_new[j] += val[j]; 
                    }
                }
            }
            else{
                // cyclic sum
                for(size_t k=from; k<dims[d]; k++){
                    size_t start = k;
                    for (size_t i=dims.size()-1; i>d; i--){
                        start *= dims[i];
                    }
                    gslice gsl(start, lengths, strides);
                    valarray <double> val = Ddata[gsl];
                    for(size_t j=0; j<val.size(); j++) {
                        if( !isnan(val[j]) && val[j] != dFILL ) val_new[j] += val[j]; 
                    }
                }
                for(size_t k=0; k<=to; k++){
                    size_t start = k;
                    for (size_t i=dims.size()-1; i>d; i--){
                        start *= dims[i];
                    }
                    gslice gsl(start, lengths, strides);
                    valarray <double> val = Ddata[gsl];
                    for(size_t j=0; j<val.size(); j++) {
                        if( !isnan(val[j]) && val[j] != dFILL ) val_new[j] += val[j]; 
                    }
                }
            }
			result->Ddata[gsl_new] = val_new;
		}
		else if(is_int()){
			int iFILL = get_iFILL();
			valarray <int> val_new(0, newArrSize*seqLen) ;

            if(to >= from){
                for(size_t k=from; k<=to; k++){
                    size_t start = k;
                    for (size_t i=dims.size()-1; i>d; i--){
                        start *= dims[i];
                    }
                    gslice gsl(start, lengths, strides);
                    valarray <int> val = Idata[gsl];
                    for(size_t j=0; j<val.size(); j++) {
                        if(val[j] != iFILL) val_new[j] += val[j];
                    }
                }
            }
            else{
                // cyclic sum
                for(size_t k=from; k<dims[d]; k++){
                    size_t start = k;
                    for (size_t i=dims.size()-1; i>d; i--){
                        start *= dims[i];
                    }
                    gslice gsl(start, lengths, strides);
                    valarray <int> val = Idata[gsl];
                    for(size_t j=0; j<val.size(); j++) {
                        if(val[j] != iFILL) val_new[j] += val[j];
                    }
                }
                for(size_t k=0; k<=to; k++){
                    size_t start = k;
                    for (size_t i=dims.size()-1; i>d; i--){
                        start *= dims[i];
                    }
                    gslice gsl(start, lengths, strides);
                    valarray <int> val = Idata[gsl];
                    for(size_t j=0; j<val.size(); j++) {
                        if(val[j] != iFILL) val_new[j] += val[j];
                    }
                }
            }
			result->Idata[gsl_new] = val_new;
		}
	}
	else{
		if(is_dbl()){
            if(to >= from){
                for(size_t k=from; k<=to; k++){
                    size_t start = k;
                    for (size_t i=dims.size()-1; i>d; i--){
                        start *= dims[i];
                    }
				
                    gslice gsl(start, lengths, strides);
				
                    result->Ddata += valarray <double> (Ddata[gsl]);
                }
            }
            else{
                // cyclic sum
                for(size_t k=from; k<dims[d]; k++){
                    size_t start = k;
                    for (size_t i=dims.size()-1; i>d; i--){
                        start *= dims[i];
                    }
                    
                    gslice gsl(start, lengths, strides);
                    
                    result->Ddata += valarray <double> (Ddata[gsl]);
                }
                for(size_t k=0; k<=to; k++){
                    size_t start = k;
                    for (size_t i=dims.size()-1; i>d; i--){
                        start *= dims[i];
                    }
                    
                    gslice gsl(start, lengths, strides);
                    
                    result->Ddata += valarray <double> (Ddata[gsl]);
                }
            }
		}
		else if(is_int()){
            if(to >= from){
                for(size_t k=from; k<=to; k++){
                    size_t start = k;
                    for (size_t i=dims.size()-1; i>d; i--){
                        start *= dims[i];
                    }
                    gslice gsl(start, lengths, strides);
				
                    result->Idata += valarray <int> (Idata[gsl]);
                }
            }
            else{
                // cyclic sum
                for(size_t k=from; k<dims[d]; k++){
                    size_t start = k;
                    for (size_t i=dims.size()-1; i>d; i--){
                        start *= dims[i];
                    }
                    gslice gsl(start, lengths, strides);
                       
                    result->Idata += valarray <int> (Idata[gsl]);
                }
                for(size_t k=0; k<=to; k++){
                    size_t start = k;
                    for (size_t i=dims.size()-1; i>d; i--){
                        start *= dims[i];
                    }
                    gslice gsl(start, lengths, strides);
                       
                    result->Idata += valarray <int> (Idata[gsl]);
                }
            }
		}
	}
	result->copy_xrefs_from(*this);

	
	for(size_t i=d+1; i<=dims.size(); i++){
		DvString dep_i = "DEPEND_";
		dep_i.append(i);		
		result->delete_xref(dep_i); // safe if xref does not exist
		// get next and move down in i
		DvString dep_ip1 = "DEPEND_";
		dep_ip1.append(i+1);
		if(result->xref_exists(dep_ip1)) {
			DvObject_var dxref = result->get_xref(dep_ip1);
			result->change_xref(dep_i, dxref);
		}
	}
	
	for(size_t i=d+1; i<=dims.size(); i++){
		DvString lab_i = "LABEL_";
		lab_i.append(i);
		result->delete_xref(lab_i); // safe if xref does not exist
		// get next and move down in i
		DvString lab_ip1 = "LABEL_";
		lab_ip1.append(i+1);
		if(result->xref_exists(lab_ip1)) {
			DvObject_var lxref = result->get_xref(lab_ip1);
			result->change_xref(lab_i, lxref);
		}
	}

	for(size_t i=d+1; i<=dims.size(); i++){
		DvString rep_i = "REPRESENTATION_";
		rep_i.append(i);
		result->delete_xref(rep_i); // safe if xref does not exist
		// get next and move down in i
		DvString rep_ip1 = "REPRESENTATION_";
		rep_ip1.append(i+1);
		if(result->xref_exists(rep_ip1)) {
			DvObject_var rxref = result->get_xref(rep_ip1);
			result->change_xref(rep_i, rxref);
		}
	}

	DvString lab_str= "sum " ;
	bool need_lab_i = true;
	DvString sumRange(" (");
	DvString depi = "DEPEND_";
	depi.append(d+1);		
	if( xref_exists(depi) ){
		DvObject_var depObj = get_xref(depi);
		sumRange.append(depObj->asDouble(from));
		sumRange += ",";
		sumRange.append(depObj->asDouble(to));
		sumRange += ")";
		if(depObj->xref_exists(LABLAXIS)){
			lab_str += depObj->getXrefText(LABLAXIS);
			need_lab_i=false;
		}
	}
	
	DvString lab_i = "LABEL_";
	lab_i.append(d+1);
	if(need_lab_i){
		if (xref_exists(lab_i)) lab_str += getXrefText(lab_i);
		else if (xref_exists(LABLAXIS)) lab_str += getXrefText(LABLAXIS);
	}
	lab_str += sumRange;
	result->change_xref(LABLAXIS, lab_str);


	DvString mods_str = "Array summed over index ";
    mods_str.append(d);
    mods_str += " (see LABLAXIS)"; 
	result->change_xref(MODS, mods_str);

	if (xref_exists(UNITS)){
		DvString units_str = getXrefText(UNITS);
		units_str += " (summed)"; 
		result->change_xref(UNITS, units_str);
    }

	// change the object SCALEMAX
	if(xref_exists(SCALEMAX))
	{
		DvObject_var obj_max = get_xref(SCALEMAX);
		obj_max = new DvObject(*obj_max);
		*obj_max *= (double) (to - from + 1);
		result->change_xref(SCALEMAX, obj_max);
	}

	return result;
}

// ---------------------

DvObject_var DvObject::averageDim(size_t d, size_t from, size_t to, bool stripBad)
{
	if(is_str() || is_time() || is_event() ){
		DvObject_var ret = new DvObject();
		ret->error("[ave dimension] Not valid for this data type");
		return ret;
	}
	
	if(from == to) return sliceDim(d, from);
	
	if(d >= dims.size()){
		DvObject *ret = new DvObject();
		ret->error("[ave dimension] requested dim > rank");
		return ret;
	}

	if( to >= dims[d] || from >= dims[d] ){
		DvObject_var ret = new DvObject();
		ret->error("[ave dimension] requested range outside dimension");
		return ret;
	}
	
	vector <size_t> newDims = dims;
	newDims.erase(newDims.begin()+d);
	DvObject_var result = this->create(newDims, seqLen); 
	
	size_t newArrSize = result->arraySize();

	valarray <size_t> lengths(dims.size()+1); // +1 to include record steps
	valarray <size_t> strides(dims.size()+1);
	size_t arrLen = arraySize();
	
	lengths[dims.size()] = dims[dims.size()-1];
	strides[dims.size()] = 1;
	for (size_t i=dims.size()-1; i>0; i--){
		lengths[i] = dims[i-1];
		strides[i] = dims[i]*strides[i+1];
	}
	lengths[d+1] = 1; // now correct for subset dimension, d+1 now we include records
	lengths[0] = seqLen;
	strides[0] = arrLen;

	size_t start_new = 0;
	valarray <size_t> lengths_new(newDims.size()+1);
	valarray <size_t> strides_new(newDims.size()+1);
	if(newDims.size() > 0){
		lengths_new[newDims.size()] = newDims[newDims.size()-1];
		strides_new[newDims.size()] = 1;
	}
	for (size_t i=newDims.size(); i>1; i--){
		lengths_new[i-1] = newDims[i-2];
		strides_new[i-1] = newDims[i-1]*strides_new[i];
	}
	lengths_new[0] = seqLen;
	strides_new[0] = newArrSize;

	gslice gsl_new(start_new, lengths_new, strides_new);
	valarray <int> binCounts(0, newArrSize*seqLen);
	// code here is duplicated for speed under different options
	if(stripBad){
		if(is_dbl()){
			double dFILL = get_dFILL(); 
			valarray <double> val_new(0.0, newArrSize*seqLen);
            if(to >= from){
                for(size_t k=from; k<=to; k++){
                    size_t start = k;
                    for (size_t i=dims.size()-1; i>d; i--){
                        start *= dims[i];
                    }
                    gslice gsl(start, lengths, strides);
                    valarray <double> val = Ddata[gsl];
                    for(size_t j=0; j<val.size(); j++) {
                        if( !isnan(val[j]) && val[j] != dFILL ) {
                            val_new[j] += val[j];
                            binCounts[j] +=1;
                        }
                    }
                }
            }
            else{
                // cyclic average
                for(size_t k=from; k<dims[d]; k++){
                    size_t start = k;
                    for (size_t i=dims.size()-1; i>d; i--){
                        start *= dims[i];
                    }
                    gslice gsl(start, lengths, strides);
                    valarray <double> val = Ddata[gsl];
                    for(size_t j=0; j<val.size(); j++) {
                        if( !isnan(val[j]) && val[j] != dFILL ) {
                            val_new[j] += val[j];
                            binCounts[j] +=1;
                        }
                    }
                }
                for(size_t k=0; k<=to; k++){
                    size_t start = k;
                    for (size_t i=dims.size()-1; i>d; i--){
                        start *= dims[i];
                    }
                    gslice gsl(start, lengths, strides);
                    valarray <double> val = Ddata[gsl];
                    for(size_t j=0; j<val.size(); j++) {
                        if( !isnan(val[j]) && val[j] != dFILL ) {
                            val_new[j] += val[j];
                            binCounts[j] +=1;
                        }
                    }
                }
            }
            
            for(size_t j=0; j<val_new.size(); j++){
                if(binCounts[j] > 0) val_new[j] /= binCounts[j];
                else val_new[j] = 0;
            }
            
			result->Ddata[gsl_new] = val_new;

		}
		else if(is_int()){
			int iFILL = get_iFILL();
			valarray <int> val_new(0, newArrSize*seqLen) ;
		
            if(to >= from){
                for(size_t k=from; k<=to; k++){
                    size_t start = k;
                    for (size_t i=dims.size()-1; i>d; i--){
                        start *= dims[i];
                    }
                    gslice gsl(start, lengths, strides);
                    valarray <int> val = Idata[gsl];
                    for(size_t j=0; j<val.size(); j++) {
                        if(val[j] != iFILL)  {
                            val_new[j] += val[j];
                            binCounts[j] +=1;
                        }
                    }
                }
            }
            else{
                // cyclic average
                for(size_t k=from; k<dims[d]; k++){
                    size_t start = k;
                    for (size_t i=dims.size()-1; i>d; i--){
                        start *= dims[i];
                    }
                    gslice gsl(start, lengths, strides);
                    valarray <int> val = Idata[gsl];
                    for(size_t j=0; j<val.size(); j++) {
                        if(val[j] != iFILL)  {
                            val_new[j] += val[j];
                            binCounts[j] +=1;
                        }
                    }
                }
                for(size_t k=0; k<=to; k++){
                    size_t start = k;
                    for (size_t i=dims.size()-1; i>d; i--){
                        start *= dims[i];
                    }
                    gslice gsl(start, lengths, strides);
                    valarray <int> val = Idata[gsl];
                    for(size_t j=0; j<val.size(); j++) {
                        if(val[j] != iFILL)  {
                            val_new[j] += val[j];
                            binCounts[j] +=1;
                        }
                    }
                }
            }
			for(size_t j=0; j<val_new.size(); j++){
				if(binCounts[j] > 0) val_new[j] /= binCounts[j];
				else val_new[j] = 0;
			}
			result->Idata[gsl_new];
		}
	}
	else{
		if(is_dbl()){
           
            if(to >= from){
                for(size_t k=from; k<=to; k++){
                    size_t start = k;
                    for (size_t i=dims.size()-1; i>d; i--){
                        start *= dims[i];
                    }
                    gslice gsl(start, lengths, strides);
				
                    result->Ddata += valarray <double>(Ddata[gsl]);
                }
                result->Ddata /= (to - from + 1);
            }
            else{
                for(size_t k=from; k<dims[d]; k++){
                    size_t start = k;
                    for (size_t i=dims.size()-1; i>d; i--){
                        start *= dims[i];
                    }
                    gslice gsl(start, lengths, strides);
                    
                    result->Ddata += valarray <double>(Ddata[gsl]);
                }
                for(size_t k=0; k<=to; k++){
                    size_t start = k;
                    for (size_t i=dims.size()-1; i>d; i--){
                        start *= dims[i];
                    }
                    gslice gsl(start, lengths, strides);
                    
                    result->Ddata += valarray <double>(Ddata[gsl]);
                }
                result->Ddata /= (to + 1 + dims[d] - from);
                // cyclic average
            }
            
		}
		else if(is_int()){
            if(to >= from){
                for(size_t k=from; k<=to; k++){
                    size_t start = k;
                    for (size_t i=dims.size()-1; i>d; i--){
                        start *= dims[i];
                    }
                    gslice gsl(start, lengths, strides);
                    result->Idata += valarray <int> (Idata[gsl]);
                }
                result->Idata /= (to - from + 1);
            }
            else{
                // cyclic average
                for(size_t k=from; k<dims[d]; k++){
                    size_t start = k;
                    for (size_t i=dims.size()-1; i>d; i--){
                        start *= dims[i];
                    }
                    gslice gsl(start, lengths, strides);
                    result->Idata += valarray <int> (Idata[gsl]);
                }
                for(size_t k=0; k<=to; k++){
                    size_t start = k;
                    for (size_t i=dims.size()-1; i>d; i--){
                        start *= dims[i];
                    }
                    gslice gsl(start, lengths, strides);
                    result->Idata += valarray <int> (Idata[gsl]);
                }
                result->Idata /= (to + 1 + dims[d] - from);
            }
		}
	}

	
	result->copy_xrefs_from(*this);

	
	for(size_t i=d+1; i<=dims.size(); i++){
		DvString dep_i = "DEPEND_";
		dep_i.append(i);		
		result->delete_xref(dep_i); // safe if xref does not exist
		// get next and move down in i
		DvString dep_ip1 = "DEPEND_";
		dep_ip1.append(i+1);
		if(result->xref_exists(dep_ip1)) {
			DvObject_var dxref = result->get_xref(dep_ip1);
			result->change_xref(dep_i, dxref);
		}
	}
	
	for(size_t i=d+1; i<=dims.size(); i++){
		DvString lab_i = "LABEL_";
		lab_i.append(i);
		result->delete_xref(lab_i); // safe if xref does not exist
		// get next and move down in i
		DvString lab_ip1 = "LABEL_";
		lab_ip1.append(i+1);
		if(result->xref_exists(lab_ip1)) {
			DvObject_var lxref = result->get_xref(lab_ip1);
			result->change_xref(lab_i, lxref);
		}
	}

	for(size_t i=d+1; i<=dims.size(); i++){
		DvString rep_i = "REPRESENTATION_";
		rep_i.append(i);
		result->delete_xref(rep_i); // safe if xref does not exist
		// get next and move down in i
		DvString rep_ip1 = "REPRESENTATION_";
		rep_ip1.append(i+1);
		if(result->xref_exists(rep_ip1)) {
			DvObject_var rxref = result->get_xref(rep_ip1);
			result->change_xref(rep_i, rxref);
		}
	}

	DvString lab_str= "ave " ;
	bool need_lab_i = true;
	DvString aveRange(" (");
	DvString depi = "DEPEND_";
	depi.append(d+1);		
	if( xref_exists(depi) ){
		DvObject_var depObj = get_xref(depi);
		aveRange.append(depObj->asDouble(from));
		aveRange += ",";
		aveRange.append(depObj->asDouble(to));
		aveRange += ")";
		if(depObj->xref_exists(LABLAXIS)){
			lab_str += depObj->getXrefText(LABLAXIS);
			need_lab_i=false;
		}
	}
	
	DvString lab_i = "LABEL_";
	lab_i.append(d+1);
	if(need_lab_i){
		if (xref_exists(lab_i)) lab_str += getXrefText(lab_i);
		else if (xref_exists(LABLAXIS)) lab_str += getXrefText(LABLAXIS);
	}
	lab_str += aveRange;
	result->change_xref(LABLAXIS, lab_str);


	DvString mods_str = "Array averaged over index ";
    mods_str.append(d);
    mods_str += " (see LABLAXIS)"; 
	result->change_xref(MODS, mods_str);

	if (xref_exists(UNITS)){
		DvString units_str = getXrefText(UNITS);
		units_str += " (averaged)"; 
		result->change_xref(UNITS, units_str);
    }


	return result;
}

// ---------------------

DvObject_var DvObject::unwrapDim(size_t d, size_t from, size_t to)
{	
	if(from == to) return sliceDim(d, from);
	
	if(d >= dims.size()){
		DvObject_var ret = new DvObject();
		ret->error("[unwrap dimension] requested dim > rank");
		return ret;
	}

	if(to >= dims[d] || from >= dims[d] ){
		DvObject_var ret = new DvObject();
		ret->error("[unwrap dimension] requested range outside dimension");
		return ret;
	}
	
	size_t newDimLen = 1; // changed later
	
	vector <size_t> newDims = dims;
	newDims.erase(newDims.begin()+d);

	DvObject_var result;	
	
	if(to > from){ 
		newDimLen = to - from + 1;
		size_t seqLen_new = seqLen * newDimLen;
	
		result = this->create(newDims, seqLen_new); 
		
		valarray <size_t> lengths(dims.size()+1); // +1 to include record steps
		valarray <size_t> strides(dims.size()+1);
		size_t arrLen = arraySize();

		valarray <size_t> lengths_new(newDims.size()+1); // +1 to include record steps
		valarray <size_t> strides_new(newDims.size()+1);
		size_t arrLen_new = result->arraySize();

		if(dims.size() > 0){
			lengths[dims.size()] = dims[dims.size()-1];
			strides[dims.size()] = 1;
		}
		
		size_t subArrLen = 1;
		for (size_t i=dims.size(); i>1; i--){
			lengths[i-1] = dims[i-2];
			strides[i-1] = dims[i-1]*strides[i];
			if(d < i-1) subArrLen *= dims[i-1];
		}
		lengths[d+1] = 1; // now correct for subset dimension, d+1 now we include records
		lengths[0] = seqLen;
		strides[0] = arrLen;

		lengths_new[newDims.size()] = newDims[newDims.size()-1];
		strides_new[newDims.size()] = 1;
		for (size_t i=newDims.size(); i>1; i--){
			lengths_new[i-1] = newDims[i-2];
			strides_new[i-1] = newDims[i-1]*strides_new[i];
		}
		lengths_new[0] = seqLen; // old seq len as loop k does dim d
		strides_new[0] = arrLen_new * newDimLen;
	
		size_t k_new = 0;

		for(size_t k=from; k<=to; k++){
			size_t start = k*subArrLen;
			size_t start_new = k_new*arrLen_new;

			gslice gsl(start, lengths, strides);
			gslice gsl_new(start_new, lengths_new, strides_new);
 		
			if(is_dbl() )result->Ddata[gsl_new] = Ddata[gsl];
			else if(is_int() )result->Idata[gsl_new] = Idata[gsl];
			else if(is_str() )result->Sdata[gsl_new] = Sdata[gsl];
			else if(is_time() )result->Tdata[gsl_new] = Tdata[gsl];
			else if(is_event() )result->Edata[gsl_new] = Edata[gsl];
			k_new++;
		}
	}
	else { 
		// Handles cyclic array slice
		// if from > to, get (from, N)+(0, to)
		newDimLen = dims[d] - from + to + 1;
		size_t seqLen_new = seqLen * newDimLen;
	
		result = this->create(newDims, seqLen_new); 
	
		valarray <size_t> lengths(dims.size()+1); // +1 to include record steps
		valarray <size_t> strides(dims.size()+1);
		size_t arrLen = arraySize();

		valarray <size_t> lengths_new(newDims.size()+1); // +1 to include record steps
		valarray <size_t> strides_new(newDims.size()+1);
		size_t arrLen_new = result->arraySize();
	
		size_t k_new = 0;
		
		if(dims.size()>0){
			lengths[dims.size()] = dims[dims.size()-1];
			strides[dims.size()] = 1;
		}
		size_t subArrLen = 1;
		for (size_t i=dims.size(); i>1; i--){
			lengths[i-1] = dims[i-2];
			strides[i-1] = dims[i-1]*strides[i];
			if(d < i-1) subArrLen *= dims[i-1];
		}
		lengths[d+1] = 1; // now correct for subset dimension, d+1 now we include records
		lengths[0] = seqLen;
		strides[0] = arrLen;

		lengths_new[newDims.size()] = newDims[newDims.size()-1];
		strides_new[newDims.size()] = 1;
		for (size_t i=newDims.size(); i>1; i--){
			lengths_new[i-1] = newDims[i-2];
			strides_new[i-1] = newDims[i-1]*strides_new[i];
		}
		lengths_new[0] = seqLen; // old seq len as loop k does dim d
		strides_new[0] = arrLen_new * newDimLen;
		
		// start with (from, N) range
		for(size_t k=from; k<dims[d]; k++){
			size_t start = k*subArrLen;

			size_t start_new = k_new*arrLen_new;

			gslice gsl(start, lengths, strides);
			gslice gsl_new(start_new, lengths_new, strides_new);
 		
			if(is_dbl() )result->Ddata[gsl_new] = Ddata[gsl];
			else if(is_int() )result->Idata[gsl_new] = Idata[gsl];
			else if(is_str() )result->Sdata[gsl_new] = Sdata[gsl];
			else if(is_time() )result->Tdata[gsl_new] = Tdata[gsl];
			else if(is_event() )result->Edata[gsl_new] = Edata[gsl];
			k_new++;
		}
		
		// now insert (0, to) range
		for(size_t k=0; k<=to; k++){
			size_t start = k*subArrLen;

			size_t start_new = k_new*arrLen_new;

			gslice gsl(start, lengths, strides);
			gslice gsl_new(start_new, lengths_new, strides_new);
 		
			if(is_dbl() )result->Ddata[gsl_new] = Ddata[gsl];
			else if(is_int() )result->Idata[gsl_new] = Idata[gsl];
			else if(is_str() )result->Sdata[gsl_new] = Sdata[gsl];
			else if(is_time() )result->Tdata[gsl_new] = Tdata[gsl];
			else if(is_event() )result->Edata[gsl_new] = Edata[gsl];
			k_new++;
		}

	}
	
	result->copy_xrefs_from(*this);
	
	for(size_t i=d+1; i<=dims.size(); i++){
		DvString dep_i = "DEPEND_";
		dep_i.append(i);		
		result->delete_xref(dep_i); // safe if xref does not exist
		// get next and move down in i
		DvString dep_ip1 = "DEPEND_";
		dep_ip1.append(i+1);
		if(result->xref_exists(dep_ip1)) {
			DvObject_var dxref = result->get_xref(dep_ip1);
			result->change_xref(dep_i, dxref);
		}
	}
	
	for(size_t i=d+1; i<=dims.size(); i++){
		DvString lab_i = "LABEL_";
		lab_i.append(i);
		result->delete_xref(lab_i); // safe if xref does not exist
		// get next and move down in i
		DvString lab_ip1 = "LABEL_";
		lab_ip1.append(i+1);
		if(result->xref_exists(lab_ip1)) {
			DvObject_var lxref = result->get_xref(lab_ip1);
			result->change_xref(lab_i, lxref);
		}
	}

	for(size_t i=d+1; i<=dims.size(); i++){
		DvString rep_i = "REPRESENTATION_";
		rep_i.append(i);
		result->delete_xref(rep_i); // safe if xref does not exist
		// get next and move down in i
		DvString rep_ip1 = "REPRESENTATION_";
		rep_ip1.append(i+1);
		if(result->xref_exists(rep_ip1)) {
			DvObject_var rxref = result->get_xref(rep_ip1);
			result->change_xref(rep_i, rxref);
		}
	}

	DvString lab_str= "unwrapped " ;
	bool need_lab_i = true;
	DvString sumRange(" (");
	DvString depi = "DEPEND_";
	depi.append(d+1);		
	if( xref_exists(depi) ){
		DvObject_var depObj = get_xref(depi);
		sumRange.append(depObj->asDouble(from));
		sumRange += ",";
		sumRange.append(depObj->asDouble(to));
		sumRange += ")";
		if(depObj->xref_exists(LABLAXIS)){
			lab_str += depObj->getXrefText(LABLAXIS);
			need_lab_i=false;
		}
	}
	
	DvString lab_i = "LABEL_";
	lab_i.append(d+1);
	if(need_lab_i){
		if (xref_exists(lab_i)) lab_str += getXrefText(lab_i);
		else if (xref_exists(LABLAXIS)) lab_str += getXrefText(LABLAXIS);
	}
	lab_str += sumRange;
	result->change_xref(LABLAXIS, lab_str);


	DvString mods_str = "Array unwrapped (aliased to time) over index ";
    mods_str.append(d);
    mods_str += " (see LABLAXIS)"; 
	result->change_xref(MODS, mods_str);

	if (xref_exists(UNITS)){
		DvString units_str = getXrefText(UNITS);
		units_str += " (unwrapped)"; 
		result->change_xref(UNITS, units_str);
    }
	
	result->unwrapXrefs(seqLen, newDimLen);

	return result;
}

// -----------------------------------------

void DvObject::unwrapXrefs(size_t oldSeqLen,  size_t nElem){
	// THIS is ONLY USED as a utility call by unwrapDim
	
	// Do DEPEND_0 differently as data spread along sequence
	DvObject_var D0 = getDep0();

	if( D0.is_ok() && D0->seqSize() == oldSeqLen) { 

		DvObject_var D0_out = D0->create(D0->Dims(), oldSeqLen*nElem);
			
			DvObject_var DP = D0->get_xref(DELTA_PLUS);
			DvObject_var DM = D0->get_xref(DELTA_MINUS);
			
			DvObject_var Dstart;
			DvObject_var Dend;
			
			if(DP.is_ok()) Dend = D0 + DP;
			else Dend = D0;
			
			if(DM.is_ok()) Dstart = D0 - DM;
			else Dstart = D0;
			
			D0_out->copy_xrefs_from(D0);
			if( D0_out->xref_exists(DEPEND_0) )	D0_out->delete_xref(DEPEND_0);
			
			double delta=1.0;
			bool calcDelta = true;
			
			if(Dstart == Dend) {
				calcDelta = false;
				delta = D0->get_spacing() / nElem;
			}
			
			for(size_t i=0; i<oldSeqLen; i++){
				if(calcDelta) delta = (Dend->asDouble(i) - Dstart->asDouble(i)) / nElem;
				for(size_t k=0; k<nElem; k++){
					// smear along axis
					D0_out[i*nElem + k] = Dstart[i];
					D0_out[i*nElem + k] += (double) (k*delta);
				}
			}			
			
			change_xref(DEPEND_0, D0_out);
			D0_out->unwrapXrefs(oldSeqLen, nElem);
			
			// now subdivide delta plus/minus for Depend 0
			D0_out->adjustDeltas(nElem); // if a sequence it had records added already
	}
	
	// repeat entries along other sequence type xrefs
  	DvNode *xref = first_xref();
	while(xref){	
		
		if(xref->name() != DEPEND_0) {  
		
			DvObject_var xref_i = xref->obj();
	
			if( xref_i.is_ok() && oldSeqLen == xref_i->seqSize()) {
			  
				DvObject_var xref_i_out = xref_i->create(xref_i->Dims(), oldSeqLen*nElem);
				for(size_t j=0; j<oldSeqLen; j++){
					for(size_t k=0; k<nElem; k++){
						// duplicate along axis
						xref_i_out[j*nElem + k] = xref_i[j];
					}
				}			
				xref_i_out->copy_xrefs_from(xref_i);
				if( xref_i_out->xref_exists(DEPEND_0) ) xref_i_out->delete_xref(DEPEND_0);
		
				xref->setObj(xref_i_out);
				xref_i_out->unwrapXrefs(oldSeqLen, nElem);
			}
		}
		xref = xref->next;
	}
		
	
	return;
}	  

void DvObject::adjustDeltas(size_t nElem)
{
	// USED only as support utility for unwrapDim, and called from unwrapXrefs
	// divide deltas by number of elements unwrapped along time axis	
	
	DvObject_var dp = get_xref(DELTA_PLUS);
	if(dp.is_ok() ){
			DvObject_var dp_new = new DvObject(dp);
			*dp_new  /= (double) nElem;
			change_xref(DELTA_PLUS, dp_new);
	}
	
	DvObject_var dm = get_xref(DELTA_MINUS);
	if(dm.is_ok() ){
			DvObject_var dm_new = new DvObject(dm);
			*dm_new  /= (double) nElem;
			change_xref(DELTA_MINUS, dm_new);
	}
	
}

// ---------------------

DvObject_var DvObject::reverseDim(size_t d, size_t from, size_t to)
{	
	if(from == to) return sliceDim(d, from);
	
	if(d >= dims.size()){
		DvObject_var ret = new DvObject();
		ret->error("[reverse dimension] requested dim > rank");
		return ret;
	}

	if(to < from){
		DvObject_var ret = new DvObject();
		ret->error("[reverse dimension] requested start after end");
		return ret;
	}

	if(to >= dims[d] ||  from >= dims[d] ){
		DvObject_var ret = new DvObject();
		ret->error("[reverse dimension] requested range outside dimension");
		return ret;
	}
	
	size_t newDimLen = to - from + 1; 
	
	vector <size_t> newDims = dims;
	newDims[d] = newDimLen;

	DvObject_var result = this->create(newDims, seqLen); 
		
	valarray <size_t> lengths(dims.size()+1); // +1 to include record steps
	valarray <size_t> strides(dims.size()+1);
	size_t arrLen = arraySize();

	valarray <size_t> lengths_new(newDims.size()+1); // +1 to include record steps
	valarray <size_t> strides_new(newDims.size()+1);
	size_t arrLen_new = result->arraySize();
	
	if(dims.size() > 0){
		lengths[dims.size()] = dims[dims.size()-1];
		strides[dims.size()] = 1;
	}
	for (int i=dims.size(); i>1; i--){
		lengths[i-1] = dims[i-2];
		strides[i-1] = dims[i-1]*strides[i];
	}
	lengths[d+1] = 1; // now correct for subset dimension, d+1 now we include records
	lengths[0] = seqLen;
	strides[0] = arrLen;

	lengths_new[newDims.size()] = newDims[newDims.size()-1];
	strides_new[newDims.size()] = 1;
	for (int i=newDims.size(); i>1; i--){
		lengths_new[i-1] = newDims[i-2];
		strides_new[i-1] = newDims[i-1]*strides_new[i];
	}
	lengths_new[d+1] = 1; // now correct for subset dimension, d+1 now we include records
	lengths_new[0] = seqLen; 
	strides_new[0] = arrLen_new;
	

	size_t k_new = newDimLen - 1;
	for(size_t k=from; k<=to; k++){
	
		size_t start = k;
		size_t start_new = k_new;
		for (size_t i=d+1; i<dims.size(); i++){
			start *= dims[i];
			start_new *= newDims[i];
		}


		gslice gsl(start, lengths, strides);
		gslice gsl_new(start_new, lengths_new, strides_new);

		if(is_dbl() )result->Ddata[gsl_new] = Ddata[gsl];
		else if(is_int() )result->Idata[gsl_new] = Idata[gsl];
		else if(is_str() )result->Sdata[gsl_new] = Sdata[gsl];
		else if(is_time() )result->Tdata[gsl_new] = Tdata[gsl];
		else if(is_event() )result->Edata[gsl_new] = Edata[gsl];
		
		k_new--;
	}
	
	result->copy_xrefs_from(*this);
	
	DvString dep_i = "DEPEND_";
	dep_i.append(d+1);		
	if(result->xref_exists(dep_i)) {
		DvObject_var dxref = result->get_xref(dep_i);
		DvObject_var dxref_new = dxref->reverseDim(0, from, to);
		result->change_xref(dep_i, dxref_new);
	}
	
	DvString lab_i = "LABEL_";
	lab_i.append(d+1);
	if(result->xref_exists(lab_i)) {
		DvObject_var lxref = result->get_xref(lab_i);
		DvObject_var lxref_new = lxref->reverseDim(0, from, to);
		result->change_xref(lab_i, lxref_new);
	}

	DvString rep_i = "REPRESENTATION_";
	rep_i.append(d+1);
	if(result->xref_exists(rep_i)) {
		DvObject_var rxref = result->get_xref(rep_i);
		DvObject_var rxref_new = rxref->reverseDim(0, from, to);
		result->change_xref(rep_i, rxref_new);
	}

	DvString lab_str= "order reversed " ;
	bool need_lab_i = true;
	DvString sumRange(" (");
	DvString depi = "DEPEND_";
	depi.append(d+1);		
	if( xref_exists(depi) ){
		DvObject_var depObj = get_xref(depi);
		sumRange.append(depObj->asDouble(from));
		sumRange += ",";
		sumRange.append(depObj->asDouble(to));
		sumRange += ")";
		if(depObj->xref_exists(LABLAXIS)){
			lab_str += depObj->getXrefText(LABLAXIS);
			need_lab_i=false;
		}
	}
	
	if(need_lab_i){
		if (xref_exists(lab_i)) lab_str += getXrefText(lab_i);
		else if (xref_exists(LABLAXIS)) lab_str += getXrefText(LABLAXIS);
	}
	lab_str += sumRange;
	result->change_xref(LABLAXIS, lab_str);


	DvString mods_str = "Array reversed over index ";
    mods_str.append(d);
    mods_str += " (see LABLAXIS)"; 
	result->change_xref(MODS, mods_str);

	
	return result;
}

DvString DvObject::getCommonText(){

	if(this->is_nil()) return DvString("");
	
	DvString txt = this->asStr(0);
	for(size_t i=1; i<totalSize(); i++){
		DvString txt2 = this->asStr(i);
		txt = txt.common(txt2);
	}
	return txt;
}

void DvObject::cycleMod(int nstart, int n, int modn){
	// can also cycle scalar sequence
	// move modn elements cyclicly by n starting at nstart
	
	
	if(is_dbl()){
		valarray<double>copy(0.0, modn);
		copy = Ddata[slice(nstart, modn, 1)];
		
		for(int i=0; i< modn; i++){
			int ito = i + n;
			if(ito > modn) ito -= modn;
			if(ito < 0) ito += modn;

			Ddata[ito+nstart] = copy[i];
		}
	}
	else if(is_int()){
		valarray<int>copy(0, modn);
		copy = Idata[slice(nstart, modn, 1)];
		
		for(int i=0; i< modn; i++){
			int ito = i + n;
			if(ito > modn) ito -= modn;
			if(ito < 0) ito += modn;

			Idata[ito+nstart] = copy[i];
		}
	}
	else if(is_str()){
		valarray<DvString>copy(modn);
		copy = Sdata[slice(nstart, modn, 1)];
		
		for(int i=0; i< modn; i++){
			int ito = i + n;
			if(ito > modn) ito -= modn;
			if(ito < 0) ito += modn;

			Sdata[ito+nstart] = copy[i];
		}
	}
	else if(is_time()){
		valarray<DvTime>copy(modn);
		copy = Tdata[slice(nstart, modn, 1)];
		
		for(int i=0; i< modn; i++){
			int ito = i + n;
			if(ito > modn) ito -= modn;
			if(ito < 0) ito += modn;

			Tdata[ito+nstart] = copy[i];
		}
	}
	else if(is_event()){
		valarray<DvEvent>copy(modn);
		copy = Edata[slice(nstart, modn, 1)];
		
		for(int i=0; i< modn; i++){
			int ito = i + n;
			if(ito > modn) ito -= modn;
			if(ito < 0) ito += modn;

			Edata[ito+nstart] = copy[i];
		}
	}
	
}


// only needed for outer product below

void outerProdLoopB(DvObject_var &res, size_t &i_A, size_t &j_B, size_t &ij, DvObject *objA, DvObject *objB, vector <size_t> B_dims){
    
    if(B_dims.size() > 0){
        vector <size_t> L_dims;
        for(size_t m=1; m< B_dims.size(); m++) L_dims.push_back(B_dims[m]);
        
        for(size_t i=0; i<B_dims[0]; i++){
            outerProdLoopB(res, i_A, j_B, ij, objA, objB, L_dims);
        }
    }
    else {
        if(objB->seqSize() > 1){
            for (size_t n=0; n<objA->seqSize(); n++){
                res->dbl(n,ij) = objA->dbl(n,i_A) * objB->dbl(n,j_B);
            }
        }
        else{
            for (size_t n=0; n<objA->seqSize(); n++){
                res->dbl(n,ij) = objA->dbl(n,i_A) * objB->dbl(0,j_B);
            }
        }
        
        j_B++;
        ij++;
    }
    return;
}

void outerProdLoopA(DvObject_var &res, size_t &i_A, size_t &ij, DvObject *objA, DvObject *objB, vector <size_t> A_dims){
    if(A_dims.size() > 0){
        vector <size_t> L_dims;
        for(size_t m=1; m< A_dims.size(); m++) L_dims.push_back(A_dims[m]);
        
        for(size_t i=0; i<A_dims[0]; i++){
            outerProdLoopA(res, i_A, ij, objA, objB, L_dims);
        }
    }
    else {
        size_t i_B=0;
        vector <size_t> B_dims;
        for(size_t m=1; m<objB->nDims(); m++) B_dims.push_back(objB->Dims()[m]);
        
        for(size_t j=0; j<objB->Dims()[0]; j++){
            outerProdLoopB(res, i_A, i_B, ij, objA, objB, B_dims);
        }
        i_A++;
    }
    return;
}

DvObject_var DvObject::outerProduct(DvObject_var &obj){
    
    if( !sameFrame(*obj) ) {
        DvObject_var res = new DvObject();
        res->error("[outer product] frames differ");
        return res;
    }
    
    if(this->seqSize() == 1 && obj->seqSize() > 1){
        // make 1st object into sequence
        DvObject_var thisSeq = new DvObject(this, obj->seqSize());
        DvObject_var tt = obj->getDep0();
        if(tt.is_ok()) thisSeq->change_xref(DEPEND_0, tt);
        DvObject_var res = thisSeq->outerProduct(obj);
        
        return res;
    }
        
    //  Rij= Ai * Bj (i= i1,i2,...   j=j1,j2,...)
                            
    if( is_dbl() && obj->is_dbl() ){
        size_t len = seqSize();
        // set up new dimension
        vector <size_t> newDims;
        for(size_t i=0; i< this->nDims(); i++) newDims.push_back(this->dims[i]);
        for(size_t j=0; j< obj->nDims(); j++) newDims.push_back(obj->Dims()[j]);
       
        DvObject_var res = new DvObject(0.0, newDims, len); 
        res->copy_xrefs_from(*this);
                    
            
        DvObject_var arg = obj;
        if(obj->seqSize() > 1) arg = obj->Join(*this);  // may be unchanged

        size_t ij=0;
        size_t i_A=0;
        // (A) this dims
        vector <size_t> A_dims;
        for(size_t m=1; m< this->nDims(); m++) A_dims.push_back(this->dims[m]);
                
        for(size_t i=0; i<dims[0]; i++){
            outerProdLoopA(res, i_A, ij, this, arg.ptr(), A_dims);
                                
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
        
        // copy deps from B
        int i = (int)nDims();
        for(int j=1; j<= obj->nDims(); j++){
            DvString depj = "DEPEND_";
            depj += j;
            DvString depij = "DEPEND_";
            depij += i+j;
            DvObject_var dep_xref = obj->get_xref(depj);
            res->change_xref(depij, dep_xref);
            
        }
        
        return res;
    }
    else  {
        DvObject_var res = new DvObject();
        res->error("[outer product] not possible");
        return res;
    }
}


// ---------------------------------------------------------------------


// Eigenvector calls (after: Janet Barnes, Stein Haaland, Eugen Sorballo)

const double kEpsilon = 1E-20;  //used to see if a float variable is zero.

double e_house(valarray <double> &x, int start, int n, valarray<double> &v)
{
	double beta = 0., sigma = 0.;

	for (int i = 1; i < n; i++) {
		v[i] = x[start+i];
		sigma += x[start+i] * x[start+i];
	}
	v[0] = 1.;

	if (sigma != 0.) {
		double mu = std::sqrt(x[start] * x[start] + sigma);
		if (x[start] < 0.) {
			v[0] = x[start] - mu;
		}
		else {
			v[0] = -sigma / (x[start] + mu);
		}
		beta = 2 * v[0] * v[0] / (sigma + v[0] * v[0]);
		for (int i = 1; i < n; i++) {
			v[i] = v[i] / v[0];
		}
		v[0] = 1.;
	}
	return beta;
}

void e_tridiag(valarray<double> &a, valarray<double> &q, int n)
{
	int i, j, k, sub_size;
	double beta, coeff, norm;
	valarray<double> v(0.0, n - 1); // Container for Householder vectors
	valarray<double> p(0.0, n - 1); // Container for Householder matrices
	valarray<double> w(0.0, n - 1);
	valarray<double> beta_arr(0.0, n - 1);
	valarray<double> q_temp(0.0, n * n);

	for (k = 0; k < n; k++) {
		q[k * n + k] = 1.;
	}

	for (k = 0; k < n - 2; k++) {
		sub_size = n - k - 1;

		beta = e_house(a, k*n+k+1, sub_size, v);
		beta_arr[k] = beta;
		//p=b*a(k+1:n, k+1:n)*v;

		for (i = 0; i < sub_size; i++) {
			p[i] = 0.;
			for (j = 0; j < sub_size; j++) {
				p[i] += a[(i + k + 1) * n + j + k + 1] * v[j];
			}
			p[i] *= beta;
		}

		//w=p-(b*p'*v/2)*v;
		coeff = 0.;
		norm = 0.;
		for (i = 0; i < sub_size; i++) {
			coeff += p[i] * v[i];
			norm += a[(i + k + 1) * n + k] * a[(i + k + 1) * n + k];
		}
		coeff = coeff * beta / 2.;
		for (i = 0; i < sub_size; i++) {
			w[i] = p[i] - coeff * v[i];
		}

		//a(k+1, k)=norm(a(k+1:n, k)); a(k, k+1)=a(k+1, k);
		norm = std::sqrt(norm);
		a[(k + 1) * n + k] = a[k * n + k + 1] = norm;

		for (i = 0; i < sub_size; i++) {
			for (j = 0; j <= i; j++) {
				a[(i + k + 1) * n + j + k + 1] = a[(j + k + 1) * n + i + k + 1] = a[(i + k + 1) * n + j + k + 1] - v[i]
						* w[j] - w[i] * v[j];
			}
			// v storage for Q computation
			if (i != 0) {
				a[(i + k + 1) * n + k] = v[i];
			}
		}
	}

	for (k = n - 3; k >= 0; k--) {
		sub_size = n - k - 1;
		for (i = 0; i < sub_size; i++) {
			if (i == 0)
				v[i] = 1;
			else {
				v[i] = a[(i + k + 1) * n + k];
				a[(i + k + 1) * n + k] = a[k * n + i + k + 1] = 0.;
			}
		}
		// Copy of q in q_temp
		for (i = 0; i < sub_size; i++) {
			for (j = 0; j < sub_size; j++) {
				q_temp[(i + k + 1) * n + j + k + 1] = q[(i + k + 1) * n + j + k + 1];
			}
		}
		for (i = 0; i < sub_size; i++) {
			for (j = 0; j < sub_size; j++) {
				coeff = 0.;
				for (int l = 0; l < sub_size; l++) {
					coeff += v[l] * q_temp[(l + k + 1) * n + j + k + 1];
				}
				q[(i + k + 1) * n + j + k + 1] -= beta_arr[k] * coeff * v[i];
			}
		}
	}

}

void e_givens(double a, double b, double * c, double * s)
{
	if (b == 0) {
		*c = 1.;
		*s = 0.;
	}
	else {
		double tau;
		if (fabs(b) > fabs(a)) {
			tau = -a / b;
			*s = 1. / std::sqrt(1 + tau * tau);
			*c = (*s) * tau;
		}
		else {
			tau = -b / a;
			*c = 1. / std::sqrt(1 + tau * tau);
			*s = (*c) * tau;
		}
	}
}

void apply_e_givens(valarray<double> &a, int size, double c, double s, int i, int k, int p, int q)
{
	int j;
	double t1, t2;

	//one side of the diagonal
	for (j = p; j < q; j++) {
		t1 = a[i * size + j];
		t2 = a[k * size + j];
		a[i * size + j] = c * t1 - s * t2;
		a[k * size + j] = s * t1 + c * t2;
	}

	//same procedure to the other side of the diagonal.
	for (j = p; j < q; j++) {
		t1 = a[j * size + i];
		t2 = a[j * size + k];
		a[j * size + i] = c * t1 - s * t2;
		a[j * size + k] = s * t1 + c * t2;
	}
}

void e_symmqr_step(valarray<double> &t, valarray<double> &Q, int size, int p, int q)
{
	double d = (t[(q - 2) * size + q - 2] - t[(q - 1) * size + q - 1]) / 2.;
	double tq2 = t[(q - 1) * size + q - 2] * t[(q - 1) * size + q - 2];
	//we know that off-diagonal elements are not zero, so tq2<>0
	double sg = (d >= 0.0 ? 1.0 : -1.0);
	// d can be zero, so 'sg' should not be zero, otherwise we have Inf.
	double mu = t[(q - 1) * size + q - 1] - tq2 / (d + sg * std::sqrt(d * d + tq2));

	double x = t[p * size + p] - mu;
	double z = t[(p + 1) * size + p];
	int k, j;
	double c, s;
	double t1, t2;
	for (k = p; k < q - 1; k++) {
		e_givens(x, z, &c, &s);
		apply_e_givens(t, size, c, s, k, k + 1, p, q);

		for (j = 0; j < size; j++) {
			t1 = Q[j * size + k];
			t2 = Q[j * size + k + 1];
			Q[j * size + k] = c * t1 - s * t2;
			Q[j * size + k + 1] = s * t1 + c * t2;
		}

		if (k < q - 2) {
			x = t[(k + 1) * size + k];
			z = t[(k + 2) * size + k];
		}
	}
}

void e_display_mat(valarray<double> &t, int size1)
{
	for (int i = 0; i < size1; ++i) {
		for (int j = 0; j < size1; ++j)
			std::cout << t[i * size1 + j] << " ";
		std::cout << std::endl;
	}
}

void e_symmqr(valarray<double> &t, valarray<double> &Q, int size)
{
	int p, q, i, end_loop, iter;
	p = 0;
	q = size - 1;
	end_loop = 0;

	//The hash is used to compute a Check Sum of the matrix.
	//If check sums after two iterations are equal, then the procedure must give up.
	double hash[20];
	hash[0] = 1;
	for (int i = 1; i < 20; i++)
		hash[i] = hash[i - 1] * 10;
	double prev_sum = 0, new_sum = 1;
	int iteration = 500;

	while ((q > p) && (new_sum != prev_sum) && (iteration > 0)) {
		//some code for checking if we don't loop too much
		iteration--;
		prev_sum = new_sum;
		new_sum = 0;
		for (int ii = 0; ii < size; ++ii)
			for (int j = 0; j < size; ++j)
				new_sum += t[ii * size + j] * hash[ii * size + j];

		//if off-diagonal elements are almost zero, make them zero.
		for (i = 0; i < size - 1; i++) {
			if (fabs(t[i * size + i + 1]) <= kEpsilon) {
				t[i * size + i + 1] = t[(i + 1) * size + i] = 0.;
			}
		}

		//Find portions where we have non-diagonal elements:
		end_loop = 0;
		iter = q - 1;  //look from the up
		while (!end_loop && iter > 0) {
			if (fabs(t[(iter + 1) * size + iter]) < kEpsilon) {
				q = iter;
				iter--;
			}
			else {
				end_loop = 1;
			}
		}

		end_loop = 0;
		iter = p;  //look from bottom
		while (!end_loop && iter < size - 1) {
			if (fabs(t[(iter + 1) * size + iter]) < kEpsilon) {
				p = ++iter;
			}
			else {
				end_loop = 1;
			}
		}

		if (q > p) {
			e_symmqr_step(t, Q, size, p, q + 1);
//			e_display_mat(t, size, size);
		}

	}
}

void e_sort(valarray<double> &a, valarray<double> &q, int size)
{
	double temp;
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size - i - 1; j++) {
			if (a[j * size + j] > a[(j + 1) * size + j + 1]) {
				temp = a[(j + 1) * size + j + 1];
				a[(j + 1) * size + j + 1] = a[j * size + j];
				a[j * size + j] = temp;

				for (int k = 0; k < size; k++) {
					temp = q[k * size + j + 1];
					q[k * size + j + 1] = q[k * size + j];
					q[k * size + j] = temp;
				}

			}
		}
	}
}

void d_eigen(valarray<double> &a, valarray<double> &v, size_t n)
{
	// eigenvectors are columns
	
	valarray<double> q(0.0, n*n);

	e_tridiag(a, q, (int) n);


	e_symmqr(a, q, (int) n);


	e_sort(a, q, (int) n);


	for (size_t i = 0; i < n; i++)
		v[i] = a[i * n + i];

		
	a = q; // copy valarray to return

}

void c_eigen(std::vector<std::complex<double> > &array, std::vector<double> &eigvalues, size_t size)
{
	using namespace std;
	if (size * size == array.size())
	{
		//if array = A + iB, then the new array is [[A, -B]; [B, A]]
		valarray<double> equiv_array(0.0, 4*size*size);
		for (size_t i = 0; i < size; ++i)
			for (size_t j = 0; j < size; ++j)
			{
				equiv_array[i * (2*size) + j] = array[i * size + j].real();  //A
				equiv_array[i * (2*size) + j + size] = -array[i * size + j].imag();  //-B
				equiv_array[(i + size) * (2*size) + j] = array[i * size + j].imag();  //B
				equiv_array[(i + size) * (2*size) + j + size] = array[i * size + j].real();  //A
			}
			
		valarray<double> evalues(0.0, 2 * size);
		d_eigen(equiv_array, evalues, 2*size);

		//Compound the solution
		eigvalues.empty();
		eigvalues.resize(size, 0.0);
		for (size_t col = 0; col < size; ++col)
		{
			//we assume here that the matrix was hermitian => eigenvalues are real!!!
			eigvalues[col] = evalues[2 * col];
			for (size_t row = 0; row < size; ++row)
				array[row*size + col] = complex<double>(equiv_array[row * (2*size) + 2*col],
														equiv_array[(row+size) * (2*size) + 2*col]);
		}

	}
}
