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


// ----------------  General translation ---------------- 

DvObject_var DvObject::minusStart(){

	// convert data to value - value of first record
	// Time tags become seconds since start of sequence
	
	
	if(this->is_time()){
		DvObject_var ret = new DvObject(0.0, this->Dims(), this->seqSize()); 
		size_t arrSize = this->arraySize();
		for (size_t n=0; n<this->seqSize(); n++){
			for(size_t i=0; i<arrSize; i++) ret->Ddata[n*arrSize+i] = this->Tdata[n*arrSize+i] - this->Tdata[0];
		}
		return ret;	
	}
	else{
		DvObject_var ret = new DvObject(*this);
		DvRecord rec0 = (*this)[0];
		for (size_t n=0; n<this->seqSize(); n++){
			(*ret)[n] -= rec0;
		}
		return ret;	
	}
		
	return new DvObject();
}

// --------------------- Integration ------------------------------

DvObject_var DvObject::integrate(DvString gapMethod, DvString integralMethod){


	if( !is_dbl() ){
		DvObject_var ret = new DvObject();
		ret->error("[Integral] Not available for this data type");
		return ret;
	}
	
	size_t length = this->seqSize();
	if(length < 2){
		DvObject *ret = new DvObject();
		ret->error("[Integral] Sequence length too small");
		return ret;
	}
	
	double deltaTop = 0.0;
	
	DvObject_var d0 = getDep0();
	if(d0.is_nil() || d0->seqSize() != length) {
		DvObject *ret = new DvObject();
		ret->error("[Integral] No useable DEPEND_0");
		return ret;
	}

	DvObject_var x = new DvObject(0.0, length);
	DvObject_var xm = new DvObject(0.0, length);
	DvObject_var xp = new DvObject(0.0, length);
	
	if(d0->xref_exists(DELTA_PLUS) && d0->xref_exists(DELTA_MINUS) && gapMethod != DV_NN){

		// deltas exist and gaps not NN (so use bin width from deltas)
		DvObject_var dp = d0->get_xref(DELTA_PLUS);
		DvObject_var dm = d0->get_xref(DELTA_MINUS);
		
		if(is_time()) {
			dp = dp->changeUnitsTo("1.0>s", "s");
			dm = dm->changeUnitsTo("1.0>s", "s");
		}
		
		xm->dbl(0) = 0.0;
		x->dbl(0) = dm->dbl(0);
		xp->dbl(0) = x->dbl(0) + dp->dbl(0);
		for(size_t i=1; i<length; i++){
			// done this way to get scalar from time
			x->dbl(i) = d0->asDouble(i) - d0->asDouble(i-1) + x->dbl(i-1);
			xm->dbl(i) = x->dbl(i) - dm->dbl(i);
			xp->dbl(i) = x->dbl(i) + dp->dbl(i);
		}
		deltaTop = dp->dbl(length-1);
			
	}
	else if(gapMethod == DV_NN || gapMethod != DV_REMOVE ){
		// Fake deltas to mid points
		x->dbl(0) = 0.0;
		xm->dbl(0) = 0.0;
		for(size_t i=1; i<length; i++){
			// done this way to get scalar from time
			xm->dbl(i) = x->dbl(i-1) + 0.5*(d0->asDouble(i) - d0->asDouble(i-1));
			x->dbl(i) = x->dbl(i-1) + d0->asDouble(i) - d0->asDouble(i-1);
			xp->dbl(i-1) = xm->dbl(i);
		}
		xp->dbl(length-1) = x->dbl(length-1) + deltaTop;
	
	}
	else{
		DvObject_var ret = new DvObject();
		ret->error("[Integral] Needs Delta_plus/minus for IGNORE GAPS option");
		return ret;
	}
		
	DvObject_var retObj;
	valarray <double> val(arraySize());
		
	// Note step lengths may not be constant, so Trapezium, Simpson's rule etc don't work.
	
	if(gapMethod == DV_NN || gapMethod == DV_REMOVE){
		// histogram
		// same method for ignore gaps or histogram as difference is only in xm and xp bin boundaries
		if(integralMethod == DV_DEF_INTEGRAL){
	
			// definite integral
			retObj = new DvObject(0.0, Dims(), 1);
			for(size_t i=0; i<length; i++){
				val = this->record(i);
				val *= (xp->dbl(i) - xm->dbl(i));
				(*retObj)[0] += val;
			}
			
		}
		else{
			// indefinite integral
			retObj = new DvObject(0.0, Dims(), seqSize());

			val = this->record(0);
			val *= (xp->dbl(0) - xm->dbl(0));
			(*retObj)[0] = val;
			for(size_t i=1; i<length; i++){
				val = this->record(i);
				val *= (xp->dbl(i) - xm->dbl(i));
			
				(*retObj)[i] = (*retObj)[i-1];
				(*retObj)[i] += val;
			}
		}
	}
	else if(gapMethod == DV_LINEAR){
		// linear
	
		if(integralMethod == DV_DEF_INTEGRAL){
	
			// definite integral
			retObj = new DvObject(0.0, Dims(), 1);
			
			val = this->record(0);
			val *= (xp->dbl(0) - xm->dbl(0));
			(*retObj)[0] = val;
			double xi = x->dbl(0);
			double xpi = xp->dbl(0);

			for(size_t i=1; i<length; i++){
				// does gap filling linearly
				double xpi1 = xpi;
				xpi = xp->dbl(i);
				double xmi = xm->dbl(i);
				double h = xpi - xmi;
				double xgap = (xmi + xpi1) * 0.5;
				double hgap = xmi - xpi1; // zero if no gap
				double xi1 = xi;
				xi = x->dbl(i);
				hgap /= (xi - xi1);
				double scale = (xgap - xi1) * hgap; // zero if no gap as hgap zero
				double scale1 = (xi - xgap) * hgap; // zero if no gap as hgap zero
				
				val = this->record(i) * scale;
				val += this->record(i-1) * scale1;
				val += this->record(i)*h;

				(*retObj)[0] += val;
			}
		}
		else{
			// indefinite integral
			retObj = new DvObject(0.0, Dims(), length);

			val = this->record(0);
			val *= (xp->dbl(0) - xm->dbl(0));
			(*retObj)[0] = val;
			double xi = x->dbl(0);
			double xpi = xp->dbl(0);

			for(size_t i=1; i<length; i++){
				// does gap filling linearly
				double xpi1 = xpi;
				xpi = xp->dbl(i);
				double xmi = xm->dbl(i);
				double h = xpi - xmi;
				double xgap = (xmi + xpi1) * 0.5;
				double hgap = xmi - xpi1; // zero if no gap
				double xi1 = xi;
				xi = x->dbl(i);
				hgap /= (xi - xi1);
				double scale = (xgap - xi1) * hgap; // zero if no gap as hgap zero
				double scale1 = (xi - xgap) * hgap; // zero if no gap as hgap zero
				
				val = this->record(i) * scale;
				val += this->record(i-1) * scale1;
				val += this->record(i)*h;

				(*retObj)[i] = (*retObj)[i-1];
				(*retObj)[i] += val;
			}			
			
		}

	}

	if( retObj.is_nil() ) {
		retObj->error("[Integral] Integral or Gap method not understood");
		return retObj;
	}


  	// set xrefs
	
	retObj->copy_xrefs_from(*this);
	
	if(integralMethod == DV_DEF_INTEGRAL) retObj->delete_xref(DEPEND_0);
	
	// Units
	
	DvString attr_ptr = d0->getXrefText(UNITS);
	DvString attr_ptr_new = " ";
	if( !attr_ptr.empty() ){
     	if( !d0->is_time() ) attr_ptr_new += attr_ptr;
		else attr_ptr_new += "s"; // DvTime does time arithmetic in seconds whatever the metadata says
	
		attr_ptr = this->getXrefText(UNITS);
		if( !attr_ptr.empty() ){
			attr_ptr += attr_ptr_new;
			retObj->change_xref(UNITS, attr_ptr);
		}
	} 
	
	// SI_conversion
    attr_ptr = d0->getXrefText(SI_CONVERSION);
   
    if( !attr_ptr.empty() ){
	    DvObject_var dummy_ptr = new DvObject(1.0); // this is irrelevant as we only need metadata
		if( !d0->is_time() ) dummy_ptr->change_xref(SI_CONVERSION, attr_ptr);
		else dummy_ptr->change_xref(SI_CONVERSION, "1.0 > s");// DvTime does time arithmetic in seconds
		
		attr_ptr = this->getSICProduct(*dummy_ptr);
		if ( !attr_ptr.empty() ) {
			retObj->change_xref(SI_CONVERSION, attr_ptr);
 		}
	}
	
	if(this->xref_exists(FIELDNAM) ){
		DvObject_var xref = this->get_xref(FIELDNAM);
		DvString fName = "Integral(";
		fName +=  xref->asStr();
		if( d0->is_time() ) fName += ") dt";
		else fName += ") ds";
		retObj->change_xref(FIELDNAM, fName);
	}
	DvString method = "Histogram sum (gaps excluded)";
	if(gapMethod == DV_NN) method = "Histogram sum (gaps filled)";
	else if(gapMethod == DV_LINEAR) method = "Linear fit (gaps filled)";
	retObj->change_xref("Method", method);


	return retObj;
}

// --------------------- Differentiation ------------------------------

DvObject_var DvObject::differentiate(DvString derivMethod, double gapSize){
	if( !is_dbl() ){
		DvObject_var ret = new DvObject();
		ret->error("[Derivative] Not available for this data type");
		return ret;
	}
	
	size_t length = this->seqSize();
	if(length < 2){
		DvObject_var ret = new DvObject();
		ret->error("[Derivative] Input sequence too short");
		return ret;
	}
	
	// find gap size if not provided (1.5 times spacing)
	if(gapSize < 0) isRegular(gapSize);
	gapSize *= 1.5;
	
	DvObject_var d0 = getDep0();
	if(d0.is_nil() || d0->seqSize() != length) {
		DvObject *ret = new DvObject();
		ret->error("[Derivative] No useable DEPEND_0");
		return ret;
	}

	valarray <double> val(arraySize());
	
	DvObject_var h = new DvObject(0.0, length);
	h->dbl(0) = 0.0;
	for(size_t i=1; i<length; i++){
		h->dbl(i) = d0->asDouble(i) - d0->asDouble(i-1);
	}
	
	
	DvObject_var retObj;
	
	if(derivMethod == DV_DERIV_3PT){
	
		// 3 point estimate	
		retObj=new DvObject(0.0, Dims(), length);
		// do first entry as 2-point estimate
		val = (record(1) - record(0) ) / h->dbl(1);
		(*retObj)[0] = val;
		
		for(size_t i=1; i<length-1; i++){
			// use 2-point measure if gap
			if(h->dbl(i+1) > gapSize) {
				val = (record(i) - record(i-1) ) / h->dbl(i);
				(*retObj)[i] = val;
				i++;
				// and set start of next block with 2-point estimate
				if(i<length-1){
					 val = (record(i+1) - record(i) ) / h->dbl(i+1);
					(*retObj)[i] = val;
				}
			}
			else {
				val = (record(i+1) - record(i-1) ) / (h->dbl(i) + h->dbl(i+1)) ;
				(*retObj)[i] = val;
			}
		}
		// do last entry as 2-point estimate
		val = (record(length-1) - record(length-2) ) / h->dbl(length-1) ;
		(*retObj)[length-1] = val;

	}
	else {
		// 5 point estimate

		DvObject_var jobj;
		DvObject_var tags;
		
		// First test for regular data
		double hReg;
		bool isReg;
		isReg = isRegular(hReg);
		
		double h12 = 12.* hReg;
		double h2 = 2.* hReg;
		
		if(!isReg){

			// make D0 regular and multijoin
			
			double range = d0->asDouble(length-1) - d0->asDouble(0);
			length = (size_t) (range/hReg + 1.0);

			tags = new DvObject(d0.ptr(), length);
			tags->assign_regular(hReg);
			tags->copy_xrefs_from(d0);
			
			jobj = this->Join(tags, false);
			d0 = tags;
			
		}
		else {
			jobj = this;
		}
		
		retObj=new DvObject(0.0, Dims(), length);
		// Do derivative
					
		// do first entry as 2-point estimate
		val = (jobj->record(1) - jobj->record(0)) / hReg;
		(*retObj)[0] = val;
		
		// do second entry as 3-point estimate
		val = (jobj->record(2) - jobj->record(0)) / h2;
		(*retObj)[1] = val;
			
		for(size_t i=2; i<length-2; i++){
			// 5-point estimate
			val = (-jobj->record(i+2) + 8.*jobj->record(i+1) - 8.*jobj->record(i-1) + jobj->record(i-2) ) / h12; 
			(*retObj)[i] = val;
		}
		
		// do last entries as 2-point and 3-point estimate
		val = (jobj->record(length-1) - jobj->record(length-3) ) / h2;
		(*retObj)[length-2] = val;
		val = (jobj->record(length-1) - jobj->record(length-2) ) / hReg;
		(*retObj)[length-1] = val;
		
	}

  	// set xrefs
	
	retObj->copy_xrefs_from(*this);
	retObj->change_xref(DEPEND_0, d0);

	// Units
	
	DvString attr_ptr = d0->getXrefText(UNITS);
	DvString string_buf="";
	if( !attr_ptr.empty() ){
     	string_buf = " ";
     	if( !d0->is_time() ) string_buf += attr_ptr;
		else string_buf += "s"; // DVOS does time arithmetic in seconds whatever the metadata says
     	string_buf += "^-1";
	
		attr_ptr = getXrefText(UNITS);
		if( !attr_ptr.empty() ){
	
			attr_ptr += string_buf;
			retObj->change_xref(UNITS, attr_ptr);
		}
	} 
	
	// SI_conversion
    attr_ptr = d0->getSICInverse();
    if( !attr_ptr.empty() ){
	
	    DvObject_var inv_ptr = new DvObject(1.0); // this is irrelevant as we only need metadata
		
		if( !d0->is_time()) inv_ptr->change_xref(SI_CONVERSION, attr_ptr);
		else inv_ptr->change_xref(SI_CONVERSION, "1.0 > s^-1");
		
		attr_ptr = getSICRatio(*inv_ptr);
		if ( !attr_ptr.empty() ) {
			retObj->change_xref(SI_CONVERSION, attr_ptr);
 		}
	}
	
	if(xref_exists("FIELDNAM") ){
		DvString xref = getXrefText("FIELDNAM");
		if(!xref.empty()){
			DvString fName = "d ";
			fName +=  xref;
			if(d0->is_time()) fName += "/ dt";
			else fName += "/ ds";
			retObj->change_xref("FIELDNAM", fName);
		}
	}

	retObj->change_xref("Method", derivMethod);

	return retObj;

}

DvObject_var DvObject::regression(){
	
	if(this->is_nil()) return new DvObject();

	DvObject_var dep = this->getDep0();
	if(dep.is_nil()) return new DvObject();

	return this->regression(dep);
}

DvObject_var DvObject::regression(DvObject_var &dep){
	
	if(dep->is_nil() || dep->arraySize() > 1 || (dep->not_dbl() && dep->not_time() )  ) {
		DvObject_var ret = new DvObject();
		ret->error("[regression] X data type not valid");
		return ret;
	}
	
	if(this->is_nil() || this->arraySize() > 1 || this->not_dbl()  ) {
		DvObject_var ret = new DvObject();
		ret->error("[regression] Y data type not valid");
		return ret;
	}
	
	DvString dep_units_str = dep->getXrefText(UNITS);	
	DvString dep_SIc_str = dep->getXrefText(SI_CONVERSION);	
	
	DvObject_var xobj = dep.ptr();
	
	if(dep->is_time()){
		xobj = new DvObject(0.0, dep->seqSize());

		for(size_t i=0; i<dep->seqSize(); i++){
			xobj->dbl(i) = dep->time(i) - dep->time(0);
		}
		dep_units_str = DvString("s");
		dep_SIc_str = DvString("1>s");
		
		xobj->change_xref(SI_CONVERSION, dep_SIc_str);
	}
	
	
	if(this->seqSize() != xobj->seqSize()) {
		DvObject_var ret = new DvObject();
		ret->error("[regression] X and Y data different lengths");
		return ret;
	}
	
	int n = this->seqSize();
	double sigxy = 0;
	double sigxx = 0;
	double sigx = 0;
	double sigy = 0;
	double sigyy = 0;
	double x, y;
	
	double min = xobj->dbl(0);
	double max = xobj->dbl(xobj->seqSize()-1);
	
	for(int i=0; i<n; i++){
		x = xobj->dbl(i);
		if(min > x) min = x;
		if(max < x) max = x;
		
		y = this->dbl(i);
		sigxy += y * x;
		sigxx += x * x;
		sigx  += x;
		sigy  += y;
		sigyy += y * y;
	}
	double m = (n*sigxy - sigx * sigy) /(n*sigxx - sigx * sigx);
	double b = (sigy - m * sigx) / n;
	double r = (n*sigxy - sigx * sigy) / std::sqrt( (n*sigxx - sigx*sigx) * (n*sigyy - sigy * sigy) );
	
	DvObject_var ret = new DvObject(0.0, 2);
	ret->dbl(0) = m*min + b;
	ret->dbl(1) = m*max + b;

	DvObject_var dep0;
	if(dep->is_time()){	
		dep0 = new DvObject(dep->time(0)+min, 2); 
		dep0->time(1) = dep->time(0)+max;
	}
	else{	
		dep0 = new DvObject(min, 2); 
		dep0->dbl(1) = max; 
	}
	if(!dep_units_str.empty()) dep0->change_xref(UNITS, dep_units_str);
	if(!dep_SIc_str.empty()) dep0->change_xref(SI_CONVERSION, dep_SIc_str);	
	
	DvObject_var mObj = new DvObject(m);
	DvObject_var bObj = new DvObject(b);
	DvObject_var rObj = new DvObject(r);

	if(this->xref_exists(UNITS)) {
		DvObject_var units_xref = this->get_xref(UNITS);
		ret->change_xref(UNITS, units_xref);
		bObj->change_xref(UNITS, units_xref);
		DvString units_str = units_xref->asStr(0);
		
		if(!dep_units_str.empty() && !units_str.empty() ) {
			DvString slope_units = units_str;
			slope_units += "/";
			slope_units += dep_units_str;
			mObj->change_xref(UNITS, slope_units);
		}
			
	}
	if(this->xref_exists(SI_CONVERSION)) {
		DvObject_var SIc_xref = this->get_xref(SI_CONVERSION);
		ret->change_xref(SI_CONVERSION, SIc_xref);
		bObj->change_xref(SI_CONVERSION, SIc_xref);
		
		DvString slope_SIc = this->getSICRatio(*xobj);
		mObj->change_xref(SI_CONVERSION, slope_SIc);
	}
		
	ret->change_xref(DEPEND_0, dep0);
	ret->change_xref("Slope", mObj);
	ret->change_xref("Intercept", bObj);
	ret->change_xref("CoeffCorrel", rObj);

	if(this->xref_exists("FIELDNAM") ){
		DvString xref = this->getXrefText(FIELDNAM);
		xref += " linear fit";
		ret->change_xref("FIELDNAM", xref);
	}
	 
	return ret;

}

DvObject_var DvObject::mergeWith(DvObject_var &dobj){

	// if either object is empty, just return the other
	if(this->is_nil())  return new DvObject(dobj.ptr());
	if(dobj.is_nil())  return new DvObject(*this);
	
	// reject object with no DEPEND_0 as this is merge on ordered dependency
	DvObject_var d0_1 = this->getDep0();
	if(d0_1.is_nil()) 
	{
		DvObject_var retObj = new DvObject(dobj);
		retObj->error("[merge] this object has no Depend_0");
		return retObj;
	}
	DvObject_var d0_2 = dobj->getDep0();
	if(d0_2.is_nil())  
	{
		DvObject_var retObj = new DvObject(this);
		retObj->error("[merge] target object has no Depend_0");
		return retObj;
	}
	
	int nD1=0, nD2=0;
	nD1 = d0_1->seqSize();
	nD2 = d0_2->seqSize();
		
	if( !this->sameUnits(*dobj) )
	{
		DvObject_var retObj = new DvObject(this);
		retObj->error("[merge] units different");
		return retObj;
	}
	
	if( !this->sameFrame(*dobj))
	{
		DvObject_var retObj = new DvObject(this);
		retObj->error("[merge] frames different");
		return retObj;
	}
	
	if( this->arraySize() != dobj->arraySize())
	{
		DvObject_var retObj = new DvObject(this);
		retObj->error("[merge] dimensions different");
		return retObj;
	}

	if( (this->is_dbl()   && dobj->not_dbl())  || 
		(this->is_int()   && dobj->not_int())  || 
		(this->is_str()   && dobj->not_str())  || 
		(this->is_time()  && dobj->not_time()) || 
		(this->is_event() && dobj->not_event()) )
	{
		DvObject_var retObj = new DvObject(this);
		retObj->error("[merge] data types different");
		return retObj;
	}
	

	// create mask for final object with input source for each record
	// 0= empty record due to earlier duplicates
	// 1= dobj1
	// 2= dobj2
	// 3= duplicate, take either but increment both on read
	vector <int> pairMask;
	pairMask.assign(nD1 + nD2, 0); // initialise to zero in case of dropped duplicates
	
	int n1=0;
	int n2=0;
	int nMask=0;

	while(n1<nD1 && n2<nD2)
	{
		if(d0_1->asDouble(n1) == d0_2->asDouble(n2)) 
		{
			pairMask[nMask] = 3;
			n1++;
			n2++;
			nMask++;
		}
		else if(d0_1->asDouble(n1) < d0_2->asDouble(n2))
		{
			pairMask[nMask] = 1;
			n1++;
			nMask++;
		}
		else
		{
			pairMask[nMask] = 2;
			n2++;
			nMask++;
		}
	}

	// finish off dobj2 if dobj1 ran out first
	while(n2<nD2)
	{
		pairMask[nMask] = 2;
		n2++;
		nMask++;
	}
	// finish off dobj1 if dobj2 ran out first
	while(n1<nD1)
	{
		pairMask[nMask] = 1;
		n1++;
		nMask++;
	}
	
	// we now have a mask with the order of records to use, even if one object is nested inside the other.
	// Do merge for object and t dependent xrefs 
	
	
	return this->mergeByMask(nMask, pairMask, dobj);
}

// -------------------------------------------------------------------

DvObject_var DvObject::mergeByMask(size_t nMask, vector <int> &pairMask, DvObject_var &dobj)
{
	if(this->is_nil()) return new DvObject(dobj.ptr());
	if(dobj.is_nil()) return new DvObject(this);
	
	int n1=0;
	int n2=0;
	
	// create new object from inputs according to mask
	
	DvObject_var retObj = new DvObject(this, nMask); // constructs object with nMask records like first of this
	
	// assemble record by record			
	for(size_t i=0; i<nMask; i++)
	{
		if(pairMask[i] == 1)
		{
			(*retObj)[i] = (*this)[n1];
			n1++;
		}
		else if(pairMask[i] == 2) 
		{
			(*retObj)[i] = dobj[n2];
			n2++;
		}
		else 
		{
			(*retObj)[i] = (*this)[n1];
			n1++;
			n2++;
		}
	}

	retObj->copy_xrefs_from(*this);
	
	// modify t-dep xrefs (inc time tags)
	DvNode *Xnode = this->first_xref();
	while(Xnode)
	{
		if( Xnode->obj()->seqSize() > 1 &&  Xnode->obj()->seqSize() == this->seqSize() ){
			DvObject_var xref_obj2 = dobj->get_xref(Xnode->name());
			if(xref_obj2.is_ok()){
				DvObject_var merged_xref_obj = Xnode->obj()->mergeByMask(nMask, pairMask, xref_obj2);
			
				retObj->change_xref(Xnode->name(), merged_xref_obj);
			}
		}
		Xnode = Xnode->next;
	}
	
	return retObj;
}

// -------------------------------------------------------------------

void DvObject::stripDuplicates(double tolerance)
{
	// removes duplicate DEPEND_0 values (time tags or scalar) (within t1 - t2 < tolerance)

	
	DvMask msk(seqSize()); // initialises to true

	DvObject_var d0 = getDep0();
	if(d0.is_nil())
	{
		error("stripDuplicates  No DEPEND_0 found");
		return;
	}
	
	if(d0->seqSize() != seqSize())
	{
		error("stripDuplicates  Depend_0 and data sequence lengths differ");
		return;
	}
	
	// remove duplicates
		
	for(size_t i=1; i<d0->seqSize(); i++)
	{
		if(fabs(d0->asDouble(i) - d0->asDouble(i-1)) <= tolerance )
		{
			msk[i] = false;
		}
	}

	this->apply_mask(msk); // does object and xrefs
	

	return;
}

// -------------------------------------------------------------------

DvObject_var DvObject::makeMonotonic()
{
	// shuffles elements so tt are monotonic
	if( seqSize() <= 1) return new DvObject(*this);
	
	// create sequence holding current ordering
	valarray <size_t> order((size_t)0, seqSize());
	
	DvObject_var d0 = getDep0();
	if(d0.is_nil())
	{
		DvObject_var retObj = new DvObject(*this);
		retObj->error("makeMonotonic:  No DEPEND_0 found");
		return retObj;
	}
	
	if(d0->seqSize() != seqSize())
	{
		DvObject_var retObj = new DvObject(*this);
		retObj->error("makeMonotonic:  Depend_0 and data sequence lengths differ");
		return retObj;
	}
	
	if(d0.is_ok()){	
		// sort order, algorithm is optimised for nearly monotonic
		order[0] = 0;
		for(int i=1; i<seqSize(); i++)
		{
			// compare next d0 with latest time found
			double diff = d0->asDouble(i) - d0->asDouble(order[i-1]);
			if(diff >=  0) order[i] = i;
			else
			{
				order[i] = order[i-1];
				order[i-1] = i;
				// place in correct place
				for(int j=i-1; j>0; j--)
				{
					if(d0->asDouble(order[j]) - d0->asDouble(order[j-1]) < 0.0) 
					{	
						int hold = order[j];
						order[j] = order[j-1];
						order[j-1] = hold;
					}
					else break;
				}
			}
		}
	}

		
	// the vector order now holds the indices of the elements in time order
	// re-order the data object and all its sequence ordered metadata

	DvObject_var retObj = this->applyOrder(order);

	
	retObj->copy_xrefs_from(*this);
	
	retObj->change_xref("Re-ordered", "Made Monotonic in QSAS");


	DvNode *Xnode = this->first_xref();
	while(Xnode)
	{
		if( Xnode->obj()->seqSize() > 1 &&  Xnode->obj()->seqSize() == this->seqSize() ){
			DvObject_var newXref = Xnode->obj()->applyOrder(order);
			
			retObj->change_xref(Xnode->name(), newXref);
		}
		Xnode = Xnode->next;
	}
	
	return retObj;
}

DvObject_var DvObject::applyOrder(valarray <size_t> &order)
{
	
	// create new object from inputs according to mask
	DvObject_var newObj = new DvObject(*this);

	for(size_t i=0; i<seqSize(); i++)
	{
	  if(i != order[i]) (*newObj)[i] = (*this)[order[i]]; // uses DvRecord operator
	}		

	return newObj;

	
}

