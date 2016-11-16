//
//  DvComponent.cc
//
//  Dvos -  Data variable object classes
//
//  Vector components 
//
//  Tony Allen
//  Nov 2012
//

#include <iostream>

#include "DvObj.h"
#include <valarray>
#include <vector>
#include <cmath>




vector <DvString> DvObject::available(DvString constraint){
	
	vector <DvString> ret;

	if(constraint == DV_SCALAR || constraint == DV_NONE){
		if(isVectorXYZ()) {
			ret.push_back(DvString(DV_RMAG));
			ret.push_back(DvString(DV_X));
			ret.push_back(DvString(DV_Y));
			ret.push_back(DvString(DV_Z));
			ret.push_back(DvString(DV_LAT_DEG));
			ret.push_back(DvString(DV_THETA_DEG));
			ret.push_back(DvString(DV_PHI_DEG));
			ret.push_back(DvString(DV_LAT_RAD));
			ret.push_back(DvString(DV_THETA_RAD));
			ret.push_back(DvString(DV_PHI_RAD));
			if(getFrame() == "gse") ret.push_back(DvString(DV_LOCAL_T));
		}
		else{
			if(arraySize() > 1) {
				if(ret.size() > 0) ret.push_back(DvString(DV_SEPARATOR));
				ret.push_back(DvString(DV_ARRAY_REDUCE));
			}
		}
	}
	if(constraint == DV_TIME || constraint == DV_NONE){
		if(ret.size() > 0) ret.push_back(DvString(DV_SEPARATOR));

	}
	
	if(seqSize()>1){
		if(ret.size() > 0) ret.push_back(DvString(DV_SEPARATOR)); 
		ret.push_back(DvString(DV_SUBSEQUENCE)); 
	}
	
	if(ret.size() > 0) ret.push_back(DvString(DV_SEPARATOR)); 
	ret.push_back(DvString(DV_METADATA)); // always possible

	return ret;
}

// ------------------------

DvObject_var DvObject::getObjComponent(DvString comp){
	
	if( comp.empty())  return new DvObject();
  
	DvObject_var dobj;

	DvString compStr = comp.betweenIncl("[", "]");

	if(compStr.empty()) 
		dobj = this;
	else if (compStr.contains("[=>") )
	{
		// put first in case metadata name mimics other components
		compStr = compStr.between("[=>", "]");
		compStr = compStr.trimmed();
		dobj = this->get_xref(compStr);       
	}
	else if (compStr.contains(DV_RMAG) ) 
		dobj = this->getVectorComponent(DV_RMAG );
	else if (compStr.contains(DV_X) )
		dobj = this->getVectorComponent(DV_X);
	else if (compStr.contains(DV_Y) )
		dobj = this->getVectorComponent(DV_Y );
	else if (compStr.contains(DV_Z) )
		dobj = this->getVectorComponent(DV_Z );
	else if (compStr.contains(DV_LAT_DEG) )
		dobj = this->getVectorComponent(DV_LAT_DEG );
	else if (compStr.contains(DV_LAT_RAD) )
		dobj = this->getVectorComponent(DV_LAT_RAD );
	else if (compStr.contains(DV_THETA_DEG) )
		dobj = this->getVectorComponent(DV_THETA_DEG );
	else if (compStr.contains(DV_THETA_RAD) )
		dobj = this->getVectorComponent(DV_THETA_RAD );
	else if (compStr.contains(DV_PHI_DEG) )
		dobj = this->getVectorComponent(DV_PHI_DEG );
	else if (compStr.contains(DV_PHI_RAD) )
		dobj = this->getVectorComponent(DV_PHI_RAD);
	else if (compStr.contains(DV_LOCAL_T) )
		dobj = this->getVectorComponent(DV_LOCAL_T);
	else if (compStr.contains(DV_T_START) )
		dobj = this->getTimeComponent(DV_T_START );
	else if (compStr.contains(DV_T_CENTRE) )
		dobj = this->getTimeComponent(DV_T_CENTRE );
	else if (compStr.contains(DV_T_END) )
		dobj = this->getTimeComponent(DV_T_END);
	else if (compStr.contains("[") )
	{
		dobj = this->getArrayComponent(compStr);
	}
	else dobj = this;
	
	compStr = comp.between(" <", ">"); // without brackets this time

  	if(!compStr.empty()) {
		if( compStr.contains("Ave") ) dobj = dobj->getTimeAve(compStr.after("Ave"));
		else dobj = dobj->getTimeSubset(compStr);
	}
	
	return dobj;
}


DvObject_var DvObject::getVecFromRTP(){

	// Reduce RTP to Cartesian (QSAS storage standard)
	// theta in range [0,pi), phi in range [0,2pi)
	
	double toRad = std::atan(1.0)/45.; // default is deg2rad
	
	if(  !is_dbl() || arraySize() != 3 ) {
		DvObject_var ret = new DvObject();
		ret->error("[Vector conversion] Not valid for this data type");
		return ret;
	}
	
	DvObject_var D1 = get_xref(DEPEND_1);
	if( D1.is_ok() && D1->is_str() ){
		for(size_t i=0; i< D1->arraySize(); i++) if( D1->str(i).contains("rad") ) toRad = 1.0;
	}
	
	DvObject_var ret = new DvObject(this); // duplicate

	// convert elements
	for(size_t n=0; n<seqSize(); n++){
	    double xy = this->dbl(3*n) * std::sin(toRad * this->dbl(3*n+1)); // R sin(Th)
		ret->dbl(3*n) = xy * std::cos(toRad * dbl(3*n+2)); // xy cos(Phi)
		ret->dbl(3*n+1) = xy * std::sin(toRad * dbl(3*n+2)); // xy sin(Phi)
		ret->dbl(3*n+2) = dbl(3*n) * std::cos(toRad * dbl(3*n+1)); // R cos(Th)
	}
	
	// fix rtp to xyz
	DvNode *Xnode = ret->first_xref();
	while(Xnode){
		if(Xnode->obj().is_ok() && Xnode->obj()->is_str() ){
			for(size_t j=0; j<Xnode->obj()->totalSize(); j++) Xnode->obj()->str(j).replace("RTP", "xyz");
		}
	}
	
	return ret;
}

DvObject_var DvObject::getVecFromRLP(){

	// Reduce RLP to Cartesian (QSAS storage standard)
	// latitude in range [-pi/2,pi/2), phi in range [0,2pi)

	double toRad = std::atan(1.0)/45.; // default is deg2rad
	
	if(  !is_dbl() || arraySize() != 3 ) {
		DvObject_var ret = new DvObject();
		ret->error("[Vector conversion] Not valid for this data type");
		return ret;
	}
	
	DvObject_var D1 = get_xref(DEPEND_1);
	if( D1.is_ok() && D1->is_str() ){
		for(size_t i=0; i< D1->arraySize(); i++) if( D1->str(i).contains("rad") ) toRad = 1.0;
	}
	
	DvObject_var ret = new DvObject(this); // duplicate

	// convert elements
	for(size_t n=0; n<seqSize(); n++){
	    double xy = dbl(3*n) * std::cos(toRad * dbl(3*n+1)); // R cos(Lat)
		ret->dbl(3*n) = xy * std::cos(toRad * dbl(3*n+2)); // xy cos(Phi)
		ret->dbl(3*n+1) = xy * std::sin(toRad * dbl(3*n+2)); // xy sin(Phi)
		ret->dbl(3*n+2) = dbl(3*n) * std::sin(toRad * dbl(3*n+1)); // R sin(Lat)
	}
	
	// fix rlp to xyz
	DvNode *Xnode = ret->first_xref();
	while(Xnode){
		if(Xnode->obj().is_ok() && Xnode->obj()->is_str() ){
			for(size_t j=0; j<Xnode->obj()->totalSize(); j++) 
			{
				Xnode->obj()->str(j) = Xnode->obj()->str(j).replace("rlp", "xyz");
			}
		}
		Xnode = Xnode->next;
	}

	return ret;
}


DvObject_var DvObject::getVectorComponent(DvString comp){

	// Three vectors are stored inside QSAS in cartesians
	double rad2deg = 45./std::atan(1.0);
	double rad2hr = 3./std::atan(1.0);
	
	
	if( !is_dbl() ) {
		DvObject_var ret = new DvObject();
		ret->error("[Vector comp] Not valid for this data type");
		return ret;
	}
	
	size_t len = seqSize();
	valarray <double> val(0.0, len);
	valarray <double> x(Ddata[slice(0, len, 3)]);
	valarray <double> y(Ddata[slice(1, len, 3)]);
	valarray <double> z(Ddata[slice(2, len, 3)]);
	
	if(comp == DV_RMAG){
		x *= x;
		y *= y;
		z *= z;
		z += x;
		z += y;
		val = std::sqrt(z);
	}
	else if(comp == DV_LAT_RAD){
		// Obtain latitude(l) { [pi/2,-pi/2) }, from Cartesian coordinates
		// latitude = pi/2 - theta, where theta is angular displacement in radians from 
		// the positive z-axis
	
		x *= x;
		y *= y;
		y = std::sqrt(y+x);
		val = std::atan2(z, y);
	}
	else if(comp == DV_X){
		// X component
	
		val = x;
	}
	else if(comp == DV_Y){
		// X component
	
		val = y;
	}
	else if(comp == DV_Z){
		// X component
	
		val = z;
	}
	else if(comp == DV_LAT_DEG){
		// Obtain latitude(l) { [pi/2,-pi/2) }, from Cartesian coordinates
		// latitude = pi/2 - theta, where theta is angular displacement in deg from 
		// the positive z-axis
	
		x *= x;
		y *= y;
		y = std::sqrt(y+x);
		val = std::atan2(z, y);
		val *= rad2deg;
	}
	else if(comp == DV_THETA_RAD){
		// Obtain theta(t) from Cartesian
		// Theta is angular displacement in radians from the positive z-axis [0,pi]
		x *= x;
		y *= y;
		y = std::sqrt(y+x);
		val = std::atan2(y, z);
	}
	else if(comp == DV_THETA_DEG){
		// Obtain theta(t) from Cartesian
		// Theta is angular displacement in deg from the positive z-axis [0,pi]
		x *= x;
		y *= y;
		y = std::sqrt(y+x);
		val = std::atan2(y, z);
		val *= rad2deg;
	}
	else if(comp == DV_PHI_RAD){
		// Obtain phi(p) from Cartesian coordinates
		// phi_out is angular displacement in radians from the positive x-axis [0,2pi)
   				
		val = atan_yx(y, x);
	}
	else if(comp == DV_PHI_DEG){
		// Obtain phi(p) from Cartesian coordinates
		// phi_out is angular displacement in deg from the positive x-axis [0,2pi)
   				
		val = atan_yx(y, x);
		val *= rad2deg;
	}
	else if(comp == DV_LOCAL_T){
		// Obtain local time (alias for from phi(p)) from Cartesian coordinates
		// Noon is at y=0, x > 0
		// Midnight is at y=0, x<0
		// dusk flank has at y>0, dawn has y<0
   				
		val = std::atan2(y, x);
		val *= rad2hr;
		val += 12.0;

	}
	
	DvObject_var ret = new DvObject(val);
	ret->copy_xrefs_from(*this);

	// Modify xrefs
	
	DvObject_var label_axis = get_xref(LABEL_1);
	if (label_axis.is_nil()) label_axis = get_xref(LABL_PTR_1);
	if (label_axis.is_nil()) label_axis = get_xref(LABLAXIS);
	if (label_axis.is_nil()) label_axis = get_xref(FIELDNAM);
	if (label_axis.is_nil()) label_axis = get_xref("ObjName");

	DvString label = label_axis->getCommonText();

	ret->delete_xref(LABEL_1);
	ret->delete_xref(LABL_PTR_1);
	
	DvString fieldnam = getXrefText(FIELDNAM);
		 
	if(label.empty()) label = fieldnam;
	fieldnam += " (";
	fieldnam += comp;
	fieldnam += ")";
	
	 
	if(!fieldnam.empty()) ret->change_xref(FIELDNAM, fieldnam);
	 
	if(label.empty()){
		// last gasp, try ObjectName created by qie
		label = getXrefText(OBJECT_NAME);
	}

  	 
	DvString frame("scalar>na");
	ret->change_xref(FRAME,frame);
	 
	DvString si_conversion;
	DvString units;
	
	if(comp==DV_LAT_DEG || comp==DV_THETA_DEG || comp==DV_PHI_DEG)
	{
		si_conversion="1>degree";
		units="deg";

		label += "_";
		label += comp;
		label = label.replace(" (deg)", "");
	}
	else if(comp==DV_LAT_RAD || comp==DV_THETA_RAD || comp==DV_PHI_RAD)
	{
		si_conversion="1>radian";
		units="rad";

		label += "_";
		label += comp;
		label = label.replace(" (rad)", "");
	}
	else if(comp==DV_RMAG)
	{
		si_conversion = getXrefText(SI_CONVERSION);
		units = getXrefText(UNITS);
		label = DvString("|") + label;
		label += "|";

	}
	else if(comp==DV_LOCAL_T)
	{
		si_conversion="360.0>s";
		units="hr";

		label="Local Time";
	}
	else
	{
		si_conversion = getXrefText(SI_CONVERSION);
		units = getXrefText(UNITS);
		label += "_";
		label += comp;
	}
     
	if(si_conversion.empty()) si_conversion = "1>(unknown)";
	ret->change_xref(SI_CONVERSION, si_conversion);
		
	if( !units.empty() ) ret->change_xref(UNITS, units);

	if(!label.empty()) ret->change_xref(LABLAXIS, label);

	if( !isVectorXYZ() ) {
		ret->error("[Vector comp] Object not a vector");
		return ret;
	}

	return ret;
	
}

DvObject_var DvObject::getTimeComponent(DvString comp){
	
	// no component processing
	if(comp.empty() && this->is_time() )  return this;

	// Look for single time requests
	if (comp.contains("[=>") )
	{
		comp = comp.between("[=>", "]");
		comp = comp.trimmed();

		DvObject_var dobj = this->get_xref(comp); 
		return dobj;   
	}
	if(comp.contains(DV_T_START)) {
		DvEvent event = this->getTimeInterval(comp); // this will do subsequence if needed
		DvTime t(event.start_abs()); 
		return new DvObject(t);
	}
	if(comp.contains(DV_T_CENTRE)){
		DvEvent event = this->getTimeInterval(comp); // this will do subsequence if needed
		DvTime t(event.centre()); 
		return new DvObject(t);
	}
	if(comp.contains(DV_T_END)){
		DvEvent event = this->getTimeInterval(comp); // this will do subsequence if needed
		DvTime t(event.end_abs()); 
		return new DvObject(t);
	}

	// return time
	
	DvObject_var tt = getTimeTags(); // gets this if it is already time
	
	tt = tt->getTimeSubset(comp.between("<",">"));
	
	return tt;
	
}

DvObject_var DvObject::getIntervalComponent(DvString name){

	// initialise var pointer here (will be nil if resolve fails)
	if(name.empty())  return new DvObject();

	DvString ivlStr = name.between("<", ">"); 
	
	if( this->is_event() ) 
	{
		if(ivlStr.empty()) return this; 
		
		int nEvent = ivlStr.toInt();
		
		if(nEvent >= (int)this->seqSize()) nEvent = (int)this->seqSize() - 1;
		
		return  this->subSequence(nEvent);
		
	}
	
	DvEvent event = this->getTimeInterval(name.betweenIncl("<", ">"));
	return  new DvObject(event);
	
}


DvObject_var DvObject::getArrayComponent(DvString comp){

	// comp of the form, e.g.  " [ AtT(1,8) Sum(1,16) (1,31) ]"
	// must subtract 1 from everything to get (0, N-1) form used in c++
	// use var as *object gets replaced
	DvObject_var dobj = this;
	
	comp = comp.trimmed();
	
	DvString opt = comp.before("(");
	DvString next = comp.between ("(", ")");
	
	size_t d=0; // dimension (of new object) being processed
	// Note d is only incremented for 'elements' option as others remove dimension
	
	while( !next.empty()){
		size_t from;
		size_t to;
		if(next.contains(',') ){
			from = next.front(",").toInt() - 1;
			to= next.back(",").toInt() - 1;
		}
		else{
			from = next.toInt() - 1;
			to = from;
		}
		

		if(dobj->dims.size() <= d || to > dobj->dims[d] || from > dobj->dims[d] ){
			this->error("[array component] out of range");
			return this;
		}
		
		if(opt.contains("Sum") ){
			//Sum
			dobj = dobj->sumDim(d, from, to );
		}
		else if(opt.contains("Ave")){
			//Average
			dobj = dobj->averageDim(d, from, to);
		}
		else if(opt.contains("AtT")){
			// Alias to Time
			dobj = dobj->unwrapDim(d, from, to);
		}
		else if(from == to){
			// slice at
			dobj = dobj->sliceDim(d, from);
		}
		else {
			// elements in range

			dobj = dobj->subsetDim(d, from, to);
			d++; 
		}
		// get next token
		comp = comp.after(")"); // move to start of next token
		opt = comp.before("(");
		next = comp.between("(", ")");

	}
	
	return dobj;

}

DvObject_var DvObject::getTimeSubset(DvString comp){
	
	DvString from = comp.before(",");
	DvString to = comp.after(",");

	if(from.empty() ) from = comp;
	if(to.empty()) return this->subSequence(from.toInt()); // single record
	
	if(from.toInt() == 1 && to.toInt() == (int)this->seqSize())  return this;
	
	DvObject_var ret = this->subSequence(from.toInt(), to.toInt());

	return ret;

}

DvObject_var DvObject::getTimeAve(DvString comp){
	
	DvString from = comp.before(",");
	DvString to = comp.after(",");

	if(from.empty() ) from = comp;
	if(to.empty()) return this->subSequence(from.toInt()); // single record
	
	int n0 = from.toInt();
	int n1 = to.toInt();
	int len = n1 - n0;
	if(len < 0) len *= -1;
	len += 1;

	if(len == 0){
		DvObject_var ret = new DvObject();
		ret->error("Zero records in time average");
		return ret;
	}
	
	if(this->is_dbl()) {
		valarray <double> ave(0.0, this->arraySize());
		for(int n=n0; n<=n1; n++){
			for(size_t i=0; i<this->arraySize(); i++){
				ave[i] += this->dbl(n,i);
			}
		}

		ave /= (double)len;
		
		DvObject_var ret = new DvObject(ave, this->Dims(), 1);
		ret->copy_xrefs_from(*this);
	
		// Also do xrefs if sequence
		DvNode *xref = this->first_xref();
		while(xref)
		{
			if(xref->obj()->seqSize() == this->seqSize()) {
				DvObject_var xref_result = xref->obj()->getTimeAve(comp);
				ret->change_xref(xref->name(), xref_result);
			}
			xref = xref->next;
		} 

		return ret;
	}
	else if(this->is_int()) {
		valarray <int> ave(0, this->arraySize());
		for(int n=n0; n<=n1; n++){
			for(size_t i=0; i<this->arraySize(); i++){
				ave[i] += this->itg(n,i);
			}
		}

		ave /= len;
	
		DvObject_var ret = new DvObject(ave, this->Dims(), 1);
		ret->copy_xrefs_from(*this);
	
		// Also do xrefs if sequence
		DvNode *xref = this->first_xref();
		while(xref)
		{
			if(xref->obj()->seqSize() == this->seqSize()) {
				DvObject_var xref_result = xref->obj()->getTimeAve(comp);
				ret->change_xref(xref->name(), xref_result);
			}
			xref = xref->next;
		} 

		return ret;
	}
	else if(this->is_time()) {
		valarray <DvTime> ave(DvTime(), this->arraySize());
		for(int n=n0; n<=n1; n++){
			for(size_t i=0; i<this->arraySize(); i++){
				ave[i] += this->time(n,i);
			}
		}

		for(size_t i=0; i<this->arraySize(); i++){
			ave[i] /= (double) len;
		}
	
		DvObject_var ret = new DvObject(ave, this->Dims(), 1);
		ret->copy_xrefs_from(*this);
	
		// Also do xrefs if sequence
		DvNode *xref = this->first_xref();
		while(xref)
		{
			if(xref->obj()->seqSize() == this->seqSize()) {
				DvObject_var xref_result = xref->obj()->getTimeAve(comp);
				ret->change_xref(xref->name(), xref_result);
			}
			xref = xref->next;
		} 

		return ret;
	}
	else{
		DvObject_var ret = new DvObject();
		ret->error("Time average not supported for data type");
		return ret;
	}
}

DvObject_var DvObject::sum(){
	
	if(this->is_dbl()) {
		valarray <double> sumAll(0.0, this->arraySize());
		for(size_t n=0; n<seqSize(); n++){
			for(size_t i=0; i<this->arraySize(); i++){
				sumAll[i] += this->dbl(n,i);
			}
		}
		DvObject_var ret = new DvObject(sumAll, this->Dims());
		ret->copy_xrefs_from(*this);

		return ret;
	}
	else if(this->is_int()) {
		valarray <int> sumAll(0.0, this->arraySize());
		for(size_t n=0; n<seqSize(); n++){
			for(size_t i=0; i<this->arraySize(); i++){
				sumAll[i] += this->itg(n,i);
			}
		}
		DvObject_var ret = new DvObject(sumAll, this->Dims());
		ret->copy_xrefs_from(*this);

		return ret;
	}
	else if(this->is_time())  {
		// only any use if later normalised
		valarray <DvTime> sumAll(0.0, this->arraySize());
		for(size_t n=0; n<seqSize(); n++){
			for(size_t i=0; i<this->arraySize(); i++){
				sumAll[i] += this->time(n,i);
			}
		}
		DvObject_var ret = new DvObject(sumAll, this->Dims());
		ret->copy_xrefs_from(*this);

		return ret;
	}
	else{
		DvObject_var ret = new DvObject();
		ret->error("Sum not supported for data type");
		return ret;
	}
}


DvObject_var DvObject::sumRecEvery(size_t stride){

	DvObject_var ret;
	size_t newSeqLen = seqSize() / stride;
	size_t arrayLen = this->arraySize();
	
	if(this->is_dbl()) {
		ret = new DvObject(0.0, this->Dims(), newSeqLen);
		for(size_t n=0; n<newSeqLen; n++){
			for(size_t j=0; j<stride; j++){
				for(size_t i=0; i<arrayLen; i++){
					ret->dbl(n,i) += this->dbl(n*stride+j,i);
				}
			}
		}
	}
	else if(this->is_int()) {
		ret = new DvObject((int)0, this->Dims(), newSeqLen);
		for(size_t n=0; n<newSeqLen; n++){
			for(size_t j=0; j<stride; j++){
				for(size_t i=0; i<arrayLen; i++){
					ret->itg(n,i) += this->itg(n*stride+j,i);
				}
			}
		}
	}
	else if(this->is_time()){
		ret = new DvObject(DvTime(), this->Dims(), newSeqLen);
		for(size_t n=0; n<newSeqLen; n++){
			for(size_t j=0; j<stride; j++){
				for(size_t i=0; i<arrayLen; i++){
					ret->time(n,i) += this->time(n*stride+j,i);
				}
			}
		}
	}
	else{
		ret = new DvObject();
		ret->error("Sum records not supported for data type");
		return ret;
	}
	
	ret->copy_xrefs_from(*this);
	return ret;

}

DvObject_var DvObject::aveRecEvery(size_t stride){

	DvObject_var ret;
	size_t newSeqLen = seqSize() / stride;
	size_t arrayLen = this->arraySize();
	
	if(this->is_dbl()) {
		ret = new DvObject(0.0, this->Dims(), newSeqLen);
		for(size_t n=0; n<newSeqLen; n++){
			for(size_t j=0; j<stride; j++){
				for(size_t i=0; i<arrayLen; i++){
					ret->dbl(n,i) += this->dbl(n*stride+j,i);
				}
			}
		}
	}
	else if(this->is_int()) {
		ret = new DvObject((int)0, this->Dims(), newSeqLen);
		for(size_t n=0; n<newSeqLen; n++){
			for(size_t j=0; j<stride; j++){
				for(size_t i=0; i<arrayLen; i++){
					ret->itg(n,i) += this->itg(n*stride+j,i);
				}
			}
		}
	}
	else if(this->is_time()){
		ret = new DvObject(DvTime(), this->Dims(), newSeqLen);
		for(size_t n=0; n<newSeqLen; n++){
			for(size_t j=0; j<stride; j++){
				for(size_t i=0; i<arrayLen; i++){
					ret->time(n,i) += this->time(n*stride+j,i);
				}
			}
		}
	}
	else{
		ret = new DvObject();
		ret->error("Sum records not supported for data type");
		return ret;
	}
	
	*ret /= (int)stride;
	
	ret->copy_xrefs_from(*this);
	
	return ret;

}


DvObject_var DvObject::sliceRecEvery(size_t stride){

	DvObject_var ret;
	size_t newSeqLen = seqSize() / stride;
	
	if(this->is_dbl()) {
		ret = new DvObject(0.0, this->Dims(), newSeqLen);
		for(size_t n=0; n<newSeqLen; n++){
			for(size_t i=0; i<this->arraySize(); i++){
				ret->dbl(n,i) = this->dbl(n*stride,i);
			}
		}
	}
	else if(this->is_int()) {
		ret = new DvObject((int)0, this->Dims(), newSeqLen);
		for(size_t n=0; n<newSeqLen; n++){
			for(size_t i=0; i<this->arraySize(); i++){
				ret->itg(n,i) = this->itg(n*stride,i);
			}
		}
	}
	else if(this->is_time()){
		ret = new DvObject(DvTime(), this->Dims(), newSeqLen);
		for(size_t n=0; n<newSeqLen; n++){
			for(size_t i=0; i<this->arraySize(); i++){
				ret->time(n,i) = this->time(n*stride,i);
			}
		}
	}
	else{
		ret = new DvObject();
		ret->error("Slice records not supported for data type");
		return ret;
	}
	
	ret->copy_xrefs_from(*this);
	return ret;

}
  
  
DvEvent DvObject::getTimeInterval(DvString comp){
	
	if(is_time()){
		DvString recRange = comp.between("<", ">");
		if(recRange.empty()) return DvEvent(Tdata[0], Tdata[Tdata.size()-1]);
		
		DvString start = recRange.before(',');
		DvString end = recRange.after(',');
		int recS = start.toInt();
		int recE = end.toInt();
		if(recS < 0 || recE >= (int)Tdata.size()) return DvEvent(Tdata[0], Tdata[Tdata.size()-1]);
		
		return DvEvent(Tdata[recS], Tdata[recE]);
	}
	else if(is_event()){
		if(seqSize() == 1) return this->eventS(0);
		
		DvString eventAt = comp.between("<", ">");
		if(eventAt.empty()) return this->eventS(0);
		
		size_t n = 0;
		
		if( eventAt.contains(",") ) {
			DvString from = eventAt.before(",");
			DvString to = eventAt.after(",");
			DvEvent event = this->eventS(from.toInt());
			event = event.event_union( this->eventS(to.toInt()) );
			return event;
		}
		else{
			n = (size_t) eventAt.toInt();
			return this->eventS(n);
		}
		

	}
	else{
		DvObject_var tt = getDep0();
		if(tt->is_time() || tt->is_event()) return tt->getTimeInterval(comp);
	}
	
	return DvEvent(); // should not get here

}

void DvObject::ensureVectorXYZ(){
	// ensure this vector is in cartesian
	if( this->isVectorXYZ() ) return;
	
	DvString rep = getRep();
	if(rep.empty()){
		this->error("[ensureVectorXYZ] failed to get representation");
		return;
	}
	
	if(rep == "rtp") this->RTPtoXYZThis();
	else if(rep == "rlp") this->RLPtoXYZThis();
	else if(rep == "rpz") this->RPZtoXYZThis();

	return;
}

