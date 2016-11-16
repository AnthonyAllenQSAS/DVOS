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


// ---------------- Target can be Time Tags or Scalar Sequence -----------------


bool DvObject::isJoined(DvObject &target){

    DvObject_var tt_target = target.getDep0();
    if(tt_target.is_nil()) tt_target = target;
            
	DvObject_var thisDep0 = this->getDep0();
	
	if(thisDep0->get_id() == -1 || thisDep0->is_str()){ // is_str catches badly formed tt
	
		// obj has no Depend_0
		if(seqLen == 1 || seqLen == tt_target->seqLen || tt_target->seqLen == 1) return true; // allow operation without join
		else return false;
	}
	
	if(thisDep0->get_id() == tt_target->get_id()) return true;

    if( seqLen == tt_target->seqLen){

		// try fuzzy join
		for(size_t n=0; n<thisDep0->seqSize(); n++) {
			if( fabs(thisDep0->asDouble(n) - tt_target->asDouble(n) ) > DV_CLOSE_ENOUGH) 	return false;
		}
		return true;
	}
	return false;
	
}

bool DvObject::same(DvObject &target){
	
	// e.g. compare time tags or scalar tags directly
	
	if(this->get_id() == target.get_id()) return true;
	
	if( totalSize() == target.totalSize()){
		// try fuzzy join
		for(size_t n=0; n<this->totalSize(); n++) {
			if( fabs(this->asDouble(n) - target.asDouble(n) ) > DV_CLOSE_ENOUGH) 	return false;
		}
		return true;
	}
	return false;
	
}

bool DvObject::isRegular(double &spacing){

	// check if depend_0 is regularly spaced
	double err;
	double Gaptest;
   
	spacing = -1.;
    DvObject_var D0 = getDep0();
	
    if ( D0.is_nil()) D0 = this;
    
	int n_seq = D0->seqSize();

    if ( n_seq < 3 || D0->arraySize() > 1) return false;

  
	spacing = D0->asDouble(1) - D0->asDouble(0);
	
	if (spacing < DV_CLOSE_ENOUGH) return false; // too small for fuzzy equality limit
  
	for (int i=2; i<n_seq; i++){
    	Gaptest = D0->asDouble(i) - D0->asDouble(i-1);
    	err = fabs(spacing - Gaptest);        
    	if (err > DV_CLOSE_ENOUGH ) return false;
	}
	return true;

}

double DvObject::get_spacing(){
    double spacing, next;
	int repeat=0;
	
    if(this->is_nil()) return  0.;
	
	int n_seq = seqSize();

	if( n_seq == 0) return 0.;
	if( n_seq == 1) return this->asDouble(0)* 0.1; // Nominal value

	double spacing_sign = 1.0;
	
  
	next = asDouble(1) - asDouble(0);
	
	spacing = fabs(next);
	spacing_sign = next / spacing;
	
	for (int i=2; i<n_seq; i++){
    	next = fabs(asDouble(i) - asDouble(i-1));
    	if(next < spacing && next > 1.e-20) { // trap duplicate time tags, next < 1.e-20
			repeat++;
			if(repeat > 1){ // require shorter spacing more than once to avoid single glitch
				spacing = next; 
				repeat=0;
			}
		}
	}
    return spacing*spacing_sign;
}


DvObject_var DvObject::Join(DvObject_var &dobj, bool withXrefs){
    
    return this->Join(*dobj, withXrefs);
    
}

DvObject_var DvObject::Join(DvObject &dobj, bool withXrefs){

    if(isJoined(dobj)) return this;

    if(xref_exists(DV_JOIN_METHOD)){
        DvString method = this->getXrefText(DV_JOIN_METHOD);
        
        if(method == DV_NN) return nnJoin(dobj, withXrefs);
        
        else if(method == DV_BOXCAR) return boxcarJoin(dobj, withXrefs);
        
        return linearJoin(dobj, withXrefs);
    }
    
    // else use defaults
    
    double ratio = 1;

    if(dobj.seqSize() > 1 && this->seqSize() > 1){
        
        // find target spacing
        DvObject_var Target = dobj.getDep0();
        if(Target.is_nil()) Target = dobj;

        double targetSpacing = Target->get_spacing();
    
        // find this spacing
        DvObject_var d0 = getDep0();
        if(d0.is_nil()) d0 = this;
        double thisSpacing = d0->get_spacing();
        
        ratio = targetSpacing / thisSpacing;
    }
    
    // do join with sensible defaults
    
    if(this->is_dbl() || this->is_int() || this->is_time()) {
        
        // use boxcar if target cadence is significantly longer than data
        if (ratio > 2) return boxcarJoin(dobj, withXrefs);
        
        return linearJoin(dobj, withXrefs);
    }
    else {
        return nnJoin(dobj, withXrefs);
    }
    
}



DvObject_var DvObject::linearJoin(DvObject &dobj, bool withXrefs, DvMask *Gmsk){

    if(isJoined(dobj)) return this;

    bool doGaps = true;
    if(Gmsk != NULL) doGaps = false;
   
    // Decide if this or Depend_0
    DvObject_var Target = dobj.getDep0();
    if(Target.is_nil()) Target = dobj;

	// joins this onto Target
	
	if(Target.is_nil()){
        DvObject_var ret = new DvObject();
		ret->error("[LinearJoin]: no target Depend_0" );
		return ret;
	}

    // test if already joined (or join not needed)
    if( isJoined(*Target) ) return this;

    // get target DEPEND_0 as seq of doubles
    valarray <double> ssTarget = Target->getCentres();
	
	
	DvObject_var gap_opt_ptr = this->get_xref(DV_GAP_OPTION);
	DvString gap_opt;
	if(gap_opt_ptr.is_ok()){
		gap_opt = gap_opt_ptr->str();
	}
	else {
		gap_opt = DV_LINEAR;
	}

    DvObject_var end_opt_ptr = this->get_xref(DV_END_OPTION);
    DvString end_opt;
    if(end_opt_ptr.is_ok()){
        end_opt = end_opt_ptr->str();
    }
    else{
        end_opt = DV_NN;
        if(gap_opt == DV_REMOVE) end_opt = DV_REMOVE;
    }
    
	// get DEPEND_0 as a sequence of doubles
    
    valarray <double> ssThis = this->getDep0()->getCentres();

	// invert decreasing tags
	size_t nThis = ssThis.size();
	size_t nTarget = ssTarget.size();

	if(nThis < 2){
        DvObject_var ret = new DvObject();
        ret->error("[LinearJoin]: Data sequence too short for linear interp" );
        return ret;
	}
	
	if(ssThis[0] > ssThis[nThis-1]){
		if(ssTarget[0] > ssTarget[nTarget-1]){
			ssTarget *= -1.;
			ssThis  *= -1.;
		}
        else{
            DvObject_var ret = new DvObject();
            ret->error("[LinearJoin]: target Depend_0 decreasing, but data Depend_0 increasing" );
            return ret;
        }
	}
    if(ssTarget[0] > ssThis[nThis-1] || ssThis[0] > ssTarget[nTarget-1]){
        DvObject_var ret = new DvObject();
        ret->error("[LinearJoin]: Disjoint data sequences in join");
        return ret;
	}
	

	// gap tolerance in seconds
	DvObject_var gap_ptr = this->get_xref(DV_GAP_WIDTH);
	double gap;
	if(gap_ptr.is_ok())
		gap = gap_ptr->dbl();
	else{
		gap = Target->get_spacing();
		if(gap <=0) {
			gap = 4.; // default to cluster spin res
			error( "[LinearJoin]: Cannot determine target spacing" );
		}
		gap *=1.5; 
		
	}
	
    if(Gmsk == NULL) Gmsk = new DvMask(nTarget);
    
    // has only one element on first call from MultiJoin
    if(Gmsk->size() != nTarget){
		Gmsk->resize(nTarget);
	}
	
	// do join
		
	DvObject_var result = new DvObject(this, nTarget, false); // Array dim of this, seq length of target

	size_t i_target=0;
	size_t i_this=1; // ensures interp starts from 1st 2 elements

	if(ssTarget[0] < ssThis[0]){

        // data starts after target, get within gap tolerence of data
        while (ssTarget[i_target] + gap < ssThis[0] ) {
			// treat as gap
            if(end_opt == DV_NN) {
                (*result)[i_target] = (*this)[i_this-1];
            }
            else if(end_opt == DV_REMOVE) {
                (*Gmsk)[i_target]=false; // remove (done later)
            }
			else if(end_opt == DV_FILL) {
                if(is_dbl()) for(size_t n =0; n<arraySize(); n++) result->dbl(i_target, n) = DvNaN;
                else if(is_int()) for(size_t n =0; n<arraySize(); n++) result->itg(i_target, n) = 0;
                else if(is_str()) for(size_t n =0; n<arraySize(); n++) result->str(i_target, n) = "";
                else if(is_time()) for(size_t n =0; n<arraySize(); n++) result->time(i_target, n) = DvTime();
                else if(is_event()) for(size_t n =0; n<arraySize(); n++) result->event(i_target, n) = DvEvent();
			}
			else if(end_opt == DV_LINEAR) {
				double slope = (ssTarget[i_target] - ssThis[i_this-1]) / (ssThis[i_this] - ssThis[i_this-1]);
				(*result)[i_target] = (*this)[i_this];
				(*result)[i_target] -= (*this)[i_this-1] ;
				(*result)[i_target] *= slope;
				(*result)[i_target] += (*this)[i_this-1];

			}
            else {
                 // zero fill
                (*result)[i_target] = 0;
            }
			i_target++;
		}
        
		if(ssTarget[i_target] < ssThis[0]){
           
			// extrapolate (same call as interpolate)
			double slope = (ssTarget[i_target] - ssThis[i_this-1]) / (ssThis[i_this] - ssThis[i_this-1]);
			(*result)[i_target]  = (*this)[i_this];
			(*result)[i_target] -= (*this)[i_this-1];
			(*result)[i_target] *= slope;
			(*result)[i_target] += (*this)[i_this-1];
			i_target++;
		}
		// else fall through to interpolation
		
	}
		
	// ssThis[i_this] starts above ssTarget[i_target] 
	while(i_target < ssTarget.size()){

        double target = ssTarget[i_target];
		
		// find data pt above target tag
		while(ssThis[i_this] < target  &&  i_this < ssThis.size()-1 ) i_this++;

		// is this a data gap or off the end of data by more than gap?
        if( ( ssThis[i_this] - ssThis[i_this-1] > gap ) ){
           
            // data gap
            
            if(gap_opt == DV_FILL) {
                if(is_dbl()) for(size_t n =0; n<arraySize(); n++) result->dbl(i_target, n) = DvNaN;
                else if(is_int()) for(size_t n =0; n<arraySize(); n++) result->itg(i_target, n) = 0;
                else if(is_str()) for(size_t n =0; n<arraySize(); n++) result->str(i_target, n) = "";
                else if(is_time()) for(size_t n =0; n<arraySize(); n++) result->time(i_target, n) = DvTime();
                else if(is_event()) for(size_t n =0; n<arraySize(); n++) result->event(i_target, n) = DvEvent();
            }
            else if(gap_opt == DV_LINEAR) {
                
                // uses DvRecord operators
                double slope = (target - ssThis[i_this-1]) / (ssThis[i_this] - ssThis[i_this-1]);
                (*result)[i_target]  = (*this)[i_this];
                (*result)[i_target] -= (*this)[i_this-1];
                (*result)[i_target] *= slope;
                (*result)[i_target] += (*this)[i_this-1];
            }
            else if(gap_opt == DV_NN) {
                if( (ssThis[i_this] - target) < (target - ssThis[i_this-1]) ) (*result)[i_target] = (*this)[i_this];
                else (*result)[i_target] = (*this)[i_this-1];
            }
            else if(gap_opt == DV_REMOVE) {
                (*Gmsk)[i_target]=false; 
            }
            else  {
                // zero fill
                (*result)[i_target] = 0;
            }
            
        }
        else if( ( target - ssThis[i_this] > gap) ){
            
            // off end
            
            if(end_opt == DV_NN) {
                if( (ssThis[i_this] - target) < (target - ssThis[i_this-1]) ) (*result)[i_target] = (*this)[i_this];
                else (*result)[i_target] = (*this)[i_this-1];
            }
            else if(end_opt == DV_REMOVE) {
                (*Gmsk)[i_target]=false;
            }
            else if(end_opt == DV_FILL) {
                if(is_dbl()) for(size_t n =0; n<arraySize(); n++) result->dbl(i_target, n) = DvNaN;
                else if(is_int()) for(size_t n =0; n<arraySize(); n++) result->itg(i_target, n) = 0;
                else if(is_str()) for(size_t n =0; n<arraySize(); n++) result->str(i_target, n) = "";
                else if(is_time()) for(size_t n =0; n<arraySize(); n++) result->time(i_target, n) = DvTime();
                else if(is_event()) for(size_t n =0; n<arraySize(); n++) result->event(i_target, n) = DvEvent();
            }
            else if(end_opt == DV_LINEAR) {
                
                // uses DvRecord operators
                double slope = (target - ssThis[i_this-1]) / (ssThis[i_this] - ssThis[i_this-1]);
                (*result)[i_target]  = (*this)[i_this];
                (*result)[i_target] -= (*this)[i_this-1];
                (*result)[i_target] *= slope;
                (*result)[i_target] += (*this)[i_this-1];
            }
            else  {
                // zero fill
                (*result)[i_target] = 0;
            }
            
            
        }
		else{
			// Not a gap, normal interpolation
			
			// if we are off end of data then interpolation algorithm works as extrapolation
		
			double slope = (target - ssThis[i_this-1]) / (ssThis[i_this] - ssThis[i_this-1]);
			(*result)[i_target]  = (*this)[i_this];
			(*result)[i_target] -= (*this)[i_this-1];
			(*result)[i_target] *= slope;
			(*result)[i_target] += (*this)[i_this-1];
		}
		i_target++;
	}

	// keep these for future use (xrefs and propagation in operator chains)
		
	result->change_xref(DV_GAP_WIDTH, gap);
	result->change_xref(DV_JOIN_METHOD, DV_LINEAR); // so double xrefs use same
    result->change_xref(DV_GAP_OPTION, gap_opt);
    result->change_xref(DV_END_OPTION, end_opt);
	
    this->delete_xref(DV_GAP_WIDTH);
    this->delete_xref(DV_BOX_WIDTH);
    this->delete_xref(DV_JOIN_METHOD);
    this->delete_xref(DV_GAP_OPTION);
    this->delete_xref(DV_END_OPTION);

    // attach xrefs from this to result after joining to target
    if(withXrefs) result->joinXrefs(*this, *Target);
   
    result->change_xref(DEPEND_0, Target);

    if(doGaps && gap_opt == DV_REMOVE && Gmsk->size() == Target->seqSize()){
        
        // do gap mask recursively now as this was not called by MultiJoin.
        result->apply_mask(*Gmsk); // recursively applies to xrefs if attached
        delete Gmsk;
    }

	return result;

}


// ---------- BOXCAR ----------


DvObject_var DvObject::boxcarJoin(DvObject &dobj, bool withXrefs, DvMask *Gmsk){

    // set up minimum number of point in boxcar
    int minBoxcarNum = 3;
    DvObject_var min_box_ptr = this->get_xref(DV_MIN_BOXCAR);
    if(min_box_ptr.is_ok()){
        minBoxcarNum = min_box_ptr->asInt();
    }
    if(minBoxcarNum < 1) minBoxcarNum = 1; // must have at least one point in box
    
    bool doGaps = true;
    if(Gmsk != NULL) doGaps = false;

    // Decide if this or Depend_0
    DvObject_var Target = dobj.getDep0();
    if(Target.is_nil()) Target = dobj;
    
	// do join
	
	// get target DEPEND_0 as seq of doubles
    valarray <double> ssTarget = Target->getCentres();
	
    if(ssTarget.size() < 1){
        DvObject_var ret = new DvObject();
        ret->error("[BoxcarJoin]: no target Depend_0" );
        return ret;
    }

	// Boxcar does not test if already joined since target tt may remain the same

	
	DvObject_var gap_opt_ptr = this->get_xref(DV_GAP_OPTION);
	DvString gap_opt;
	if(gap_opt_ptr.is_ok()){
		gap_opt = gap_opt_ptr->str();
	}
	else {
		gap_opt = DV_LINEAR;
	}
    
    DvObject_var end_opt_ptr = this->get_xref(DV_END_OPTION);
    DvString end_opt;
    if(end_opt_ptr.is_ok()){
        end_opt = end_opt_ptr->str();
    }
    else{
        end_opt = DV_NN;
        if(gap_opt == DV_REMOVE) end_opt = DV_REMOVE;
    }
    
	
	// get DEPEND_0 as a sequence of doubles
		
    valarray <double> ssThis = this->getDep0()->getCentres();

    
	// invert decreasing tags
	size_t nThis = ssThis.size();
	size_t nTarget = ssTarget.size();
	
	if(ssThis[0] > ssThis[nThis-1]){
		if(ssTarget[0] > ssTarget[nTarget-1]){
			ssTarget *= -1.;
			ssThis  *= -1.;
		}
        else{
            DvObject_var ret = new DvObject();
            ret->error("[BoxcarJoin]: target Depend_0 decreasing, but data Depend_0 increasing" );
            return ret;
		}
	}
	
    if(ssTarget[0] > ssThis[nThis-1] || ssThis[0] > ssTarget[nTarget-1]){
        DvObject_var ret = new DvObject();
        ret->error("[BoxcarJoin]: Disjoint data sequences in join");
        return ret;
	}
	
	// boxcar in seconds
	DvObject_var width_ptr = this->get_xref(DV_BOX_WIDTH);
	double width;

	if(width_ptr.is_ok() )
		width = width_ptr->asDouble();
	else{
		width = 2 * Target->get_spacing();
		if(width <=0)  error( "[BoxcarJoin]: Cannot determine target spacing" );
	
		width = 10.;
	}

	double halfwidth = width*0.5;

    DvMask *LinMsk;
    if( gap_opt == DV_LINEAR || end_opt == DV_LINEAR) LinMsk = new DvMask(nTarget);
    if(Gmsk == NULL) Gmsk = new DvMask(nTarget);
    
    // has only one element on first call from MultiJoin
    if(Gmsk->size() != nTarget){
        Gmsk->resize(nTarget);
    }

	
	// size of boxcar record
	size_t nArray = arraySize();
	
	// do join
	DvObject_var result = new DvObject(this, nTarget, false); // Array dim of this, seq length of target

	size_t i_target=0;
	size_t i_this=1; // ensures interp starts from 1st 2 elements
	
	if(ssTarget[0] < ssThis[0]){
		// data starts after target, get within gap tolerence of data

        while (ssTarget[i_target] + halfwidth < ssThis[0] ) {
			// treat as gap
			if(end_opt == DV_FILL) {
				(*result)[i_target] = DvNaN;
			}
            else if(end_opt == DV_LINEAR) {
                (*result)[i_target] = 0.;
               (*LinMsk)[i_target]=false; // interpolate later
			}
			else if(end_opt == DV_NN) {
				(*result)[i_target] = (*this)[i_this-1];
			}
			else if(end_opt == DV_REMOVE) {
				(*Gmsk)[i_target]=false; // remove (done later)
			}
			else (*result)[i_target] = 0; // zero fill
			i_target++;
		}
		
	}
		
	// ssThis[i_this] starts above ssTarget[i_target] 
	while(i_target < ssTarget.size()){
	
		double target = ssTarget[i_target];
		
		// find first data pt in boxcar (it will be i_this - 1)
        while(ssThis[i_this-1] < target - halfwidth &&  i_this < nThis-1 ) {
            i_this++;
        }
		// accumulate boxcar here
		valarray <double> boxcar(0.0, nArray);	
		size_t count= 0;	
		
		size_t n = i_this-1; // first point in boxcar
		while(n < nThis && ssThis[n] < target + halfwidth){
			slice sl(n*nArray, nArray, 1);
			valarray <double> data = this->Ddata[sl];
			for(size_t k=0; k<nArray; k++){
				boxcar[k] += data[k];
			}
			count ++;
			n++;
		}
		if(n > i_this - 1 && count >= minBoxcarNum) {
			for(size_t k=0; k<nArray; k++){
				boxcar[k] /= count;
			}
			(*result)[i_target] = boxcar;
		}
        else if(i_this == nThis){
            
            // off end (at most 2 data points
            if(gap_opt == DV_FILL) {
                (*result)[i_target] = DvNaN;
            }
            else if(gap_opt == DV_LINEAR) {
                (*result)[i_target] = 0.;
                (*LinMsk)[i_target]=false; // interpolate later
            }
            else if(gap_opt == DV_NN) {
                (*result)[i_target] = (*this)[i_this-1];
            }
            else if(gap_opt == DV_REMOVE) {
                (*result)[i_target] = 0.; // remove (done later)
                (*Gmsk)[i_target]=false; 
            }
            else (*result)[i_target] = 0.; // zero fill record
            
        }
        else if( ( target - ssThis[i_this] > width) ){
            
            // off end
            
            if(end_opt == DV_NN) {
                if( (ssThis[i_this] - target) < (target - ssThis[i_this-1]) ) (*result)[i_target] = (*this)[i_this];
                else (*result)[i_target] = (*this)[i_this-1];
            }
            else if(end_opt == DV_REMOVE) {
                (*Gmsk)[i_target]=false;
            }
            else if(end_opt == DV_FILL) {
                if(is_dbl()) for(size_t n =0; n<arraySize(); n++) result->dbl(i_target, n) = DvNaN;
                else if(is_int()) for(size_t n =0; n<arraySize(); n++) result->itg(i_target, n) = 0;
                else if(is_str()) for(size_t n =0; n<arraySize(); n++) result->str(i_target, n) = "";
                else if(is_time()) for(size_t n =0; n<arraySize(); n++) result->time(i_target, n) = DvTime();
                else if(is_event()) for(size_t n =0; n<arraySize(); n++) result->event(i_target, n) = DvEvent();
            }
            else if(end_opt == DV_LINEAR) {
                
                // uses DvRecord operators
                double slope = (target - ssThis[i_this-1]) / (ssThis[i_this] - ssThis[i_this-1]);
                (*result)[i_target]  = (*this)[i_this];
                (*result)[i_target] -= (*this)[i_this-1];
                (*result)[i_target] *= slope;
                (*result)[i_target] += (*this)[i_this-1];
            }
            else  {
                // zero fill
                (*result)[i_target] = 0;
            }
            
            
        }
        else{
            // data gap
            if(gap_opt == DV_FILL) {
                (*result)[i_target] = DvNaN;
            }
            else if(gap_opt == DV_LINEAR) {
                (*result)[i_target] = 0.;
                (*LinMsk)[i_target]=false; // interpolate later

            }
            else if(gap_opt == DV_NN) {
                if( n != nThis && (ssThis[n] - target) < (target - ssThis[n-1]) ) (*result)[i_target] = (*this)[n];
                else (*result)[i_target] = (*this)[n-1];
            }
            else if(gap_opt == DV_REMOVE) {
                (*result)[i_target] = 0.; // remove (done later)
                (*Gmsk)[i_target]=false; 
            }
            else (*result)[i_target] = 0.; // zero fill record
            
        }

		i_target++;
	}
	
	// this object has no xrefs yet and has not had the gap mask applied	

	// keep these for future use (xrefs and propagation in operator chains)
	result->change_xref(DV_BOX_WIDTH, width);
    result->change_xref(DV_GAP_OPTION, gap_opt);
    result->change_xref(DV_END_OPTION, end_opt);
	result->change_xref(DV_JOIN_METHOD, DV_BOXCAR); // so double xrefs use same

    this->delete_xref(DV_BOX_WIDTH);
    this->delete_xref(DV_GAP_WIDTH);
    this->delete_xref(DV_JOIN_METHOD);
    this->delete_xref(DV_GAP_OPTION);
    this->delete_xref(DV_END_OPTION);

    // attach xrefs from this to result after joining to target
    if(withXrefs) result->joinXrefs(*this, *Target); // true means use boxcar
    
    if( gap_opt == DV_LINEAR || end_opt == DV_LINEAR) {
        result->applyLinearFill(*LinMsk);
        delete LinMsk;
    }
    
    result->change_xref(DEPEND_0, Target);

    if( doGaps && gap_opt == DV_REMOVE && Gmsk->size() == Target->seqSize()){
        
        // do gap mask recursively now as this was not called by MultiJoin.
        result->apply_mask(*Gmsk); // recursively applies to xrefs
        delete Gmsk;
    }
    

	return result;
	
}

// ---------- Nearest Neighbour ----------


DvObject_var DvObject::nnJoin(DvObject &dobj, bool withXrefs, DvMask *Gmsk){

    bool doGaps = true;
    if(Gmsk != NULL) doGaps = false;
    
    // Decide if this or Depend_0
    DvObject_var Target = dobj.getDep0();
    if(Target.is_nil()) Target = dobj;

	// Nearest Neighbour
	
	// get target DEPEND_0 as seq of doubles
    valarray <double> ssTarget = Target->getCentres();

    if(ssTarget.size() < 1){
        DvObject_var ret = new DvObject();
        ret->error("[nnJoin]: no target Depend_0" );
        return ret;
	}

	// test if already joined (or join not needed)
	if( isJoined(*Target) ) return this;


	DvObject_var gap_opt_ptr = this->get_xref(DV_GAP_OPTION);
	DvString gap_opt;
	if(gap_opt_ptr.is_ok()){
		gap_opt = gap_opt_ptr->str();
	}
	else {
		gap_opt = DV_NN;
	}
    
    DvObject_var end_opt_ptr = this->get_xref(DV_END_OPTION);
    DvString end_opt;
    if(end_opt_ptr.is_ok()){
        end_opt = end_opt_ptr->str();
    }
    else{
        end_opt = DV_NN;
        if(gap_opt == DV_REMOVE) end_opt = DV_REMOVE;
    }
	
	if(gap_opt == DV_LINEAR && !(this->is_dbl()) ) gap_opt = DV_NN; 
	
	// get DEPEND_0 as a sequence of doubles
		
    valarray <double> ssThis = this->getDep0()->getCentres();

	// invert decreasing tags
	size_t nThis = ssThis.size();
	size_t nTarget = ssTarget.size();
	
	if(ssThis[0] > ssThis[nThis-1]){
		if(ssTarget[0] > ssTarget[nTarget-1]){
			ssTarget *= -1.;
			ssThis  *= -1.;
		}
        else{
            DvObject_var ret = new DvObject();
            ret->error("[nnJoin]: target Depend_0 decreasing, but data Depend_0 increasing" );
            return ret;
		}
	}
	
    if(ssTarget[0] > ssThis[nThis-1] || ssThis[0] > ssTarget[nTarget-1]){
        DvObject_var ret = new DvObject();
        ret->error("[nnJoin]: Disjoint data sequences in join");
        return ret;
	}
	
	// gap tolerance in seconds
	DvObject_var gap_ptr = this->get_xref(DV_GAP_WIDTH);
	double gap;
	if(gap_ptr.is_ok())
		gap = gap_ptr->dbl();
	else{
		gap = Target->get_spacing();
		if(gap <=0) error( "[LinearJoin]: Cannot determine target spacing" );
	
		gap *=1.5; 
		
	}
	
    
    if(Gmsk == NULL) Gmsk = new DvMask(nTarget);
    
    // has only one element on first call from MultiJoin
    if(Gmsk->size() != nTarget){
        Gmsk->resize(nTarget);
    }
    
	
	// do join
		
	DvObject_var result = new DvObject(this, nTarget, false); // Array dim of this, seq length of target

	size_t i_target=0;
	size_t i_this=1; // ensures interp starts from 1st 2 elements
	
	if(ssTarget[0] < ssThis[0]){
		// data starts after target, get within gap tolerence of data
		while (ssTarget[i_target] + gap < ssThis[0] ) {
            if(end_opt == DV_FILL) {
                if(is_dbl()) for(size_t n =0; n<arraySize(); n++) result->dbl(i_target, n) = DvNaN;
                else if(is_int()) for(size_t n =0; n<arraySize(); n++) result->itg(i_target, n) = 0;
                else if(is_str()) for(size_t n =0; n<arraySize(); n++) result->str(i_target, n) = "";
                else if(is_time()) for(size_t n =0; n<arraySize(); n++) result->time(i_target, n) = DvTime();
                else if(is_event()) for(size_t n =0; n<arraySize(); n++) result->event(i_target, n) = DvEvent();
			}
			else if(end_opt == DV_LINEAR ) {
				double slope = (ssTarget[i_target] - ssThis[i_this-1]) / (ssThis[i_this] - ssThis[i_this-1]);
				(*result)[i_target] = (*this)[i_this];
				(*result)[i_target] -= (*this)[i_this-1] ;
				(*result)[i_target] *= slope;
				(*result)[i_target] += (*this)[i_this-1];
			}
			
			else if(end_opt == DV_REMOVE) {
				(*Gmsk)[i_target]=false; // remove (done later)
			}
			else if(end_opt == DV_ZERO_FILL) {
				(*result)[i_target] = 0; // zero fill
			}
			else  {
				// NN default
				(*result)[i_target] = (*this)[i_this-1];
			}
			i_target++;
		}
		
		if(ssTarget[i_target] < ssThis[0]){
			
			(*result)[i_target]  = (*this)[0];
			i_target++;
		}
		// else fall through to interpolation
		
	}
		
	// ssThis[i_this] starts above ssTarget[i_target] 
	while(i_target < ssTarget.size()){
	
		double target = ssTarget[i_target];
		
		while(ssThis[i_this] < target  &&  i_this<ssThis.size()-1 ) i_this++;
		
		// is this a data gap or off the end of data by more than gap?
		
        if( ( ssThis[i_this] - ssThis[i_this-1] > gap ) ){
            
            // data gap
            
            if(gap_opt == DV_FILL) {
                if(is_dbl()) for(size_t n =0; n<arraySize(); n++) result->dbl(i_target, n) = DvNaN;
                else if(is_int()) for(size_t n =0; n<arraySize(); n++) result->itg(i_target, n) = 0;
                else if(is_str()) for(size_t n =0; n<arraySize(); n++) result->str(i_target, n) = "";
                else if(is_time()) for(size_t n =0; n<arraySize(); n++) result->time(i_target, n) = DvTime();
                else if(is_event()) for(size_t n =0; n<arraySize(); n++) result->event(i_target, n) = DvEvent();           
            }
            else if(gap_opt == DV_LINEAR ) {
                // uses DvRecord operators
                double slope = (target - ssThis[i_this-1]) / (ssThis[i_this] - ssThis[i_this-1]);
                (*result)[i_target]  = (*this)[i_this];
                (*result)[i_target] -= (*this)[i_this-1];
                (*result)[i_target] *= slope;
                (*result)[i_target] += (*this)[i_this-1];
            }
            else if(gap_opt == DV_NN) {
                if( (ssThis[i_this] - target) < (target - ssThis[i_this-1]) ) (*result)[i_target] = (*this)[i_this];
                else (*result)[i_target] = (*this)[i_this-1];
            }
            else if(gap_opt == DV_REMOVE) {
                (*Gmsk)[i_target]=false; 
            }
            else (*result)[i_target] = 0; // zero fill record
            
            
        }
        else if( ( target - ssThis[i_this] > gap) ){
            
            // off end
            
            if(end_opt == DV_FILL) {
                if(is_dbl()) for(size_t n =0; n<arraySize(); n++) result->dbl(i_target, n) = DvNaN;
                else if(is_int()) for(size_t n =0; n<arraySize(); n++) result->itg(i_target, n) = 0;
                else if(is_str()) for(size_t n =0; n<arraySize(); n++) result->str(i_target, n) = "";
                else if(is_time()) for(size_t n =0; n<arraySize(); n++) result->time(i_target, n) = DvTime();
                else if(is_event()) for(size_t n =0; n<arraySize(); n++) result->event(i_target, n) = DvEvent();
            }
            else if(end_opt == DV_LINEAR ) {
                // uses DvRecord operators
                double slope = (target - ssThis[i_this-1]) / (ssThis[i_this] - ssThis[i_this-1]);
                (*result)[i_target]  = (*this)[i_this];
                (*result)[i_target] -= (*this)[i_this-1];
                (*result)[i_target] *= slope;
                (*result)[i_target] += (*this)[i_this-1];
            }
            else if(end_opt == DV_NN) {
                (*result)[i_target] = (*this)[i_this];
            }
            else if(end_opt == DV_REMOVE) {
                (*Gmsk)[i_target]=false; 
            }
            else (*result)[i_target] = 0; // zero fill record
            
            
        }
		else{
			// Not a gap, nn interpolation
			if( (ssThis[i_this] - target) < (target - ssThis[i_this-1]) ) (*result)[i_target] = (*this)[i_this];
			else (*result)[i_target] = (*this)[i_this-1];
		}
			
		i_target++;
	}
	
	// this object has no xrefs yet and has not had the gap mask applied	

	// keep these for future use (xrefs and propagation in operator chains)
	result->change_xref(DV_GAP_WIDTH, gap);
    result->change_xref(DV_GAP_OPTION, gap_opt);
    result->change_xref(DV_END_OPTION, end_opt);
	result->change_xref(DV_JOIN_METHOD, DV_BOXCAR); // so double xrefs use same

    this->delete_xref(DV_BOX_WIDTH);
    this->delete_xref(DV_GAP_WIDTH);
    this->delete_xref(DV_JOIN_METHOD);
    this->delete_xref(DV_GAP_OPTION);
    this->delete_xref(DV_END_OPTION);

    // attach xrefs from this to result after joining to target
    if(withXrefs) result->joinXrefs(*this, *Target); // true means use boxcar
    
    if(doGaps && gap_opt == DV_REMOVE && Gmsk->size() == Target->seqSize()){
        
        // do gap mask recursively now as this was not called by MultiJoin.
        result->apply_mask(*Gmsk); // recursively applies to xrefs
        delete Gmsk;
    }

    result->change_xref(DEPEND_0, Target);

	return result;
	
}



// ---------------------

void DvObject::joinXrefs(DvObject &obj, DvObject &Target){

	// xrefs will come from the pre-joined object (obj) (we need their original time tags)
	// and be attached to the joined object

    DvNode *next = obj.first_xref();
	while(next)
	{
		if(next->name() == DEPEND_0 ) {
			next = next->next;
			continue;
		}
		  
		if( next->obj()->seqSize() != obj.seqSize() || !next->obj()->is_dbl()) {
			// record varying double type
			this->change_xref(next->name(), next->obj()); ;
			next = next->next;
			continue;
		}
		
		// join rec varying double data xrefs
		DvObject_var xref = next->obj();
		
		DvObject_var D0 = xref->getDep0();
		if( D0.is_nil() ) {
			// sequence length is same as object, so use tt from object
			D0 = obj.getDep0();
			xref = new DvObject(xref); // take copy
			xref->change_xref(DEPEND_0, D0); 
		}

		DvObject_var gapWidth = this->get_xref(DV_GAP_WIDTH);
		xref->change_xref(DV_GAP_WIDTH, gapWidth);
		
		DvObject_var boxWidth = this->get_xref(DV_BOX_WIDTH);
		xref->change_xref(DV_BOX_WIDTH, boxWidth);
		
		DvObject_var xref_result;

		DvString join_method = this->getXrefText(DV_JOIN_METHOD);
		
		if(join_method == DV_BOXCAR){
			xref_result = xref->boxcarJoin(Target);
		}
		else if(join_method == DV_NN){
			xref_result = xref->nnJoin(Target);
		}
		else{
			xref_result = xref->linearJoin(Target);
		}
        
		this->change_xref(next->name(), xref_result);

		next = next->next;
	} 
}

// ------------------

DvObject_var DvObject::interpAt(DvTime &target){
	size_t iAbove = 0;
	DvObject_var tt = getTimeTags();
	if( tt.is_nil() || tt->not_time() ) return new DvObject();
	
	for (size_t i=0; i<seqSize(); i++) {
			iAbove = i;
			if( tt->time(i) > target ) break;
	}
	if(iAbove == 0) iAbove = 1; // extrapolate
	
	valarray<double> values(arraySize()); 
	double h1 = tt->time(iAbove) - target;
	double h2 = target - tt->time(iAbove-1);
	double h3 = h1 + h2;
	for(size_t k=0; k<arraySize(); k++){
		values[k] = ( this->asDouble(iAbove-1,k)*h1  + this->asDouble(iAbove, k)*h2  ) / h3;
	}

	DvObject_var ret = new DvObject(values, this->dims, 1);
	
	ret->copy_xref_from(SI_CONVERSION, *this);
	ret->copy_xref_from(UNITS, *this);
	ret->copy_xref_from(SCALETYP, *this);
	ret->copy_xref_from(FRAME, *this);
	ret->copy_xref_from(TOFRAME, *this);
	ret->copy_xref_from(DEPEND_1, *this);
	ret->copy_xref_from(DEPEND_2, *this);
	ret->copy_xref_from(DEPEND_3, *this);
	ret->copy_xref_from(DEPEND_4, *this);
	ret->copy_xref_from(DELTA_PLUS, *this);
	ret->copy_xref_from(DELTA_MINUS, *this);
	ret->copy_xref_from(TENSOR_ORDER, *this);
	ret->copy_xref_from(COORDINATE_SYSTEM, *this);
	ret->copy_xref_from(REPRESENTATION, *this);
	ret->copy_xref_from(REPRESENTATION_1, *this);
	ret->copy_xref_from(REPRESENTATION_2, *this);
	ret->copy_xref_from(REPRESENTATION_3, *this);

	return ret;
}


valarray <double>  DvObject::getCentres(){
    // returns centre values for DEPEND_0
    DvObject_var deltaMinus = get_xref(DELTA_MINUS);
    DvObject_var deltaPlus = get_xref(DELTA_PLUS);

    // copy without xrefs, stops loads of tests
    valarray <double> centre = this->toDouble();
    
    valarray <double> offset;
    
    if(deltaPlus.is_ok())  {
        offset.resize(deltaPlus->totalSize(), 0.0);
        offset = deltaPlus->changeUnitsTo("1>s", "s")->toDouble();
    }
    
    if(deltaMinus.is_ok()) {
        valarray <double> dm;
        dm.resize(deltaMinus->totalSize(), 0.0);
        dm = deltaMinus->changeUnitsTo("1>s", "s")->toDouble();
        
        if(offset.size() == dm.size()) {
            offset -= dm;
        }
        else if(dm.size() == 1 && offset.size() > 0) {
            offset -= dm[0];
        }
        else if(offset.size() == 1) {
            offset.resize(dm.size(), offset[0]);
            offset -= dm;
        }
        else {
            offset.resize(dm.size(), 0.0);
            offset = dm;
            offset *= -1.;
        }
    }
    
    if(offset.size() > 0){
        offset *= 0.5;
        if(offset.size() == centre.size()) {
            centre +=  offset;
        }
        else if(offset.size() == 1) {
            centre +=  offset[0];
        }
    }   
 
    return centre;
    
}

void DvObject::applyLinearFill(DvMask &msk){
    
    // used for gap filling in Boxcar as it linearly interpolates between
    // boxcar averaged values rather than input values
    
    if(msk.size() != seqSize()) {
        error("[Mask] wrong length");
        // mask wrong length
        return;
    }
    
    if( msk.resultSize() < 2){
        error("[Mask] not enough valid values");
        return;
    }
    
    valarray <double> ssThis = this->getDep0()->toDouble();

    
    if(is_dbl() || is_int()){

        for(size_t i=0; i<seqSize(); i++){
            // at least two valid values exist somewhere
            if(msk[i] == false) {
                // get value below
                size_t i1 = i-1;
                if(i == 0) {
                    // no value below
                    i1 = i;
                    while(i1 < seqSize() && msk[i1] == 0) i1++;

                }
                // find higher value value
                size_t i2 = i1+1;
                while(i2 < seqSize() && msk[i2] == 0) i2++;
                if(i2 == seqSize()){
                    i2 = i1;
                    i1--; // this exists by now
                }
                
                
                double slope = (ssThis[i] - ssThis[i1]) / (ssThis[i2] - ssThis[i1]);
                (*this)[i]  = (*this)[i2];
                (*this)[i] -= (*this)[i1];
                (*this)[i] *= slope;
                (*this)[i] += (*this)[i1];
            }
        }
    }
 
    
}


// ----------------------

DvString DvJoinList::MultiJoin( DvObject_var    &target,
								DvJoinList     	&retList)
{
	// join method, gap method and gap/width are attached to objects as xefs
		
	DvMask Gmsk(1); // resized by join if used
	
	DvNode *node = this->first();

    while( node )
    {  
        if (node->obj()->seqSize() == 1)
        {
			// replace single input with sequence of right length and all values the same
			DvObject_var ret = new DvObject( *(node->obj()), target->seqSize());
	     	ret->change_xref(DEPEND_0, target) ; // is automatically joined
			
			retList.append(node->name(), ret);
			node = node->next;
			continue;
        }
		
		DvString join_method = node->obj()->getXrefText(DV_JOIN_METHOD);
		
		if(join_method == DV_BOXCAR){
			
			DvObject_var ret = node->obj()->boxcarJoin(*target, true, &Gmsk);
			retList.append(node->name(), ret);
		}
		else if(join_method == DV_NN){
			
			DvObject_var ret = node->obj()->nnJoin(*target, true, &Gmsk);
			retList.append(node->name(), ret);
		}
		else{
			// DV_LINEAR is default
			DvObject_var ret = node->obj()->linearJoin(*target, true, &Gmsk);
			retList.append(node->name(), ret);
		}
		node = node->next;
   }
    
    // Apply gap mask to all objects
	if(Gmsk.size() == target->seqSize() ){
		DvNode *next = retList.first();
		while(next){
			// do gap mask recursively now.
			next->obj()->apply_mask(Gmsk); // recursively applies to xrefs
			next = next->next;
		}

	}
	
	return DvString("Multi Join OK");
}



