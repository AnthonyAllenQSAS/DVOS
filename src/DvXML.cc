//
//  DvObj.h
//
//  Dvos -  Data variable object classes
//
//  Tony Allen
//  May 2014
//

#include "DvObj.h"
#include <iostream>
#include <fstream> 

unsigned int XMLblock = 32000;

char * readXML(DvString fullName){
	
	if(fullName.empty()) return NULL;

	std::ifstream ifp;
	unsigned int length;
		
	// clear eof flag or any errors
	ifp.clear();
	
	// try reading to see if file of that name exists 
	// use binary to stop Windows msys trying to remove \r and messing up
	ifp.open(fullName.c_str(), std::fstream::in | std::fstream::binary);

	if(!ifp.is_open() ) return NULL;
	
    // get length of file:
    ifp.seekg (0, ifp.end);
    length = ifp.tellg();
    ifp.seekg (0, ifp.beg);

	char * xml = new char[length+1];
	// NULL terminate for parsing
	xml[length] = '\0'; // done now as length changes
			
    // read data as blocks (read() uses a signed int which can be too small for MMS):
	if(length > XMLblock) {
		char * ptr = xml;
		while(length > XMLblock){
			ifp.read(ptr, XMLblock);
			ptr += XMLblock;
			length -= XMLblock;
		}
		ifp.read(ptr, length); // finish off
	}
	else {
		ifp.read(xml, length);
	}

	return xml;
}



char* getNextTag(char *xml, size_t &start, DvString &tag){
	// returns pointer to xml beyond this tag pair
	// start is the offset to tag content (NULL terminator inserted)
	// tag is name of tag
	start = 0;
	size_t end = 0;
	tag = "";
	
	while(*(xml+start) != '<') {
		if(*(xml+start) == '\0') return NULL; // end of input
		start++;
	}
	end = ++start; // move past '<'

	while(*(xml+end) != '>') {
		if(*(xml+end) == '\0') return NULL;
		end++;
	}
	tag = DvString(xml+start, end-start);
	start = ++end;
	
	// get closing tag

	DvString close = "</";
	close += tag;
	close += ">";
	
	DvString open = "<";
	open += tag;
	open += ">";
	
	char *ptr = strstr(xml+start, close.c_str());
	char *tagptr = strstr(xml+start, open.c_str());
	
	while(tagptr != NULL && ptr != NULL && tagptr < ptr  ){
		tagptr = strstr(tagptr+open.size(), open.c_str());
		ptr = strstr(ptr+close.size(), close.c_str());
	}

	if(ptr == NULL) {
		// no closing tag
		tag = "";
		return NULL;
	}
	
	
	// NULL terminate sub-str and return ptr to beyond it
	*ptr = '\0';
	return ptr+close.length(); // points to string beyond tag
}

// ----------------- DvList -----------------

void DvList::toXML(std::ofstream &ofp, DvString ListName) 
{
	if(size() == 0) return;
	
	ofp << "  <" << ListName<< ">\n";
	
		DvNode *node = first();
		while(node){
			node->toXML(ofp);
			node = node->next;
		}
	ofp << "  </" << ListName<< ">\n";
}

void DvList::fromXML(char *xml) 
{
	size_t start;
	DvString tag;
	char *ptr = getNextTag(xml, start, tag);
	while(ptr){
		if(tag == "DV_NODE") {
			DvNode *node = makeNode();
			node->fromXML(xml+start);
		}
		xml = ptr;
		ptr = getNextTag(xml, start, tag);
	}
}

// -------------------- DvNode ---------------------------

void DvNode::toXML(std::ofstream &ofp) 
{
	ofp << "  <DV_NODE>\n";
		ofp << "   <DV_NODE_NAME>";
		ofp << _name;
		ofp << "</DV_NODE_NAME>\n";

		_obj->toXML(ofp);

	ofp << "  </DV_NODE>\n";
}

void DvNode::fromXML(char *xml) 
{
	size_t start;
	DvString tag;
	
	char *ptr = getNextTag(xml, start, tag);
	while(ptr){
		if(tag == "DV_NODE_NAME") {
			_name = xml+start;
		}
		else if(tag == "DVOS_OBJ") {
			_obj->fromXML(xml+start);
		}
		xml = ptr;
		ptr = getNextTag(xml, start, tag);
	}

}

// ---------------- DvObject ----------------

void DvObject::toXML(std::ofstream &ofp) 
{	
	ofp << "   <DVOS_OBJ>\n";
	
		ofp << "    <SEQ_LEN>";
            ofp << (int) seqLen;
		ofp << "</SEQ_LEN>\n";
	
		for(size_t i=0; i<dims.size(); i++){
			ofp << "    <DIM>";
                ofp << (int) dims[i];
			ofp << "</DIM>\n";
		}	
	
		ofp << "    <SEC_RES>";
            ofp << secRes;
            int oldRes = 3; // this value 3 should never get used
            if(is_time() && Tdata.size() > 0) {
                oldRes = Tdata[0].setSecResolution(secRes);
            }
            else if(is_event() && Edata.size() > 0) {
                oldRes = Edata[0].start().setSecResolution(secRes);
            }
		ofp << "</SEC_RES>\n";
	
		ofp << "    <DATA_TYPE>";
            if(is_dbl()) ofp << "DV_DBL";
            else if(is_int()) ofp << "DV_INT";
            else if(is_time()) ofp << "DV_TIME";
            else if(is_event()) ofp << "DV_EVENT";
            else ofp << "DV_STR";
		ofp << "</DATA_TYPE>\n";

		ofp << "    <DATA>";
			if(is_dbl()){
				for(size_t i=0; i<totalSize()-1; i++){
					ofp << Ddata[i];
					ofp << ',';
				}
				if(totalSize() > 0) ofp << Ddata[totalSize()-1];
			}
			else if(is_int()){
				for(size_t i=0; i<totalSize()-1; i++){
					ofp << Idata[i];
					ofp << ',';
				}
				if(totalSize() > 0) ofp << Idata[totalSize()-1];
			}
			else{
				for(size_t i=0; i<totalSize()-1; i++){
					ofp << asStr(i);
					ofp << ',';
				}

				if(totalSize() > 0) ofp << asStr(totalSize()-1);
			}
			
            if(is_time() && Tdata.size() > 0) {
                Tdata[0].setSecResolution(oldRes);
            }			
            else if(is_event() && Edata.size() > 0) {
                Edata[0].start().setSecResolution(oldRes);
            }	
				
		ofp << "</DATA>\n";

		xrefs.toXML(ofp, DvString("META_DATA"));
				
	ofp << "   </DVOS_OBJ>\n";

	return;

}

void DvObject::fromXML(char *xml) 
{
	size_t start;
	DvString tag;
	DvString type;
	dims.clear();
	
	object_id=++ID_counter; // unique object id
	
	char *ptr = getNextTag(xml, start, tag);
	while(ptr){

		if(tag == "SEQ_LEN") {
			seqLen = (size_t) DvString(xml+start).toInt();
		}
		else if(tag == "DIM") {
			dims.push_back((size_t) DvString(xml+start).toInt());
		}
		else if(tag == "SEC_RES") {
			secRes = DvString(xml+start).toInt();
		}
		else if(tag == "DATA_TYPE") {
			type =xml+start;
		}
		else if(tag == "DATA") {
			if(seqLen == 0 ) return; // must find length before data
			
			if(type == "DV_DBL"){
				size_t len=arraySize()*seqLen;
				Ddata.resize(len, 0.0);
				if(len == 0) Ddata.resize(1, 0.0);
				
				xml += start;
				char* end = xml;
				for(size_t i=0; i<len; i++){
					end = strchr(end, ',');
					// Increment to next value only if delimiter found
					// if we run out of data we fill valarray with last element
					if(end){
						*end = '\0';
						end++;
					}
					Ddata[i] = strtod(xml, NULL);
					xml = end;
				}
			}
			else if(type == "DV_INT"){
				Ddata.resize(0); // empty dummy valarray
				size_t len=arraySize()*seqLen;
				Idata.resize(len, 0);
				if(len == 0) Idata.resize(1, 0);
				
				xml += start;
				char* end = xml;
				for(size_t i=0; i<len; i++){
					end = strchr(end, ',');
					// Increment to next value only if delimiter found
					// if we run out of data we fill valarray with last element
					if(end){
						*end = '\0';
						end++;
					}
					Idata[i] = (int) strtol(xml, NULL, 10);
					xml = end;
				}
			}
			else if(type == "DV_STR"){
				Ddata.resize(0); // empty dummy valarray
				size_t len=arraySize()*seqLen;
				Sdata.resize(len, DvString(""));
				if(len == 0) Sdata.resize(1, DvString(""));

				xml += start;
				char* end = xml;
				for(size_t i=0; i<len; i++){
					end = strchr(end, ',');
					// Increment to next value only if delimiter found
					// if we run out of data we fill valarray with last element
					if(end){
						*end = '\0';
						end++;
					}
					Sdata[i] = xml;
					xml = end;
				}
			}
			else if(type == "DV_TIME"){
				Ddata.resize(0); // empty dummy valarray
				size_t len=arraySize()*seqLen;
				Tdata.resize(len, DvTime());
				if(len == 0) Tdata.resize(1, DvTime());
				
				xml += start;
				char* end = xml;
				for(size_t i=0; i<len; i++){
					end = strchr(end, ',');
					// Increment to next value only if delimiter found
					// if we run out of data we fill valarray with last element
					if(end){
						*end = '\0';
						end++;
					}
					Tdata[i] = DvString(xml).toTime();
					xml = end;
				}
			}
			else if(type == "DV_EVENT"){
				Ddata.resize(0); // empty dummy valarray
				size_t len=arraySize()*seqLen;
				Edata.resize(len, DvEvent());
				if(len == 0) Edata.resize(1, DvEvent());
				
				xml += start;
				char* end = xml;
				for(size_t i=0; i<len; i++){
					end = strchr(end, ',');
					// Increment to next value only if delimiter found
					// if we run out of data we fill valarray with last element
					if(end){
						*end = '\0';
						end++;
					}
					Edata[i] = DvString(xml).toEvent();
					xml = end;
				}
			}
			else return; // no idea
				
		}
		else if(tag == "META_DATA") {
			xrefs.fromXML(xml+start);
		}
		
		xml = ptr;
		ptr = getNextTag(xml, start, tag);
	}

	if(dims.size() == 0) dims.push_back(1); // old savesets may have no dims on scalars

}

