//
//  Dvar.h
//
//   Var pointer class
//

#ifndef Q_DVAR_H
#define Q_DVAR_H

#include <iostream>

class DvRecord;
class DvEvent;
class DvString;

//-----------------------------------------------------------------------------
// Class Dvar<T>
//-----------------------------------------------------------------------------
template<class T> class Dvar
{
	
  protected:
    T *p; // this should never be null or undefined
	
  public:
    // constructors (underlying object is always created)
    Dvar() {  p = new T(); p->inc_ref_count(); } // empty constructor

    Dvar( const Dvar<T>& var ) : p(var.p) // from Dvar 
    { p->inc_ref_count(); } 
	
    Dvar( T * _p ) : p(_p) { p->inc_ref_count(); } // from object ptr
	
    Dvar( T & _p ) : p(&_p) { p->inc_ref_count(); } // from object
	
    ~Dvar() { if( p->dec_ref_count() == 0 ) delete p; }
	
 	
    Dvar<T>& operator= (  T *x ) 
    { 
      x->inc_ref_count(); // increment first in case x.p is p already (var = var->fn())
	  if( p->dec_ref_count() == 0 ) {delete p;}
	  p = x; 
	  return *this;
    }
	
    Dvar<T>& operator= (const Dvar <T> &x ) 
    { 
      x.p->inc_ref_count(); // increment first in case x.p is p already (var = var->fn())
      if(  p->dec_ref_count() == 0 ) {delete p;}
	  p = x.p; 
	  return *this;
    }


    Dvar<T>& operator+= (const Dvar <T> &x ){
		*p += *(x.p);
        
		return *this;
	} 
    Dvar<T>& operator+= (const DvEvent &x ){
		*p += x;
		return *this;
	} 
    Dvar<T>& operator+= (const DvString &x ){
		*p += x;
		return *this;
	} 
    Dvar<T>& operator+= (const DvTime &x ){
		*p += x;
		return *this;
	} 
    Dvar<T>& operator+= (const double &x ){
		*p += x;
		return *this;
	} 
    Dvar<T>& operator+= (const int &x ){
		*p += x;
		return *this;
	} 


    Dvar<T>& operator*= (const Dvar <T> &x ){
		*p *= *(x.p);
		return *this;
	} 
    Dvar<T>& operator*= (const DvEvent &x ){
		*p *= x;
		return *this;
	} 
    Dvar<T>& operator*= (const DvString &x ){
		*p *= x;
		return *this;
	} 
    Dvar<T>& operator*= (const DvTime &x ){
		*p *= x;
		return *this;
	} 
    Dvar<T>& operator*= (const double &x ){
		*p *= x;
		return *this;
	} 
    Dvar<T>& operator*= (const int &x ){
		*p *= x;
		return *this;
	} 


    Dvar<T>& operator/= (const Dvar <T> &x ){
		*p /= *(x.p);
		return *this;
	}  
    Dvar<T>& operator/= (const DvEvent &x ){
		*p /= x;
		return *this;
	} 
    Dvar<T>& operator/= (const DvString &x ){
		*p /= x;
		return *this;
	} 
    Dvar<T>& operator/= (const DvTime &x ){
		*p /= x;
		return *this;
	} 
    Dvar<T>& operator/= (const double &x ){
		*p /= x;
		return *this;
	} 
    Dvar<T>& operator/= (const int &x ){
		*p /= x;
		return *this;
	}
    
	
    Dvar<T>& operator-= (const Dvar <T> &x ){
		*p -= *(x.p);
		return *this;
	}   
    Dvar<T>& operator-= (const DvEvent &x ){
		*p -= x;
		return *this;
	} 
    Dvar<T>& operator-= (const DvString &x ){
		*p -= x;
		return *this;
	} 
    Dvar<T>& operator-= (const DvTime &x ){
		*p -= x;
		return *this;
	} 
    Dvar<T>& operator-= (const double &x ){
		*p -= x;
		return *this;
	} 
    Dvar<T>& operator-= (const int &x ){
		*p -= x;
		return *this;
	}
    
	
    T* operator+ (const Dvar <T> &x ){
		return *p + *(x.p);
	} 
	
    T* operator- (const Dvar <T> &x ){
		return *p - *(x.p);
	} 

    T* operator* (const Dvar <T> &x ){
		return *p * *(x.p);
	} 
	
    T* operator/ (const Dvar <T> &x ){
		return *p / *(x.p);
	} 
	
    bool is_nil() const { return p->is_nil(); }
    bool is_ok() const { return p->is_ok(); }
	
	// Equality tests on Dvar pointers based  on same id,
	// for data content equality use e.g.  *obj1_var == *obj2_var
    bool operator==(const Dvar<T>& v) const{  return (p->get_id() == v.p->get_id());}
    bool operator!=(const Dvar<T>& v) const{ return (p->get_id() != v.p->get_id());}
	

	T* operator->() const { return p; };
	DvRecord operator[](size_t rec);
	
	//------
	
    T& operator*() { return *p;}
    const T& operator*() const { return *p;}
	
    T* ptr() { return p; }
    const T* ptr() const { return p; }

};


#endif // #ifndef Q_DVAR_H

