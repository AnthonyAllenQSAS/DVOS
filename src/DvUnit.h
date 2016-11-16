/***********************************************************************
*           See Accompanying licence file QSAS_licence                 *
***********************************************************************/


#ifndef _DV_UNIT_H_
#define _DV_UNIT_H_

#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <memory.h>

 std::string  DvUnitGetBaseSI (const char * input_str);

 std::string  DvUnitGetTidySI (const char * input_str);

 std::string  DvUnitGetProductSI (const char * SIc1, const char * SIc2);
 
 std::string  DvUnitGetPowerSI (const char * input_str, float p);

 std::string  DvUnitMakeSIconvStr (double conv, const char * SIstr);

 std::string  DvUnitUnitsStr(const char * SIc);

 bool	      DvUnitHasUnits (const char * input_str);

 bool         DvUnitAreUnitsSame(const char *SIc1, const char *SIc2);

 bool         DvUnitSameBaseSI(const char *SIc1, const char *SIc2);

 double       DvUnitConvFactor(const char * SIc);

 double       DvUnitConversion(const char * SIcSource,  const char * SIcTarget);

/** Some useful constants  **/

#define QUNIT_TRUE            1
#define QUNIT_FALSE           0

#define QUNIT_CATALOGUE_ENTRIES  28   /* size of unit description table */


/** Global var for level of error messages **/

enum DvUnitVerbosity {
	QUNIT_VERBOSITY_DEBUG,      /* chat on stdout, plus below */
	QUNIT_VERBOSITY_WARNINGS,   /* warnings to stderr, plus below */
	QUNIT_VERBOSITY_QUIET       /* errors to stderr */
};

extern enum DvUnitVerbosity QUNIT_VERBOSITY;




/**  Unit description structure for catalogue of recognized units  **/

typedef struct DvUnitUnitDescription {

	const char* symbol;
	const char*   name;
                        /* indices to express unit in terms of */
        int    m;       /* the seven official SI base units... */
        int   kg;
        int    s;
        int    A;
        int    K;
        int  mol;
        int   cd;
        int  rad;       /* ... and the two SI supplementary units */
        int   sr;
        int  deg;       /* this is a non-standard extension we want */

} qunit_unit_description_t;



/** Symbol table element to hold unit expressions and constants **/


typedef struct DvUnitUnitEntry {

    struct DvUnitUnitEntry       *next;      /*  points to next entry        */
    struct DvUnitUnitDescription *SI_unit;   /*  points to catalogue entry   */
    float  exponent;                        /*  allow non-integer exponents */
    char*  text;                            /*  for non-catalogued items    */

} qunit_unit_entry_t;



/** Structure to describe a parsed SI_conversion attribute **/

typedef struct DvUnitSIConversion {

        float  value;
        struct DvUnitUnitEntry *first_entry;

        float    m;       /* cumulative base-unit-eqvt exponents */
        float   kg;
        float    s;
        float    A;
        float    K;
        float  mol;
        float   cd;
        float  rad;
        float   sr;
        float  deg;   

}  qunit_si_conversion;


/**  qunit catalogue routines  **/

void DvUnitInitialiseCatalogue ( );
void DvUnitPrintCatalogueDescription (const struct DvUnitUnitDescription*);
struct DvUnitUnitDescription* DvUnitLookupCatalogueDescription (const char*);

/**  qunit structure creation, destruction and output routines  **/

struct DvUnitSIConversion* DvUnitCreateSIConvDescrip ( );
void DvUnitDestroySIConvDescrip (struct DvUnitSIConversion*);

struct DvUnitUnitEntry*    DvUnitCreateUnitEntry ( );
void  DvUnitDestroyUnitEntry (struct DvUnitUnitEntry*);

void DvUnitStrCatSIConversionBase (char*, struct DvUnitSIConversion*);
void DvUnitStrCatSIConversionTidy (char*, struct DvUnitSIConversion*);
void DvUnitStrCatUnitEntry    (char*, struct DvUnitUnitEntry*);
void DvUnitStrCatPowerTidy (char* str_buf, struct DvUnitSIConversion *SIc, double p);


/**  qunit structure accumulating routines  **/

struct DvUnitUnitEntry* DvUnitAccumulateUnitEntry
                          (struct DvUnitUnitEntry*, struct DvUnitUnitEntry*);

struct DvUnitSIConversion* DvUnitAccumulateSIConversion
                                               (struct DvUnitSIConversion*);


/**  parsing routine  **/

int  DvUnitSIConvParse (struct DvUnitSIConversion*);

char * qunit_findNext();

void DvUnitAddUnit(struct DvUnitSIConversion * SI_c, 
                  char * unit, 
		  float power);
		  
double Qar_pow(double x, double y);

double DvUnit_atof(const char * vtxt);

#endif
