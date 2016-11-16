/***********************************************************************
*           See Accompanying licence file QSAS_licence                 *
***********************************************************************/

#include <iostream>
#include "DvUnit.h"

/* global variables */

static char qunit_line_buf[132];
static char *qunit_line_buf_ptr;

static int _qunit_catalogue_initialised = 0;

static qunit_unit_description_t qunit_catalogue[QUNIT_CATALOGUE_ENTRIES];


/*************************************************************************/

int DvUnitSIConvParse(struct DvUnitSIConversion* SI_conv)
{
  char * ptr = NULL;
  char * unit = NULL;
  float power;
  char * token;
  
  ptr = strchr (qunit_line_buf, '>');
  
  if (ptr == NULL) {
    printf("no separator found\n"); fflush(stdout);
    return -1;
  }
  
  ptr[0] = '\0';
  SI_conv->value = atof(qunit_line_buf_ptr);
  if(SI_conv->value == 0) {
    printf("non-numeric or zero conversion factor: %s\n", ptr);  fflush(stdout);
    return -1;
  }
  
  qunit_line_buf_ptr = &(ptr[1]); /* clear '>' */
  
  while( (token = qunit_findNext()) != NULL){ 
    
    /* remove space to keep token and buffer aligned */
    while(qunit_line_buf_ptr[0] == ' ' ) qunit_line_buf_ptr++;
 
    /* test for '(' */
  
    if(token[0] == '(') {
       
	   // Ignore bracketed units but be sure to remove power
	   
      /* find power if present */
      power = 1.0;
      qunit_line_buf_ptr += strlen(token)+1; // allow for ')'
	        
      while(qunit_line_buf_ptr[0] == ' ') qunit_line_buf_ptr++;
      
     
      if(qunit_line_buf_ptr[0] == '^') {
        qunit_line_buf_ptr++;
      
	    while (isdigit(qunit_line_buf_ptr[0]) 
	    || qunit_line_buf_ptr[0] == 'e' 
	    || qunit_line_buf_ptr[0] == '.' 
        || qunit_line_buf_ptr[0] == '+' 
	    || qunit_line_buf_ptr[0] == '-') qunit_line_buf_ptr++;

      }
	  
    } /* end on (text) */
    else{ /* SI unit */
      
      unit = &(token[0]);
    
      /* find power if present */
      power = 1.0;
      qunit_line_buf_ptr += strlen(token); 
      
      while(qunit_line_buf_ptr[0] == ' ') qunit_line_buf_ptr++;
      
     
      if(qunit_line_buf_ptr[0] == '^') {
        qunit_line_buf_ptr++;
        power = atof(qunit_line_buf_ptr);
        if (power == 0 ) {
          printf("power of unit zero or not found\n"); fflush(stdout);
          return -1;
        }
	while (isdigit(qunit_line_buf_ptr[0]) 
	|| qunit_line_buf_ptr[0] == '.' 
	|| qunit_line_buf_ptr[0] == '+' 
	|| qunit_line_buf_ptr[0] == '-') qunit_line_buf_ptr++;
      }
  
	  /* pass unit and power to cumulative structure */
    
      DvUnitAddUnit(SI_conv, unit, power);
    
    } /* end else on unit */
  

    
    
  } /* end while */
  
  SI_conv = DvUnitAccumulateSIConversion(SI_conv); // unpack units list
  
  return 0;

}


/*************************************************************************/
char * qunit_findNext(){
  static char localcopy[132];
  char * lptr;
  
  strncpy(&(localcopy[0]), qunit_line_buf_ptr, 131);
  
  
  lptr = strtok(&(localcopy[0]), " ^)");
  
  return lptr; 
}

/*************************************************************************/

void DvUnitAddUnit(struct DvUnitSIConversion * SI_c, char * unit, float power){
 
struct DvUnitUnitEntry *unit_entry_ptr;

unit_entry_ptr = DvUnitCreateUnitEntry( );
unit_entry_ptr->exponent = power;

if( unit[0] == '(' ) {
  unit_entry_ptr->text = (char *) malloc(strlen(unit));
  strcpy(unit_entry_ptr->text, &(unit[1]));
}
else{
  unit_entry_ptr->SI_unit = DvUnitLookupCatalogueDescription(unit);
}

 SI_c->first_entry = DvUnitAccumulateUnitEntry (SI_c->first_entry, unit_entry_ptr);
 return;
}


/*************************************************************************/
/*************************************************************************/

struct DvUnitSIConversion* DvUnitCreateSIConvDescrip ( )

{
  struct DvUnitSIConversion *SI_conversion_ptr;


  SI_conversion_ptr = (struct DvUnitSIConversion *) malloc(sizeof(struct DvUnitSIConversion));

  SI_conversion_ptr->value = 0.0;
  SI_conversion_ptr->first_entry = NULL;

  SI_conversion_ptr->m    =  0.0;
  SI_conversion_ptr->kg   =  0.0;
  SI_conversion_ptr->s    =  0.0;
  SI_conversion_ptr->A    =  0.0;
  SI_conversion_ptr->K    =  0.0;
  SI_conversion_ptr->mol  =  0.0;
  SI_conversion_ptr->cd   =  0.0;
  SI_conversion_ptr->rad  =  0.0;
  SI_conversion_ptr->sr   =  0.0;
  SI_conversion_ptr->deg  =  0.0;

  return SI_conversion_ptr;

}

/*************************************************************************/

void DvUnitDestroySIConvDescrip  (struct DvUnitSIConversion *SI_conversion_ptr)

{

  if (SI_conversion_ptr != NULL)  {
    if (SI_conversion_ptr->first_entry != NULL) free(SI_conversion_ptr->first_entry);

    free (SI_conversion_ptr);

  }

  return;

}

/*************************************************************************/

struct DvUnitUnitEntry* DvUnitCreateUnitEntry ( )

{
  struct DvUnitUnitEntry *unit_entry_ptr;


  unit_entry_ptr = (struct DvUnitUnitEntry *) malloc(sizeof(struct DvUnitUnitEntry));

  unit_entry_ptr->next = NULL;
  unit_entry_ptr->SI_unit = NULL;
  unit_entry_ptr->exponent = 0.0;
  unit_entry_ptr->text = NULL;

  return unit_entry_ptr;

}

/*************************************************************************/

void DvUnitDestroyUnitEntry (struct DvUnitUnitEntry *unit_entry_ptr)

{

  if (unit_entry_ptr != NULL) {

      if (unit_entry_ptr->next != NULL)
          DvUnitDestroyUnitEntry (unit_entry_ptr->next);

      if (unit_entry_ptr->text !=NULL)
          free(unit_entry_ptr->text);

      free (unit_entry_ptr);

  }

  return;

}

/*************************************************************************/

struct DvUnitUnitEntry*
       DvUnitAccumulateUnitEntry (struct DvUnitUnitEntry *entry_list,
                                         struct DvUnitUnitEntry *new_entry)

{
  struct DvUnitUnitEntry** list_ptr;

  /* if this is the first call start new list */
  
  if( entry_list == NULL) {
    entry_list = new_entry;
    return entry_list;
  }
  
  /** If new_entry is an SI_unit, look for an existing entry for this
      unit and accumulate the exponent if found, then free the new_entry.
      If new_entry is an SI_unit with no existing entry, or if it is a
      text entry, append it to the existing list (and don't free it!)
  **/


  if (new_entry->SI_unit != NULL)  {    /* seek existing entry for unit */
      list_ptr = &(entry_list);

      while (    (*list_ptr != NULL)
              && ((*list_ptr)->SI_unit != new_entry->SI_unit) )
          list_ptr = &(*list_ptr)->next;

      if (*list_ptr != NULL)  {     /* found existing unit entry */
          (*list_ptr)->exponent += new_entry->exponent;
          free(new_entry);
      }
      else
            *list_ptr = new_entry;  /* not found, so append it */

  }


  else if (new_entry->text != NULL)  {    /* text entry, append it */

      list_ptr = &(entry_list->next);

      while (*list_ptr != NULL)
          list_ptr = &(*list_ptr)->next;

      *list_ptr = new_entry;

  }

  else
  {
      printf("WARNING: DvUnitAccumulateUnitEntry: fell-thru ?SI_unit/text!\a\n");
	  fflush(stdout);
  }


  return entry_list;

}

/*************************************************************************/

struct DvUnitSIConversion* DvUnitAccumulateSIConversion
                                           (struct DvUnitSIConversion *SIc)

{


 /**  Accumulate base unit counts in SIConversion from list of unit
      entries, and remove any unit entries with exponent zero from list.
 **/

  struct DvUnitUnitEntry *unit_entry = NULL;
  struct DvUnitUnitEntry *prev_unit_entry = NULL;

  float  exponent;



  if (SIc != NULL)   {

      /* first remove any zero-exponent units from head of list */

      while (   (SIc->first_entry != NULL) 
             && (SIc->first_entry->SI_unit != NULL)
             && (-0.0001 < SIc->first_entry->exponent) 
             && ( 0.0001 > SIc->first_entry->exponent) )  {

          prev_unit_entry = SIc->first_entry;
          SIc->first_entry = SIc->first_entry->next;

          if (prev_unit_entry != NULL)
              free(prev_unit_entry);

      }
 

      if (SIc->first_entry != NULL)   {

          /* now remove any zero-exponent units in remainder of list */

          prev_unit_entry = SIc->first_entry;
          unit_entry = SIc->first_entry;
 
          while (unit_entry != NULL)  {

              if (unit_entry->SI_unit != NULL)  {

                  exponent = unit_entry->exponent;

                  if ( (exponent < -0.0001) || (0.0001 < exponent) )  {

                      SIc->m    +=  exponent * unit_entry->SI_unit->m;
                      SIc->kg   +=  exponent * unit_entry->SI_unit->kg;
                      SIc->s    +=  exponent * unit_entry->SI_unit->s;
                      SIc->A    +=  exponent * unit_entry->SI_unit->A;
                      SIc->K    +=  exponent * unit_entry->SI_unit->K;
                      SIc->mol  +=  exponent * unit_entry->SI_unit->mol;
                      SIc->cd   +=  exponent * unit_entry->SI_unit->cd;
                      SIc->rad  +=  exponent * unit_entry->SI_unit->rad;
                      SIc->sr   +=  exponent * unit_entry->SI_unit->sr;
                      SIc->deg  +=  exponent * unit_entry->SI_unit->deg;

                      prev_unit_entry = unit_entry;
                      unit_entry = unit_entry->next;

                  }

                  else  {     /* remove from list if expnt = zero */
 
                      prev_unit_entry->next = unit_entry->next;
                      free(unit_entry);
                      unit_entry = prev_unit_entry->next;

                  }    /*  end of: zero-expnt test */

              }

              else   {

                  prev_unit_entry = unit_entry;
                  unit_entry = unit_entry->next;

              }   /*  end of: if (unit_entry->SI_unit != NULL)  */

          }  /*  end of: while (unit_entry != NULL)  */

      }  /*  end of: if (SIc->first_entry != NULL)  */


  }  /*  end of: if (SIc != NULL)...  */


  return SIc;

}

/*************************************************************************/

void DvUnitStrCatSIConversionBase (char* str_buf, struct DvUnitSIConversion *SIc)

{

  static char local_buf[100];

  int   ii;
  const char *symb[10] = {"m","kg","s","A","K","mol","cd", "rad","sr","deg"};
  float expnt[10];

  if (SIc != NULL)  {
  
    // convert deg to rad
      
    if(SIc->value != 0){
        
        double conv = atan(1.0) / 45.0;
        SIc->value *= pow(conv, SIc->deg);
        SIc->rad += SIc->deg; 
        SIc->deg = 0; 
    }
      
    // need conversion value as well as base units 
    sprintf(local_buf, "%.5g >", (double)SIc->value);  
    strcat(str_buf, local_buf);


    expnt[0]  = SIc->m;
    expnt[1]  = SIc->kg;
    expnt[2]  = SIc->s;
    expnt[3]  = SIc->A;
    expnt[4]  = SIc->K;
    expnt[5]  = SIc->mol;
    expnt[6]  = SIc->cd;
    expnt[7]  = SIc->rad;
    expnt[8]  = SIc->sr;
    expnt[9]  = SIc->deg;


    for (ii  =0; ii < 9; ii++)  { // deg not used
      if (!( (-0.0001 < expnt[ii]) && (expnt[ii] < 0.0001) ))  {    /* != 0 */

        if ( (0.9999 < expnt[ii]) && (expnt[ii] < 1.0001) )    /* expnt = 1 */
          sprintf(local_buf, " %s", symb[ii]);
        else                                           /* expnt !=0 and !=1 */
          sprintf(local_buf, " %s^%g", symb[ii], expnt[ii]);

        strcat(str_buf, local_buf);

      }

    }

  }
  return;

}

/*************************************************************************/


void DvUnitInitialiseCatalogue ( )

{

  /**  Inits catalogue of recognized SI units and their properties  **/


  if ( _qunit_catalogue_initialised == 1) return;

  _qunit_catalogue_initialised = 1;



  /* describe the SI base and supplementary units, and our "degree" */


  qunit_catalogue[0].symbol = "m";
  qunit_catalogue[0].name   = "metre";
  qunit_catalogue[0].m   = 1;
  qunit_catalogue[0].kg  = 0;
  qunit_catalogue[0].s   = 0;
  qunit_catalogue[0].A   = 0;
  qunit_catalogue[0].K   = 0;
  qunit_catalogue[0].mol = 0;
  qunit_catalogue[0].cd  = 0;
  qunit_catalogue[0].rad = 0;
  qunit_catalogue[0].sr  = 0;
  qunit_catalogue[0].deg = 0;

  qunit_catalogue[1].symbol = "kg";
  qunit_catalogue[1].name   = "kilogram";
  qunit_catalogue[1].m   = 0;
  qunit_catalogue[1].kg  = 1;
  qunit_catalogue[1].s   = 0;
  qunit_catalogue[1].A   = 0;
  qunit_catalogue[1].K   = 0;
  qunit_catalogue[1].mol = 0;
  qunit_catalogue[1].cd  = 0;
  qunit_catalogue[1].rad = 0;
  qunit_catalogue[1].sr  = 0;
  qunit_catalogue[1].deg = 0;

  qunit_catalogue[2].symbol = "s";
  qunit_catalogue[2].name   = "second";
  qunit_catalogue[2].m   = 0;
  qunit_catalogue[2].kg  = 0;
  qunit_catalogue[2].s   = 1;
  qunit_catalogue[2].A   = 0;
  qunit_catalogue[2].K   = 0;
  qunit_catalogue[2].mol = 0;
  qunit_catalogue[2].cd  = 0;
  qunit_catalogue[2].rad = 0;
  qunit_catalogue[2].sr  = 0;
  qunit_catalogue[2].deg = 0;

  qunit_catalogue[3].symbol = "A";
  qunit_catalogue[3].name   = "ampere";
  qunit_catalogue[3].m   = 0;
  qunit_catalogue[3].kg  = 0;
  qunit_catalogue[3].s   = 0;
  qunit_catalogue[3].A   = 1;
  qunit_catalogue[3].K   = 0;
  qunit_catalogue[3].mol = 0;
  qunit_catalogue[3].cd  = 0;
  qunit_catalogue[3].rad = 0;
  qunit_catalogue[3].sr  = 0;
  qunit_catalogue[3].deg = 0;

  qunit_catalogue[4].symbol = "K";
  qunit_catalogue[4].name   = "kelvin";
  qunit_catalogue[4].m   = 0;
  qunit_catalogue[4].kg  = 0;
  qunit_catalogue[4].s   = 0;
  qunit_catalogue[4].A   = 0;
  qunit_catalogue[4].K   = 1;
  qunit_catalogue[4].mol = 0;
  qunit_catalogue[4].cd  = 0;
  qunit_catalogue[4].rad = 0;
  qunit_catalogue[4].sr  = 0;
  qunit_catalogue[4].deg = 0;

  qunit_catalogue[5].symbol = "mol";
  qunit_catalogue[5].name   = "mole";
  qunit_catalogue[5].m   = 0;
  qunit_catalogue[5].kg  = 0;
  qunit_catalogue[5].s   = 0;
  qunit_catalogue[5].A   = 0;
  qunit_catalogue[5].K   = 0;
  qunit_catalogue[5].mol = 1;
  qunit_catalogue[5].cd  = 0;
  qunit_catalogue[5].rad = 0;
  qunit_catalogue[5].sr  = 0;
  qunit_catalogue[5].deg = 0;

  qunit_catalogue[6].symbol = "cd";
  qunit_catalogue[6].name   = "candela";
  qunit_catalogue[6].m   = 0;
  qunit_catalogue[6].kg  = 0;
  qunit_catalogue[6].s   = 0;
  qunit_catalogue[6].A   = 0;
  qunit_catalogue[6].K   = 0;
  qunit_catalogue[6].mol = 0;
  qunit_catalogue[6].cd  = 1;
  qunit_catalogue[6].rad = 0;
  qunit_catalogue[6].sr  = 0;
  qunit_catalogue[6].deg = 0;

  qunit_catalogue[7].symbol = "rad";
  qunit_catalogue[7].name   = "radian";
  qunit_catalogue[7].m   = 0;
  qunit_catalogue[7].kg  = 0;
  qunit_catalogue[7].s   = 0;
  qunit_catalogue[7].A   = 0;
  qunit_catalogue[7].K   = 0;
  qunit_catalogue[7].mol = 0;
  qunit_catalogue[7].cd  = 0;
  qunit_catalogue[7].rad = 1;
  qunit_catalogue[7].sr  = 0;
  qunit_catalogue[7].deg = 0;

  qunit_catalogue[8].symbol = "sr";
  qunit_catalogue[8].name   = "steradian";
  qunit_catalogue[8].m   = 0;
  qunit_catalogue[8].kg  = 0;
  qunit_catalogue[8].s   = 0;
  qunit_catalogue[8].A   = 0;
  qunit_catalogue[8].K   = 0;
  qunit_catalogue[8].mol = 0;
  qunit_catalogue[8].cd  = 0;
  qunit_catalogue[8].rad = 0;
  qunit_catalogue[8].sr  = 1;
  qunit_catalogue[8].deg = 0;

  qunit_catalogue[9].symbol = "deg";
  qunit_catalogue[9].name   = "degree";
  qunit_catalogue[9].m   = 0;
  qunit_catalogue[9].kg  = 0;
  qunit_catalogue[9].s   = 0;
  qunit_catalogue[9].A   = 0;
  qunit_catalogue[9].K   = 0;
  qunit_catalogue[9].mol = 0;
  qunit_catalogue[9].cd  = 0;
  qunit_catalogue[9].rad = 0;
  qunit_catalogue[9].sr  = 0;
  qunit_catalogue[9].deg = 1;



  /* describe derived units (append further entries as required) */


  qunit_catalogue[10].symbol = "Hz";
  qunit_catalogue[10].name   = "hertz";
  qunit_catalogue[10].m   =  0;
  qunit_catalogue[10].kg  =  0;
  qunit_catalogue[10].s   = -1;
  qunit_catalogue[10].A   =  0;
  qunit_catalogue[10].K   =  0;
  qunit_catalogue[10].mol =  0;
  qunit_catalogue[10].cd  =  0;
  qunit_catalogue[10].rad =  0;
  qunit_catalogue[10].sr  =  0;
  qunit_catalogue[10].deg =  0;

  qunit_catalogue[11].symbol = "N";
  qunit_catalogue[11].name   = "newton";
  qunit_catalogue[11].m   =  1;
  qunit_catalogue[11].kg  =  1;
  qunit_catalogue[11].s   = -2;
  qunit_catalogue[11].A   =  0;
  qunit_catalogue[11].K   =  0;
  qunit_catalogue[11].mol =  0;
  qunit_catalogue[11].cd  =  0;
  qunit_catalogue[11].rad =  0;
  qunit_catalogue[11].sr  =  0;
  qunit_catalogue[11].deg =  0;

  qunit_catalogue[12].symbol = "Pa";
  qunit_catalogue[12].name   = "pascal";
  qunit_catalogue[12].m   = -1;
  qunit_catalogue[12].kg  =  1;
  qunit_catalogue[12].s   = -2;
  qunit_catalogue[12].A   =  0;
  qunit_catalogue[12].K   =  0;
  qunit_catalogue[12].mol =  0;
  qunit_catalogue[12].cd  =  0;
  qunit_catalogue[12].rad =  0;
  qunit_catalogue[12].sr  =  0;
  qunit_catalogue[12].deg =  0;

  qunit_catalogue[13].symbol = "J";
  qunit_catalogue[13].name   = "joule";
  qunit_catalogue[13].m   =  2;
  qunit_catalogue[13].kg  =  1;
  qunit_catalogue[13].s   = -2;
  qunit_catalogue[13].A   =  0;
  qunit_catalogue[13].K   =  0;
  qunit_catalogue[13].mol =  0;
  qunit_catalogue[13].cd  =  0;
  qunit_catalogue[13].rad =  0;
  qunit_catalogue[13].sr  =  0;
  qunit_catalogue[13].deg =  0;

  qunit_catalogue[14].symbol = "W";
  qunit_catalogue[14].name   = "watt";
  qunit_catalogue[14].m   =  2;
  qunit_catalogue[14].kg  =  1;
  qunit_catalogue[14].s   = -3;
  qunit_catalogue[14].A   =  0;
  qunit_catalogue[14].K   =  0;
  qunit_catalogue[14].mol =  0;
  qunit_catalogue[14].cd  =  0;
  qunit_catalogue[14].rad =  0;
  qunit_catalogue[14].sr  =  0;
  qunit_catalogue[14].deg =  0;

  qunit_catalogue[15].symbol = "C";
  qunit_catalogue[15].name   = "coulomb";
  qunit_catalogue[15].m   =  0;
  qunit_catalogue[15].kg  =  0;
  qunit_catalogue[15].s   =  1;
  qunit_catalogue[15].A   =  1;
  qunit_catalogue[15].K   =  0;
  qunit_catalogue[15].mol =  0;
  qunit_catalogue[15].cd  =  0;
  qunit_catalogue[15].rad =  0;
  qunit_catalogue[15].sr  =  0;
  qunit_catalogue[15].deg =  0;

  qunit_catalogue[16].symbol = "V";
  qunit_catalogue[16].name   = "volt";
  qunit_catalogue[16].m   =  2;
  qunit_catalogue[16].kg  =  1;
  qunit_catalogue[16].s   = -3;
  qunit_catalogue[16].A   = -1;
  qunit_catalogue[16].K   =  0;
  qunit_catalogue[16].mol =  0;
  qunit_catalogue[16].cd  =  0;
  qunit_catalogue[16].rad =  0;
  qunit_catalogue[16].sr  =  0;
  qunit_catalogue[16].deg =  0;

  qunit_catalogue[17].symbol = "F";
  qunit_catalogue[17].name   = "farad";
  qunit_catalogue[17].m   = -2;
  qunit_catalogue[17].kg  = -1;
  qunit_catalogue[17].s   =  4;
  qunit_catalogue[17].A   =  2;
  qunit_catalogue[17].K   =  0;
  qunit_catalogue[17].mol =  0;
  qunit_catalogue[17].cd  =  0;
  qunit_catalogue[17].rad =  0;
  qunit_catalogue[17].sr  =  0;
  qunit_catalogue[17].deg =  0;

  qunit_catalogue[18].symbol = "ohm";
  qunit_catalogue[18].name   = "ohm";
  qunit_catalogue[18].m   =  2;
  qunit_catalogue[18].kg  =  1;
  qunit_catalogue[18].s   = -3;
  qunit_catalogue[18].A   = -2;
  qunit_catalogue[18].K   =  0;
  qunit_catalogue[18].mol =  0;
  qunit_catalogue[18].cd  =  0;
  qunit_catalogue[18].rad =  0;
  qunit_catalogue[18].sr  =  0;
  qunit_catalogue[18].deg =  0;

  qunit_catalogue[19].symbol = "S";
  qunit_catalogue[19].name   = "siemens";
  qunit_catalogue[19].m   = -2;
  qunit_catalogue[19].kg  = -1;
  qunit_catalogue[19].s   =  3;
  qunit_catalogue[19].A   =  2;
  qunit_catalogue[19].K   =  0;
  qunit_catalogue[19].mol =  0;
  qunit_catalogue[19].cd  =  0;
  qunit_catalogue[19].rad =  0;
  qunit_catalogue[19].sr  =  0;
  qunit_catalogue[19].deg =  0;

  qunit_catalogue[20].symbol = "Wb";
  qunit_catalogue[20].name   = "weber";
  qunit_catalogue[20].m   =  2;
  qunit_catalogue[20].kg  =  1;
  qunit_catalogue[20].s   = -2;
  qunit_catalogue[20].A   = -1;
  qunit_catalogue[20].K   =  0;
  qunit_catalogue[20].mol =  0;
  qunit_catalogue[20].cd  =  0;
  qunit_catalogue[20].rad =  0;
  qunit_catalogue[20].sr  =  0;
  qunit_catalogue[20].deg =  0;

  qunit_catalogue[21].symbol = "T";
  qunit_catalogue[21].name   = "tesla";
  qunit_catalogue[21].m   =  0;
  qunit_catalogue[21].kg  =  1;
  qunit_catalogue[21].s   = -2;
  qunit_catalogue[21].A   = -1;
  qunit_catalogue[21].K   =  0;
  qunit_catalogue[21].mol =  0;
  qunit_catalogue[21].cd  =  0;
  qunit_catalogue[21].rad =  0;
  qunit_catalogue[21].sr  =  0;
  qunit_catalogue[21].deg =  0;

  qunit_catalogue[22].symbol = "H";
  qunit_catalogue[22].name   = "henry";
  qunit_catalogue[22].m   =  2;
  qunit_catalogue[22].kg  =  1;
  qunit_catalogue[22].s   = -2;
  qunit_catalogue[22].A   = -2;
  qunit_catalogue[22].K   =  0;
  qunit_catalogue[22].mol =  0;
  qunit_catalogue[22].cd  =  0;
  qunit_catalogue[22].rad =  0;
  qunit_catalogue[22].sr  =  0;
  qunit_catalogue[22].deg =  0;

  qunit_catalogue[23].symbol = "lm";
  qunit_catalogue[23].name   = "lumen";
  qunit_catalogue[23].m   =  0;
  qunit_catalogue[23].kg  =  0;
  qunit_catalogue[23].s   =  0;
  qunit_catalogue[23].A   =  0;
  qunit_catalogue[23].K   =  0;
  qunit_catalogue[23].mol =  0;
  qunit_catalogue[23].cd  =  1;
  qunit_catalogue[23].rad =  0;
  qunit_catalogue[23].sr  =  1;
  qunit_catalogue[23].deg =  0;

  qunit_catalogue[24].symbol = "lx";
  qunit_catalogue[24].name   = "lux";
  qunit_catalogue[24].m   = -2;
  qunit_catalogue[24].kg  =  0;
  qunit_catalogue[24].s   =  0;
  qunit_catalogue[24].A   =  0;
  qunit_catalogue[24].K   =  0;
  qunit_catalogue[24].mol =  0;
  qunit_catalogue[24].cd  =  1;
  qunit_catalogue[24].rad =  0;
  qunit_catalogue[24].sr  =  1;
  qunit_catalogue[24].deg =  0;

  qunit_catalogue[25].symbol = "Bq";
  qunit_catalogue[25].name   = "becquerel";
  qunit_catalogue[25].m   =  0;
  qunit_catalogue[25].kg  =  0;
  qunit_catalogue[25].s   = -1;
  qunit_catalogue[25].A   =  0;
  qunit_catalogue[25].K   =  0;
  qunit_catalogue[25].mol =  0;
  qunit_catalogue[25].cd  =  0;
  qunit_catalogue[25].rad =  0;
  qunit_catalogue[25].sr  =  0;
  qunit_catalogue[25].deg =  0;

  qunit_catalogue[26].symbol = "Gy";
  qunit_catalogue[26].name   = "gray";
  qunit_catalogue[26].m   =  2;
  qunit_catalogue[26].kg  =  0;
  qunit_catalogue[26].s   = -2;
  qunit_catalogue[26].A   =  0;
  qunit_catalogue[26].K   =  0;
  qunit_catalogue[26].mol =  0;
  qunit_catalogue[26].cd  =  0;
  qunit_catalogue[26].rad =  0;
  qunit_catalogue[26].sr  =  0;
  qunit_catalogue[26].deg =  0;

  qunit_catalogue[27].symbol = "unitless";
  qunit_catalogue[27].name   = "unitless";
  qunit_catalogue[27].m   =  0;
  qunit_catalogue[27].kg  =  0;
  qunit_catalogue[27].s   =  0;
  qunit_catalogue[27].A   =  0;
  qunit_catalogue[27].K   =  0;
  qunit_catalogue[27].mol =  0;
  qunit_catalogue[27].cd  =  0;
  qunit_catalogue[27].rad =  0;
  qunit_catalogue[27].sr  =  0;
  qunit_catalogue[27].deg =  0;
 

  return;

}


/*************************************************************************/

struct DvUnitUnitDescription* DvUnitLookupCatalogueDescription
                                                 (const char* test_string)

{

  /**  Lookup a unit in the catalogue by symbol or name and        **/
  /**  return a pointer to its description, or NULL if not found.  **/


  struct DvUnitUnitDescription *matching_description;
  int icount;

  char* lc_test_string;   /* lowercase copy of test_string if reqd */



  /** First try for an exact case-sensitive match of unit symbol **/

  icount = 0;
  matching_description = NULL;

  while ( (icount < QUNIT_CATALOGUE_ENTRIES)
                                   && (NULL == matching_description) )  {

      if ( 0 == strcmp(qunit_catalogue[icount].symbol, test_string) )
          matching_description = &qunit_catalogue[icount];

      ++icount;
  }


  /**
    If symbol was not matched, try for case-insensitive name match.
    Names are stored lowercase as per SI standard, so can do this
    by testing for an exact match of a lowercase copy of test_string.
  **/

  if (NULL == matching_description)  {

      lc_test_string = (char*) malloc(strlen(test_string)+1);
      icount = 0;
      while ((lc_test_string[icount] = tolower(test_string[icount])) != '\0')
          ++icount;

      icount = 0;

      while ( (icount < QUNIT_CATALOGUE_ENTRIES)
                                       && (NULL == matching_description) )  {

          if ( 0 == strcmp(qunit_catalogue[icount].name, lc_test_string) )
              matching_description = &qunit_catalogue[icount];

          ++icount;
      }

      free(lc_test_string);

  }


  return matching_description;   /* NULL if no match was made */

}

/*************************************************************************/

void DvUnitPrintCatalogueDescription
                             (const struct DvUnitUnitDescription* descrip)

{
 
  /**  Print qunit catalogue description of an SI unit  **/


  printf("DvUnit catalogue: %s (%s) =", descrip->name, descrip->symbol);

    if (descrip->m   != 0)   printf(" m^%-2d",   descrip->m);
    if (descrip->kg  != 0)   printf(" kg^%-2d",  descrip->kg);
    if (descrip->s   != 0)   printf(" s^%-2d",   descrip->s);
    if (descrip->A   != 0)   printf(" A^%-2d",   descrip->A);
    if (descrip->K   != 0)   printf(" K^%-2d",   descrip->K);
    if (descrip->mol != 0)   printf(" mol^%-2d", descrip->mol);
    if (descrip->cd  != 0)   printf(" cd^%-2d",  descrip->cd);
    if (descrip->rad != 0)   printf(" rad^%-2d", descrip->rad);
    if (descrip->sr  != 0)   printf(" sr^%-2d",  descrip->sr);
    if (descrip->deg != 0)   printf(" deg^%-2d", descrip->deg);

  printf("\n");
  fflush(stdout);

  return;

}

/****************************************************************************/
/****************************************************************************/


std::string  DvUnitGetBaseSI (const char * input_str){

  struct DvUnitSIConversion *SI_conv = DvUnitCreateSIConvDescrip();
  
  char ret_buf[300];
 
  DvUnitInitialiseCatalogue ( );
  
  qunit_line_buf[0] = '\0';
  strncat(qunit_line_buf, input_str, 131);
  qunit_line_buf_ptr = qunit_line_buf;

  if ( 0 == DvUnitSIConvParse (SI_conv) )  {    /* CALL PARSER */
 
	ret_buf[0] = '\0';
	DvUnitStrCatSIConversionBase (ret_buf, SI_conv);

  }
  else
  {
	printf("DvUnitGetBaseSI: Error parsing /%s/, giving up.\n", input_str);
	fflush(stdout);
  }

  DvUnitDestroySIConvDescrip(SI_conv);

  return std::string(&(ret_buf[0]));

}

/****************************************************************************/


std::string DvUnitGetTidySI (const char * input_str)
{
  struct DvUnitSIConversion *SI_conv;

  SI_conv = DvUnitCreateSIConvDescrip();

  char ret_buf[300];

  DvUnitInitialiseCatalogue ( );

  qunit_line_buf[0] = '\0';
  strncat(qunit_line_buf, input_str, 131);
  qunit_line_buf_ptr = qunit_line_buf;

  if ( 0 == DvUnitSIConvParse (SI_conv) )  {    /* CALL PARSER */

     ret_buf[0] = '\0';
     DvUnitStrCatSIConversionTidy (&(ret_buf[0]), SI_conv);

  }

  else
  {
     printf("DvUnitGetTidySI: Error parsing /%s/, giving up.\n", qunit_line_buf);
     fflush(stdout);
  }

  DvUnitDestroySIConvDescrip(SI_conv);
     
  return std::string(&(ret_buf[0]));

} /* end DvUnitGetTidySI */

/*************************************************************************/


void DvUnitStrCatUnitEntry ( char* str_buf, struct DvUnitUnitEntry *ue )

{
	static char local_buf[100];

	if (ue != NULL)  {

		if (ue->SI_unit != NULL)  {

			if( (ue->exponent < 0.0001) && (-.0001 < ue->exponent) ){
				/* skip cancelled units */
				local_buf[0] = '\0';
			}
			else if(strcmp(ue->SI_unit->symbol, "unitless") == 0 ){
				/* skip unitless */
				local_buf[0] = '\0';
			}
			else if ( (ue->exponent < 0.9999) || (1.0001 < ue->exponent) )
				sprintf(local_buf, "%s^%g", ue->SI_unit->symbol, ue->exponent);
			else
				sprintf(local_buf, "%s", ue->SI_unit->symbol);

			strcat(str_buf, local_buf);
		}

		else if (ue->text != NULL){
			strcat(str_buf, "(");
			strcat(str_buf, ue->text);
			strcat(str_buf, ")");
		}

		else
			strcat (str_buf, "WARNING: DvUnitStrCatUnitEntry called with empty arg!");

	}

	else
		strcat(str_buf, "WARNING: DvUnitStrCatUnitEntry called with NULL arg!");


  return;

}

/*************************************************************************/

void DvUnitStrCatSIConversionTidy (char* str_buf, struct DvUnitSIConversion *SIc)

{

  struct DvUnitUnitEntry *next_unit_entry;
  static char local_buf[100];


  if (SIc != NULL)  {


    next_unit_entry = SIc->first_entry;

    sprintf(local_buf, "%.5g >", SIc->value);
    strcat(str_buf, local_buf);

    while (next_unit_entry != NULL)  {

      strcat(str_buf," ");
      DvUnitStrCatUnitEntry(str_buf, next_unit_entry);
      next_unit_entry = next_unit_entry->next;

    }

    /* if it all cancelled out show as unitless */
    char *ptr0 = strrchr(str_buf, '>');
    char *ptr = &(ptr0[1]);
    while (ptr[0] == ' ') ptr ++;
    if(ptr[0] == '\0') {
      ptr0[1] = '\0';
      strcat(str_buf," (unitless)");
    }

  }
  return;

}
/****************************************************************************/

bool   DvUnitHasUnits (const char * input_str){

  struct DvUnitSIConversion *SI_conv;
  
  SI_conv = DvUnitCreateSIConvDescrip();
  
  const char * ptr;
  float expnt = 0;
  
  /* trap empty input SI_conversion string */
  
  if(input_str == NULL) return false;
  ptr = &(input_str[0]);
  while( (ptr[0] == ' ') || (ptr[0] == '\n')) ptr ++;
  if(ptr[0] == '\0') return false;
  
  /* string is not empty, so look for units */
  
  DvUnitInitialiseCatalogue ( );

    qunit_line_buf[0] = '\0';
    strncat(qunit_line_buf, input_str, 131);
    qunit_line_buf_ptr = qunit_line_buf;

    if ( 0 == DvUnitSIConvParse (SI_conv) )  {    /* CALL PARSER */

      /* resolve units to base units */
      
      SI_conv = DvUnitAccumulateSIConversion(SI_conv);
      
      expnt  += fabs(SI_conv->m);
      expnt  += fabs(SI_conv->kg);
      expnt  += fabs(SI_conv->s);
      expnt  += fabs(SI_conv->A);
      expnt  += fabs(SI_conv->K);
      expnt  += fabs(SI_conv->mol);
      expnt  += fabs(SI_conv->cd);
      expnt  += fabs(SI_conv->rad);
      expnt  += fabs(SI_conv->sr);
      expnt  += fabs(SI_conv->deg);

      
        if ( expnt > 0.0001 )  {    /* != 0 */
          /* if any unit found return yes */
          DvUnitDestroySIConvDescrip(SI_conv);
          return true;
        }

    }

    else
	{
		printf("DvUnitHasUnits:Error parsing /%s/, giving up.\n", qunit_line_buf);
		fflush(stdout);
	}

    DvUnitDestroySIConvDescrip(SI_conv);  

  return false;

}

/****************************************************************************/

std::string  DvUnitGetProductSI (const char * SIc1, const char * SIc2)
{
  
  char local_buf[300];
  if(SIc1 == NULL && SIc2 == NULL)  return std::string("");
  if(SIc1 == NULL || strlen(SIc1) == 0)  return SIc2;
  if(SIc2 == NULL || strlen(SIc2) == 0)  return SIc1;

  double Conv1 = DvUnitConvFactor(SIc1);
  double Conv2 = DvUnitConvFactor(SIc2);
  double newConv = Conv1*Conv2;
  sprintf(local_buf, "%g >", newConv);
  std::string working = std::string(local_buf);
  
  working += DvUnitUnitsStr(SIc1) + " ";
  working += DvUnitUnitsStr(SIc2);

  return DvUnitGetTidySI(working.c_str()) ;

}

/****************************************************************************/

std::string  DvUnitMakeSIconvStr (double conv, const char * SIstr)
{
  
  char local_buf[300];
  if(SIstr == NULL || strlen(SIstr) == 0)  return std::string("");

  sprintf(local_buf, "%g >", conv);
  std::string working = std::string(local_buf);
  
  working += std::string(SIstr);

  return DvUnitGetTidySI(working.c_str()) ;

}

/****************************************************************************/

std::string  DvUnitGetPowerSI (const char * input_str, float p)
{

  struct DvUnitSIConversion* SI_conv;
  SI_conv = DvUnitCreateSIConvDescrip();
  
  char ret_buf[300];
  ret_buf[0] = '\0';
  if(input_str == NULL || strlen(input_str) == 0)  return std::string("");

  DvUnitInitialiseCatalogue ( );

  qunit_line_buf[0] = '\0';
  strncat(qunit_line_buf, input_str, 131);
  qunit_line_buf_ptr = qunit_line_buf;


  if ( 0 == DvUnitSIConvParse (SI_conv) )  {    /* CALL PARSER */
   
     DvUnitStrCatPowerTidy (&(ret_buf[0]), SI_conv, (double) p);

  }

  else
  {
     printf("DvUnitGetPowerSI: Error parsing %s, giving up.\n", input_str);
     fflush(stdout);
  }

  DvUnitDestroySIConvDescrip(SI_conv);

  return std::string(&(ret_buf[0]));

} /* end DvUnitGetPowerSI */


/*************************************************************************/


void DvUnitStrCatPowerTidy (char* str_buf, 
                           struct DvUnitSIConversion *SIc, 
                           double p)

{

  struct DvUnitUnitEntry *next_unit_entry;
  static char local_buf[100];
  double base_value;
  double new_value;

  base_value = (double) SIc->value;

  if (SIc != NULL)  {
    if( SIc->value <= 0) return;
    new_value = Qar_pow(base_value, p);
    
    next_unit_entry = SIc->first_entry;

    sprintf(local_buf, "%g >", new_value);
    strcat(str_buf, local_buf);

    while (next_unit_entry != NULL)  {
      next_unit_entry->exponent *= p;
      strcat(str_buf," ");
      DvUnitStrCatUnitEntry(str_buf, next_unit_entry);
      next_unit_entry = next_unit_entry->next;

    }
    
    /* if it all cancelled out show as ratio */
    char *ptr0 = strrchr(str_buf, '>');
    char *ptr = &(ptr0[1]);
    while (ptr[0] == ' ') ptr ++;
    if(ptr == '\0') {
      ptr0[1] = '\0';
      strcat(str_buf," (ratio)");
    }
       
  }
  return;

}
double Qar_pow(double x, double y){

 double value;

   if (x == 0) value = (double) 0.0;
   if (x <  0) value = pow((double)x, (int)y);
   else value = (double) pow((double)x, (double)y);

 return value;
}

/****************************************************************************/

bool DvUnitAreUnitsSame(const char *SIc1, const char *SIc2){
  
  if(SIc1 == NULL || SIc2 == NULL) return false;
  
  // if parsing fails, the input string is returned from DvUnitGetBaseSI() so safe to continue
  std::string base_str1 = DvUnitGetBaseSI(SIc1);  
  std::string base_str2 = DvUnitGetBaseSI(SIc2);

  if(base_str1 == base_str2) return true;

  return false;
  
} 

bool DvUnitSameBaseSI(const char *SIc1, const char *SIc2){
 
  char * ptr1, * ptr2;
 
  if(SIc1 == NULL || SIc2 == NULL) return false;

  std::string base_str1 = DvUnitGetBaseSI(SIc1); 
  std::string base_str2 = DvUnitGetBaseSI(SIc2);
  
  ptr1 = (char*)strstr(base_str1.c_str(), ">");
  ptr2 = (char*)strstr(base_str2.c_str(), ">");

  if(ptr1 == NULL || ptr2 == NULL || strcmp(ptr1, ptr2) !=0 ) return false;
  
  return true;
  
}


double DvUnitConvFactor(const char * SIc){

  size_t ptr;
  std::string newSIc(SIc);
  
  double value=0.0;
  ptr = newSIc.find( ">");
  if (ptr != std::string::npos){
	
    value = (double) DvUnit_atof(newSIc.substr(0, ptr).c_str());
  }
  
  return value;
  
}

double DvUnit_atof(const char * vtxt){

	// corrects for use of tex syntax on LHS of SI_CONVERSION
  double value=0.0;
  char *ptr = (char*)strstr(vtxt, "10^");
  if (ptr != NULL){
    double power = atof(&(ptr[3]));
    ptr[0] = '\0'; 
    value = (double) atof(vtxt);
    if(value < 1e-13) value = 1.0; // if vtxt of form "10^6" value is zero since there is no mantissa
    value *= pow(10., power);
  }
  else
      value = (double) atof(vtxt);

  return value;
  
}


double DvUnitConversion(const char * SIcSource,  const char * SIcTarget){

  double conversion;
  double convTarget;
  
  if(SIcSource == NULL || SIcTarget == NULL) return 0;
 
  if( !DvUnitSameBaseSI(SIcSource, SIcTarget) ) return 0;
  
  convTarget = DvUnitConvFactor( DvUnitGetBaseSI(SIcTarget).c_str() );
  if( convTarget < 1.0e-60 ) return 0; 
  
  // Conversion factor
  conversion = DvUnitConvFactor( DvUnitGetBaseSI(SIcSource).c_str() );
  conversion  /=  convTarget;

  return conversion;
  
} 

 std::string  DvUnitUnitsStr(const char * SIc){
   
     std::string newSIc(SIc);
	 size_t ptr = newSIc.find( ">");
	 if (ptr == std::string::npos) return std::string("");
	 
	 return newSIc.substr(ptr+1);

 }
 

