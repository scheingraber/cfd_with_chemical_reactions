#include <ctype.h>
#include <errno.h>
#include <stdio.h>
#include "helper.h"
#include "matrix.h"

/* ----------------------------------------------------------------------- */
/*                             auxiliary functions                         */
/* ----------------------------------------------------------------------- */

/* ----------------------------------------------------------------------- */
/*                         local auxiliary functions                       */
/* ----------------------------------------------------------------------- */

clock_t last_timer_reset;



/* ----------------------------------------------------------------------- */
/*                             read datafile                               */
/* ----------------------------------------------------------------------- */

void errhandler( int nLine, const char *szFile, const char *szString )
{
    int err = errno;

    fprintf( ERROUT, "%s:%d Error : %s", szFile, nLine, szString );
    fprintf( ERROUT, "\n" );
    
    /* if an error within the c-library occured, an error code can be   */
    /* found in the global variable err                                 */
    if( err != 0 )
    {
    fprintf( ERROUT, "C-Lib   errno    = %d\n", err);
    fprintf( ERROUT, "C-Lib   strerror = %s\n", strerror( err ) );
    }
    exit(1);
}


/*  for comfort */
#define READ_ERROR(szMessage, szVarName, szFileName, nLine) \
  { char szTmp[80]; \
    if( nLine ) \
    sprintf( szTmp, " %s  File: %s   Variable: %s  Line: %d", szMessage, szFileName, szVarName, nLine ); \
    else \
    sprintf( szTmp, " %s  File: %s   Variable: %s ", szMessage, szFileName, szVarName); \
    ERROR( szTmp ); \
  }
    

/* --------------------------------------------------------------------------*/
/* The function searches the datafile fh for the line defining the variable  */
/* szVarName and returns the respctive string including the value of the     */
/* variable. If there's no appropriate line within the datafile, the program */
/* stops with an error messsage.                                             */
/* ATTENTION: The pointer returned refers to a static variable within the    */
/* function. To maintain the string over several program calls, it has to be */
/* copied!!!                                                                 */
/*                                                                           */
char* find_string( const char* szFileName, const char *szVarName )
{ 
    int nLine = 0;
    unsigned int i;
    FILE *fh = NULL;
    
    static char szBuffer[MAX_LINE_LENGTH];  /* containes the line read  */
                                               /* from the datafile        */

    char* szLine = szBuffer;
    char* szValue = NULL;
    char* szName = NULL;

    /* open file */
    fh = fopen( szFileName, "rt" );
    if( fh == 0 ) 
    READ_ERROR("Could not open file", szVarName, szFileName, 0);

    /* searching */
    while( ! feof(fh) )
    {
    fgets( szLine, MAX_LINE_LENGTH, fh );
    ++nLine;

    /* remove comments */
    for( i = 0; i < strlen(szLine); i++)
        if( szLine[i] == '#' )
        {
        szLine[i] = '\0'; /* Stringende setzen */
        break;
        }

    /* remove empty lines */
    while( isspace( (int)*szLine ) && *szLine) ++szLine;
    if( strlen( szLine ) == 0) continue; 

    /* now, the name can be extracted */
    szName = szLine;
    szValue = szLine;
    while( (isalnum( (int)*szValue ) || *szValue == '_') && *szValue) ++szValue;
    
    /* is the value for the respective name missing? */
    if( *szValue == '\n' || strlen( szValue) == 0)  
        READ_ERROR("wrong format", szName, szFileName, nLine);
    
    *szValue = 0;       /* complete szName! at the right place */
    ++szValue;
        
    /* read next line if the correct name wasn't found */
    if( strcmp( szVarName, szName)) continue;

    /* remove all leading blnkets and tabs from the value string  */
    while( isspace( (int)*szValue) ) ++szValue;
    if( *szValue == '\n' || strlen( szValue) == 0)  
        READ_ERROR("wrong format", szName, szFileName, nLine);
    
    fclose(fh);
    return szValue;
    }  
   
    READ_ERROR("variable not found", szVarName, szFileName, nLine);
    
    return NULL;        /* dummy to satisfy the compiler  */
} 

void read_string( const char* szFileName, const char* szVarName, char*   pVariable)
{
    char* szValue = NULL;   /* string containg the read variable value */

    if( szVarName  == 0 )  ERROR("null pointer given as variable name" );
    if( szFileName == 0 )  ERROR("null pointer given as filename" );
    if( pVariable  == 0 )  ERROR("null pointer given as variable" );

    if( szVarName[0] == '*' )
    szValue = find_string( szFileName, szVarName +1 );
    else
    szValue = find_string( szFileName, szVarName );
    
    if( sscanf( szValue, "%s", pVariable) == 0)
    READ_ERROR("wrong format", szVarName, szFileName,0);

}

void read_int( const char* szFileName, const char* szVarName, int* pVariable)
{
    char* szValue = NULL;   /* string containing the read variable value */

    if( szVarName  == 0 )  ERROR("null pointer given as varable name" );
    if( szFileName == 0 )  ERROR("null pointer given as filename" );
    if( pVariable  == 0 )  ERROR("null pointer given as variable" );

    if( szVarName[0] == '*' )
    szValue = find_string( szFileName, szVarName +1 );
    else
    szValue = find_string( szFileName, szVarName );

    if( sscanf( szValue, "%d", pVariable) == 0)
    READ_ERROR("wrong format", szVarName, szFileName, 0);

}

void read_uint( const char* szFileName, const char* szVarName, unsigned int* pVariable)
{
    char* szValue = NULL;   /* string containing the read variable value */

    if( szVarName  == 0 )  ERROR("null pointer given as varable name" );
    if( szFileName == 0 )  ERROR("null pointer given as filename" );
    if( pVariable  == 0 )  ERROR("null pointer given as variable" );

    if( szVarName[0] == '*' )
    szValue = find_string( szFileName, szVarName +1 );
    else
    szValue = find_string( szFileName, szVarName );

    if( sscanf( szValue, "%u", pVariable) == 0)
    READ_ERROR("wrong format", szVarName, szFileName, 0);

}

void read_double( const char* szFileName, const char* szVarName, double* pVariable)
{
    char* szValue = NULL;   /* String mit dem eingelesenen Variablenwert */

    if( szVarName  == 0 )  ERROR("null pointer given as varable name" );
    if( szFileName == 0 )  ERROR("null pointer given as filename" );
    if( pVariable  == 0 )  ERROR("null pointer given as variable" );

    if( szVarName[0] == '*' )
    szValue = find_string( szFileName, szVarName +1 );
    else
    szValue = find_string( szFileName, szVarName );
    
    if( sscanf( szValue, "%lf", pVariable) == 0)
    READ_ERROR("wrong format", szVarName, szFileName, 0);

}


/* ----------------------------------------------------------------------- */
/*                   write matrices to a file                              */
/* ----------------------------------------------------------------------- */

void write_matrix( const char* szFileName,       /* filename */
           double **m,             /* matrix */
           int nrl,            /* first column */
           int nrh,            /* last column */
           int ncl,            /* first row */
           int nch,            /* last row */
         double xlength,           /* size of the geometry in */
                                               /* x-direction */
         double ylength,           /* size of the geometry in */
                                               /* y-direction  */
           int fFirst )            /* 0 == append, else overwrite*/
{
   int i, j;
   FILE * fh = 0;
   int nSize = (nrh-nrl+1) * (nch-ncl+1);
   float *tmp = (float *)malloc( (size_t)(nSize * sizeof(float)));
   int k = 0;

   if( fFirst )             /* first call of the function ? */
   {
       fh = fopen( szFileName, "w");    /* overwrite file/write new file */
       if( fh == NULL )         /* opening failed ? */
       {
       char szBuff[80];
       sprintf( szBuff, "Outputfile %s cannot be created", szFileName );
       ERROR( szBuff );
       }
       
/*       fprintf( fh,"%f\n%f\n%d\n%d\n%d\n%d\n", xlength, ylength, nrl, nrh, ncl, nch ); */
   }
   else
   {
       fh = fopen( szFileName ,"a");    /* append to the file */
       if( fh == NULL )         /* opening failed ? */
       {
       char szBuff[80];
       sprintf( szBuff, "Outputfile %s cannot be opened", szFileName );
       ERROR( szBuff );
       }
   } 

   for( j = ncl; j <= nch; j++)
       for( i = nrl; i <= nrh; i++)
       tmp[k++] = (float)m[i][j];

   fwrite( tmp, sizeof(float), nSize, fh);

   if( fclose(fh) )
   {
       char szBuff[80];
       sprintf( szBuff, "Outputfile %s cannot be closed", szFileName );
       ERROR( szBuff );
   };

   free( tmp );
}


void read_matrix( const char* szFileName,       /* filename */
           double **m,             /* matrix */
           int nrl,            /* first column */
           int nrh,            /* last column */
           int ncl,            /* first row */
           int nch             /* last row */
                  )       
{
   int i, j;
   FILE * fh = 0;
   int nSize = (nrh-nrl+1) * (nch-ncl+1);
   float *tmp = (float *)malloc( (size_t)(nSize * sizeof(float)));
   int k = 0;

       fh = fopen( szFileName, "r");    /* overwrite file/write new file */
       if( fh == NULL )         /* opening failed ? */
       {
       char szBuff[80];
       sprintf( szBuff, "Can not read file %s !!!", szFileName );
       ERROR( szBuff );
       }


   fread( tmp, sizeof(float), nSize, fh);

   for( j = ncl; j <= nch; j++)
       for( i = nrl; i <= nrh; i++)
       m[i][j]=tmp[k++];

   if( fclose(fh) )
   {
       char szBuff[80];
       /*orig bug:
       sscanf( szBuff, "Inputfile %s cannot be closed", szFileName );*/
       sprintf( szBuff, "Inputfile %s cannot be closed", szFileName );
       ERROR( szBuff );
   };

   free( tmp );
}


/* ----------------------------------------------------------------------- */
/*                      general matrix functions                           */
/* ----------------------------------------------------------------------- */


void init_matrix( double **m, int nrl, int nrh, int ncl, int nch, double a)
{
   int i,j;
   for( i = nrl; i <= nrh; i++)
       for( j = ncl; j <= nch; j++)
       m[i][j] = a;
}


void init_imatrix( int **m, int nrl, int nrh, int ncl, int nch, int a)
{
   int i,j;
   for( i = nrl; i <= nrh; i++)
       for( j = ncl; j <= nch; j++)
       m[i][j] = a;
}

