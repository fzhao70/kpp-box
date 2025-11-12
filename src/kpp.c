/******************************************************************************

  KPP - The Kinetic PreProcessor
        Builds simulation code for chemical kinetic systems

  Copyright (C) 1995-1996 Valeriu Damian and Adrian Sandu
  Copyright (C) 1997-2005 Adrian Sandu

  KPP is free software; you can redistribute it and/or modify it under the
  terms of the GNU General Public License as published by the Free Software
  Foundation (http://www.gnu.org/copyleft/gpl.html); either version 2 of the
  License, or (at your option) any later version.

  KPP is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
  details.

  You should have received a copy of the GNU General Public License along
  with this program; if not, consult http://www.gnu.org/copyleft/gpl.html or
  write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330,
  Boston, MA  02111-1307,  USA.

  Adrian Sandu
  Computer Science Department
  Virginia Polytechnic Institute and State University
  Blacksburg, VA 24060
  E-mail: sandu@cs.vt.edu

******************************************************************************/


#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include <sys/stat.h>
#include "gdata.h"
#include "scan.h"

char *eqFileName;
char *rootFileName = "ff";
char Home[ MAX_PATH ] = ""; 

short int linStru[ MAX_SPECIES ];
short int colStru[ MAX_SPECIES ];
short int bestStru[ MAX_SPECIES ];
short int *Stru;

enum stru_criteria { UNSORT, LINSORT, COLSORT, BESTSORT };

void EqCopy( EQ_VECT e1, EQ_VECT e2 )
{
int i;

  for( i = 0; i < EqnNr; i++ ) e2[i] = e1[i];
}

int NoSort( const void *p1, const void *p2 )
{
  return -1;
}

int CodeCmp( const void *p1, const void *p2 )
{
CODE *c1, *c2;

  c1 = (CODE*)p1;  
  c2 = (CODE*)p2;  
  
  if ( *c1 < *c2 ) return -1;
  if ( *c1 > *c2 ) return 1;
  return 0;
}

int CodeRCmp( const void *p1, const void *p2 )
{
int rc1, rc2;
CODE *c1, *c2;

  c1 = (CODE*)p1;  
  c2 = (CODE*)p2;  

  rc1 = Reactive[ ReverseCode[ *c1 ] ]; 
  rc2 = Reactive[ ReverseCode[ *c2 ] ];
  if ( rc1 > rc2 ) return -1; 
  if ( rc1 < rc2 ) return 1;
  if ( *c1 < *c2 ) return -1;
  if ( *c1 > *c2 ) return 1;
  return 0;
}

int CodeSSCmp( const void *p1, const void *p2 )
{
  return -CodeRCmp(p1,p2);
}

int CodeSCmp( const void *p1, const void *p2 )
{
CODE *c1, *c2;
short int sc1, sc2;

  c1 = (CODE*)p1;  
  c2 = (CODE*)p2;  
  
  sc1 = Stru[ ReverseCode[ *c1 ] ];
  sc2 = Stru[ ReverseCode[ *c2 ] ];
  
  if ( sc1 > sc2 ) return 1; 
  if ( sc1 < sc2 ) return -1;
  if ( *c1 < *c2 ) return 1;
  if ( *c1 > *c2 ) return -1;
  return 0;
}

void UpdateStructJ()
{
int i,j,k;

  for ( i=0; i<VarNr; i++ )
    for ( j=0; j<VarNr; j++ )
      structJ[i][j]=(i==j)?1:0;
              
  for (i = 0; i < VarNr; i++)
    for (j = 0; j < VarNr; j++)
      for (k = 0; k < EqnNr; k++)
        if ( Stoich[i][k]*((Stoich_Left[j][k])?1:0) != 0.0 )
          structJ[i][j]=1;

  for ( i=0; i<VarNr; i++ )
    for ( j=0; j<VarNr; j++ )
      LUstructJ[i][j]=structJ[i][j];
              
}

int ComputeLUStructJ()
{
int i,j,k;
int nu,nl;

  for (j = 0; j < VarNr-1; j++) {
    for (i = j+1; i < VarNr; i++) {
      if( LUstructJ[i][j] ) {
        for (k = j; k < VarNr; k++) 
         /*   LUstructJ[i][k] += LUstructJ[j][k]; */
	   if ( LUstructJ[j][k] != 0 )
              LUstructJ[i][k] = 1; 
      }   
    }
  }  

  
  nu = 0; nl = 0;
  for (i = 0; i < VarNr; i++) 
    for (j = 0; j < VarNr; j++)
      if( LUstructJ[i][j] ) {
        if(i > j) nl++;
        if(i <= j) nu++;     
      }          

  return nu+nl;
}

int LUnonZero()
{
CODE v[MAX_SPECIES];
CODE *var;
int i,j,k;
int nu,nl;

  var = v;
  if( Stru != bestStru ) {
    for( i=0; i<VarNr; i++ )
      var[i] = Code[i];
    qsort( (void*)var, VarNr, sizeof(CODE), CodeSCmp );
  } else {
    var = bestStru;
  }
          
  for (i = 0; i < VarNr; i++)
    for (j = 0; j < VarNr; j++)
      LUstructJ[i][j] = structJ[ ReverseCode[var[i]] ][ ReverseCode[var[j]] ];  

  return ComputeLUStructJ();
}

void LinColSparsity()
{
int i,j,k;
int nlin, ncol;
FILE * fff;

  for ( i=0; i<VarNr; i++ )
    for ( j=0; j<VarNr; j++ )
      structJ[i][j]=(i==j)?1:0;
              
  for (i = 0; i < VarNr; i++)
    for (j = 0; j < VarNr; j++)
      for (k = 0; k < EqnNr; k++)
        if ( Stoich[i][k]*((Stoich_Left[j][k])?1:0) != 0.0 )
          structJ[i][j]=1;

  for ( i=0; i<VarNr; i++ ) {
    linStru[i] = 0;
    for (j = 0; j < VarNr; j++)
      linStru[i] += structJ[i][j];
  }

  for ( i=0; i<VarNr; i++ ) {
    colStru[i] = 0;
    for (j = 0; j < VarNr; j++)
      colStru[i] += structJ[j][i];
    colStru[i] *= linStru[i];  
  }

  Stru = linStru;
  nlin = LUnonZero();
  Stru = colStru;
  ncol = LUnonZero();
  if( nlin <= ncol ) { 
    Stru = linStru;
    LUnonZero();
  } 
}

void BestSparsity()
{
int i,j,k;
int cnz, lnz;
int best, crt;
int best_i;
int tmp;
int s;

  UpdateStructJ();
      
  for ( i=0; i<VarNr; i++ )
    bestStru[i] = Code[i];

  for ( s=0; s<VarNr-1; s++ ) {
    best = MAX_SPECIES*MAX_SPECIES; best_i = 0;
    for ( i=s; i<VarNr; i++ ) {
      cnz = 0;lnz = 0;
      for (j = s; j < VarNr; j++) {
        cnz += (LUstructJ[i][j]?1:0);
        lnz += (LUstructJ[j][i]?1:0);
      }
      crt = (cnz-1)*(lnz-1);
      if( crt < best ) {
        best = crt;
        best_i = i;
      }    
    }
    for ( i=0; i<VarNr; i++ ) {
      tmp = LUstructJ[s][i];
      LUstructJ[s][i] = LUstructJ[best_i][i];
      LUstructJ[best_i][i] = tmp;      
    }   
    for ( i=0; i<VarNr; i++ ) {
      tmp = LUstructJ[i][s];
      LUstructJ[i][s] = LUstructJ[i][best_i];
      LUstructJ[i][best_i] = tmp;      
    }   
    tmp = bestStru[s];
    bestStru[s] = bestStru[best_i];
    bestStru[best_i] = tmp;
    
    for (i = s+1; i < VarNr; i++) {
      if( LUstructJ[i][s] ) {
        for (k = s; k < VarNr; k++) 
          LUstructJ[i][k] += LUstructJ[s][k];
      }   
    }
  }

  Stru = bestStru;
}

void ReorderSpecies( int criteria )
{
CODE *var;
CODE *fix;
CODE *dummy;
EQ_VECT *tmpStoich_Left;
EQ_VECT *tmpStoich_Right;
EQ_VECT *tmpStoich;
CODE *tmpCode;
CODE *tmpReact;
int i, k;
int new;
int (*cmpVar)(const void *, const void *);
int (*cmpFix)(const void *, const void *);
int dummyNr;

  cmpVar = CodeCmp;
  cmpFix = CodeCmp;

  switch( criteria ) {
    case UNSORT:   cmpVar = useJacobian ? CodeRCmp : CodeCmp;
                   break;  
    case LINSORT:  cmpVar = useJacobian ? CodeSCmp : CodeCmp;
                   Stru = linStru;  
                   break;  
    case COLSORT:  cmpVar = useJacobian ? CodeSCmp : CodeCmp;
                   Stru = colStru;
                   break;  
    case BESTSORT: cmpVar = useJacobian ? NoSort : CodeCmp;
                   break;  
  }

  VarNr = 0;
  VarActiveNr = 0;
  FixNr = 0;
  dummyNr = 0;

  var = (CODE*)malloc( SpcNr * sizeof(CODE) );
  fix = (CODE*)malloc( SpcNr * sizeof(CODE) );
  dummy = (CODE*)malloc( 5 * sizeof(CODE) );
  tmpStoich_Left = (EQ_VECT*)malloc( SpcNr * sizeof(EQ_VECT) );
  tmpStoich_Right = (EQ_VECT*)malloc( SpcNr * sizeof(EQ_VECT) );
  tmpStoich = (EQ_VECT*)malloc( SpcNr * sizeof(EQ_VECT) );
  tmpCode = (CODE*)malloc( SpcNr * sizeof(CODE) );
  tmpReact = (CODE*)malloc( SpcNr * sizeof(CODE) );

  for( i = 0; i < SpcNr; i++ ) {
    switch( SpeciesTable[ Code[i] ].type ) {
      case VAR_SPC:  var[ VarNr++ ] = Code[ i ];
		     break;
      case FIX_SPC:  fix[ FixNr++ ] = Code[ i ]; 
		     break;
      case DUMMY_SPC:dummy[ dummyNr++ ] = Code[ i ]; 
		     break;
    }
  }

  if( Stru != bestStru ) {
    qsort( (void*)var, VarNr, sizeof(CODE), cmpVar );
  } else {
    for( i = 0; i < SpcNr; i++ )
      var[i] = bestStru[i];
  }  
  qsort( (void*)fix, FixNr, sizeof(CODE), cmpFix );

  for( i = 0; i < SpcNr; i++ ) {
    EqCopy( Stoich_Left[i], tmpStoich_Left[i] );
    EqCopy( Stoich_Right[i], tmpStoich_Right[i] );
    EqCopy( Stoich[i], tmpStoich[i] );
    tmpCode[i] = Code[i];
    tmpReact[i] = Reactive[i];
  }

  SpcNr -= dummyNr;
  dummyNr = 0;

  k = 0;
  for( i = 0; i < VarNr; i++ ) {
    new = ReverseCode[ var[i] ];
    EqCopy( tmpStoich_Left[ new ], Stoich_Left[ k ] );
    EqCopy( tmpStoich_Right[ new ], Stoich_Right[ k ] );
    EqCopy( tmpStoich[ new ], Stoich[ k ] );
    Code[ k ] = tmpCode[ new ];
    Reactive[ k ] = tmpReact[ new ];
    if( Reactive[ k ] ) VarActiveNr++; 
    k++;
  }
  for( i = 0; i < FixNr; i++ ) {
    new = ReverseCode[ fix[i] ];
    EqCopy( tmpStoich_Left[ new ], Stoich_Left[ k ] );
    EqCopy( tmpStoich_Right[ new ], Stoich_Right[ k ] );
    EqCopy( tmpStoich[ new ], Stoich[ k ] );
    Code[ k ] = tmpCode[ new ];
    Reactive[ k ] = tmpReact[ new ];
    k++;
  }
  for( i = 0; i < dummyNr; i++ ) {
    new = ReverseCode[ dummy[i] ];
    EqCopy( tmpStoich_Left[ new ], Stoich_Left[ k ] );
    EqCopy( tmpStoich_Right[ new ], Stoich_Right[ k ] );
    EqCopy( tmpStoich[ new ], Stoich[ k ] );
    Code[ k ] = tmpCode[ new ];
    Reactive[ k ] = tmpReact[ new ];
    k++;
  }


  for( i = 0; i < SpcNr+dummyNr; i++ )
    ReverseCode[ Code[i] ] = i;

  free( tmpReact );
  free( tmpCode );
  free( tmpStoich );
  free( tmpStoich_Right );
  free( tmpStoich_Left );
  free( dummy );
  free( fix );
  free( var );   

  fflush(stdout);
}

/* Allocate Internal Arrays */
void  AllocInternalArrays( void )
{
int i;

if ( (Stoich_Left =(float**)calloc(MAX_SPECIES,sizeof(float*)))==NULL ) 
    FatalError(-30,"Cannot allocate Stoich_Left.\n");

for (i=0; i<MAX_SPECIES; i++)    
    if ( (Stoich_Left[i] = (float*)calloc(MAX_EQN,sizeof(float)))==NULL ) {
        FatalError(-30,"Cannot allocate Stoich_Left[%d]",i,MAX_SPECIES);
    }

if ( (Stoich_Right = (float**)calloc(MAX_SPECIES,sizeof(float*)))==NULL ) 
    FatalError(-30,"Cannot allocate Stoich_Right.\n");

for (i=0; i<MAX_SPECIES; i++)    
    if ( (Stoich_Right[i] = (float*)calloc(MAX_EQN,sizeof(float)))==NULL ) {
        FatalError(-30,"Cannot allocate Stoich_Right[%d].",i);
    }

if ( (Stoich = (float**)calloc(MAX_SPECIES,sizeof(float*)))==NULL ) 
    FatalError(-30,"Cannot allocate Stoich.\n");

for (i=0; i<MAX_SPECIES; i++)    
    if ( (Stoich[i] = (float*)calloc(MAX_EQN,sizeof(float)))==NULL ) {
        FatalError(-30,"Cannot allocate Stoich[%d].",i);
    }

}


/* Allocate Structure Arrays */
void  AllocStructArrays( void )
{
int i;


if ( (structB = (int**)calloc(EqnNr,sizeof(int*)))==NULL ) 
    FatalError(-30, "Cannot allocate structB.");

for (i=0; i<EqnNr; i++)    
    if ( (structB[i] =(int*) calloc(SpcNr,sizeof(int)))==NULL )
        FatalError(-30, "Cannot allocate structB[%d].\n",i);
    
if ( (structJ = (int**)calloc(SpcNr,sizeof(int*)))==NULL ) 
    FatalError(-30, "Cannot allocate structJ.");

for (i=0; i<SpcNr; i++)    
    if ( (structJ[i] =(int*) calloc(SpcNr,sizeof(int)))==NULL ) 
        FatalError(-30, "Cannot allocate structJ[%d].\n",i);
    
if ( (LUstructJ = (int**)calloc(SpcNr,sizeof(int*)))==NULL ) 
    FatalError(-30, "Cannot allocate LUstructJ.");

for (i=0; i<SpcNr; i++)    
    if ( (LUstructJ[i] = (int*)calloc(SpcNr,sizeof(int)))==NULL ) 
        FatalError(-30, "Cannot allocate LUstructJ[%d].\n",i);

}

/*******************************************************************/                    
int Postprocess( char * root )
{
char buf[ 200 ];
char cmd[500];
char cmdexe[500];
static char tmpfile[] = "kppfile.tmp";
FILE * fp;

  if ( useLang == MATLAB_LANG ) {
 /*  Add rate function definitions as internal functions to the Update_RCONST file*/
    sprintf( buf, "cat %s_Update_RCONST.m %s_Rates.m > tmp; mv tmp %s_Update_RCONST.m;", 
        root, root, root ); 
    system( buf );	 
  }

/*    Postprocessing to replace parameter names by values in the declarations
  strcpy( cmd, "sed " );
  sprintf( cmd, "%s -e 's/(NVAR)/(%d)/g'", cmd, VarNr );  
  sprintf( cmd, "%s -e 's/(NFIX)/(%d)/g'", cmd, FixNr );  
  sprintf( cmd, "%s -e 's/(NSPEC)/(%d)/g'", cmd,SpcNr );  
  sprintf( cmd, "%s -e 's/(NREACT)/(%d)/g'", cmd, EqnNr );  
  sprintf( cmd, "%s -e 's/(NONZERO)/(%d)/g'", cmd, Jac_NZ );  
  sprintf( cmd, "%s -e 's/(LU_NONZERO)/(%d)/g'", cmd, LU_Jac_NZ );  
  sprintf( cmd, "%s -e 's/(NHESS)/(%)/g'", cmd, Hess_NZ );  
   
  sprintf( buf, "%s_Function", rootFileName );  
  switch( useLang ) { 
    case F77_LANG: sprintf( buf, "%s.f", buf );
                 break;
    case F90_LANG: sprintf( buf, "%s.f90", buf );
                 break;
    case C_LANG: sprintf( buf, "%s.c", buf );
                 break;
    case MATLAB_LANG: sprintf( buf, "%s.m", buf );
                 break;
    default: printf("\n Language '%d' not implemented!\n",useLang); 
                 exit(1);
  }
  sprintf( cmdexe, "%s %s > %s; mv %s %s;", cmd, buf, tmpfile, tmpfile, buf );  
  printf("\n\nCMDEXE='%s'\n",cmdexe);
  system( cmdexe );
*/
}

/*******************************************************************/
void PrintUsage( void )
{
  printf("\nKPP - The Kinetic PreProcessor version %s\n", KPP_VERSION);
  printf("Builds simulation code for chemical kinetic systems\n\n");

  printf("Usage:\n");
  printf("  kpp <equations_file> [output_file]     Generate code from KPP input file\n");
  printf("  kpp --help, -h                         Show this help message\n");
  printf("  kpp --version, -v                      Show version information\n");
  printf("  kpp --list-models                      List available models\n");
  printf("  kpp --list-integrators                 List available integrators\n");
  printf("  kpp --list-drivers                     List available drivers\n\n");

  printf("Examples:\n");
  printf("  kpp small.kpp                          Process small.kpp, output to small.*\n");
  printf("  kpp mymodel.kpp output                 Process mymodel.kpp, output to output.*\n\n");

  printf("Environment:\n");
  printf("  KPP_HOME    Path to KPP installation directory (currently: %s)\n\n",
         Home[0] ? Home : "not set");

  printf("For more information, see the KPP User Manual.\n\n");
}

/*******************************************************************/
void PrintVersion( void )
{
  printf("KPP version %s\n", KPP_VERSION);
  printf("Copyright (C) 1995-1996 Valeriu Damian and Adrian Sandu\n");
  printf("Copyright (C) 1997-2005 Adrian Sandu\n");
  printf("License: GPL v2 or later\n\n");
}

/*******************************************************************/
void ListDirectory( const char* label, const char* path, const char* extension )
{
  DIR *dir;
  struct dirent *entry;
  int count = 0;
  char fullpath[MAX_PATH];

  if (!path || !path[0]) {
    printf("Error: %s path not set (check KPP_HOME environment variable)\n", label);
    return;
  }

  dir = opendir(path);
  if (!dir) {
    printf("Error: Cannot open %s directory: %s\n", label, path);
    return;
  }

  printf("\nAvailable %s (in %s):\n", label, path);
  printf("-----------------------------------------------\n");

  while ((entry = readdir(dir)) != NULL) {
    if (entry->d_name[0] == '.') continue;

    if (extension && extension[0]) {
      char *dot = strrchr(entry->d_name, '.');
      if (!dot || strcmp(dot, extension) != 0) continue;

      /* Print name without extension */
      int len = dot - entry->d_name;
      printf("  %.*s\n", len, entry->d_name);
      count++;
    } else {
      /* Check if it's a regular file */
      snprintf(fullpath, MAX_PATH, "%s/%s", path, entry->d_name);
      struct stat st;
      if (stat(fullpath, &st) == 0 && S_ISREG(st.st_mode)) {
        printf("  %s\n", entry->d_name);
        count++;
      }
    }
  }

  closedir(dir);

  if (count == 0) {
    printf("  (none found)\n");
  } else {
    printf("\nTotal: %d\n", count);
  }
}

/*******************************************************************/
void ListModels( void )
{
  char path[MAX_PATH];
  DIR *dir;

  if (!Home[0]) {
    printf("Error: KPP_HOME environment variable not set.\n");
    printf("Please set KPP_HOME to your KPP installation directory.\n");
    return;
  }

  /* Try demo/models first, then models */
  snprintf(path, MAX_PATH, "%s/demo/models", Home);
  dir = opendir(path);
  if (!dir) {
    snprintf(path, MAX_PATH, "%s/models", Home);
  } else {
    closedir(dir);
  }

  ListDirectory("Models", path, ".def");

  printf("\nTo use a model, add this line to your .kpp file:\n");
  printf("  #INCLUDE <model_name>\n\n");
}

/*******************************************************************/
void ListIntegrators( void )
{
  char path[MAX_PATH];

  if (!Home[0]) {
    printf("Error: KPP_HOME environment variable not set.\n");
    printf("Please set KPP_HOME to your KPP installation directory.\n");
    return;
  }

  snprintf(path, MAX_PATH, "%s/int", Home);
  ListDirectory("Integrators", path, ".def");

  printf("\nTo use an integrator, add this line to your .kpp file:\n");
  printf("  #INTEGRATOR <integrator_name>\n\n");
}

/*******************************************************************/
void ListDrivers( void )
{
  char path[MAX_PATH];

  if (!Home[0]) {
    printf("Error: KPP_HOME environment variable not set.\n");
    printf("Please set KPP_HOME to your KPP installation directory.\n");
    return;
  }

  snprintf(path, MAX_PATH, "%s/drv", Home);
  ListDirectory("Drivers", path, NULL);

  printf("\nTo use a driver, add this line to your .kpp file:\n");
  printf("  #DRIVER <driver_name>\n\n");
}

/*******************************************************************/
int main( int argc, char * argv[] )
{
int status;
char name[ 200 ];
char *p;
int i,j;

  AllocInternalArrays();

  p = getenv("KPP_HOME");
  if( p ) strcpy( Home, p );

  /* Handle command-line flags */
  if( argc == 2 ) {
    if( strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0 ) {
      PrintUsage();
      exit(0);
    }
    if( strcmp(argv[1], "--version") == 0 || strcmp(argv[1], "-v") == 0 ) {
      PrintVersion();
      exit(0);
    }
    if( strcmp(argv[1], "--list-models") == 0 ) {
      ListModels();
      exit(0);
    }
    if( strcmp(argv[1], "--list-integrators") == 0 ) {
      ListIntegrators();
      exit(0);
    }
    if( strcmp(argv[1], "--list-drivers") == 0 ) {
      ListDrivers();
      exit(0);
    }
  }

  /* Handle normal file processing */
  switch( argc ) {
    case 3: eqFileName = argv[1];
            rootFileName = argv[2];
            break;
    case 2: eqFileName = argv[1];
            strcpy( name, eqFileName );
            p = name + strlen(name);
            while( p > name ) {
              if( *p == '.') {
                *p = '\0';
                break;
              }
              p--;
            }
	    rootFileName = name;
	    break;
    case 1: printf("\nError: Missing required argument <equations_file>\n");
            PrintUsage();
            exit(1);
    default: printf("\nError: Too many arguments\n");
            PrintUsage();
            exit(1);
  }

  /* Validate KPP_HOME is set */
  if (!Home[0]) {
    printf("\n%s\n", "========================================");
    printf("%s\n", "ERROR: KPP_HOME not set");
    printf("%s\n", "========================================");
    printf("\nThe KPP_HOME environment variable is not set.\n");
    printf("KPP needs this to locate models, integrators, and drivers.\n\n");
    printf("To fix this, either:\n");
    printf("  1. Run the installation script:\n");
    printf("       cd %s\n", getenv("PWD") ? getenv("PWD") : ".");
    printf("       ./install.sh\n\n");
    printf("  2. Or manually set KPP_HOME:\n");
    printf("       export KPP_HOME=/path/to/kpp-box\n");
    printf("       export PATH=$PATH:$KPP_HOME/bin\n\n");
    exit(1);
  }

  /* Validate input file exists */
  FILE *testfile = fopen(eqFileName, "r");
  if (!testfile) {
    printf("\n%s\n", "========================================");
    printf("ERROR: Cannot open input file\n");
    printf("%s\n", "========================================");
    printf("\nThe file '%s' does not exist or cannot be read.\n\n", eqFileName);
    printf("Please check:\n");
    printf("  - File path is correct\n");
    printf("  - File has .kpp extension\n");
    printf("  - You have read permissions\n\n");
    printf("Example files are available in:\n");
    printf("  %s/demo/examples/\n\n", Home);
    exit(1);
  }
  fclose(testfile);

  printf("\n========================================\n");
  printf("KPP - Kinetic PreProcessor v%s\n", KPP_VERSION);
  printf("========================================\n");

  printf("\n[1/4] Parsing equation file '%s'...\n", eqFileName);
  status = ParseEquationFile( argv[1] );

  if( status ) FatalError(2,"%d errors and %d warnings encountered.",
                           nError, nWarning );
  printf("      Success: Parsed %d species and %d reactions\n", SpcNr, EqnNr);

  /* Allocate some internal data structures */
  AllocStructArrays();

  printf("\n[2/4] Computing Jacobian sparsity structure...\n");
  ReorderSpecies( UNSORT );
  if (useReorder==1){
    BestSparsity();
    ReorderSpecies( BESTSORT );
    printf("      Optimized variable ordering for sparsity\n");
    }
  UpdateStructJ();
  ComputeLUStructJ();
  printf("      Success: %d variable species, %d fixed species\n", VarNr, FixNr);

  if( initNr == -1 ) initNr = VarNr;


  printf("\n[3/4] Generating code files...\n");
  Generate( rootFileName );
  printf("      Generated code files for model '%s'\n", rootFileName);

  printf("\n[4/4] Post-processing...\n");
  Postprocess( rootFileName );
  printf("      Post-processing complete\n");

  printf("\n========================================\n");
  printf("SUCCESS: Model '%s' created successfully!\n", rootFileName);
  printf("========================================\n\n");

  if( nError ) exit(4);
  if( nWarning ) exit(5);
  
  exit(0);
}
