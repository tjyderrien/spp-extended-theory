/*
   =================================================================
  interspectus -- parser.c : used to read parameters from the input file
   =================================================================

   Copyright (C) 1998  Valerio Olevano
   Distributed under the Library General Public License version 2. 
   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330,
   Boston, MA 02111-1307, USA.
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define MAXLENSTR 80
#define MAXSTR 200
#define MAXLINELEN 200

char string[MAXSTR][MAXLENSTR];
int nstr;
char cdum[MAXLINELEN];

read_input_file_(char *filename, int *l) { read_input_file(filename,l);}
read_input_file(char *filename, int *l)
/* read a file or stdin writing it in string */
/* discard comments after the '#' symbol */
{
  FILE * fp;
  char *pc;
  char * fn;

  /* next 3lines to overcome the fortran problem to \0 terminate the strings */
  fn = (char *) calloc(*l + 1,sizeof(char));
  strncpy(fn,filename,*l);
  fn[*l] = '\0';

  nstr = 0;
  if(strcmp(fn,"stdin")) {
    fp = fopen(fn,"r");
    if(fp==NULL) { 
      fprintf(stderr,"error: input file '%s' not found!\n",fn);
      exit(2);
    }
  } 
  else fp = stdin;

  while(fscanf(fp,"%s",&string[nstr][0])!=EOF) { /* read strings */
    if((pc=strchr(&string[nstr][0],'#'))==NULL) nstr++; /* find comments */
    else {
      if(pc!=&string[nstr][0]) {*pc = '\0'; nstr++;}
      if(fgets(cdum,MAXLINELEN,fp)==NULL) { /* skip rest of line */
        if(strcmp(fn,"stdin")) fclose(fp);
        free(fn);
        return(1);
      }
    }
    if(nstr>=MAXSTR) fprintf(stderr,"error: too many strings."), exit(2);
  }  

  if(strcmp(fn,"stdin")) fclose(fp);
  /* printf("file composition:\n");
     for(i=0;i<nstr;i++) printf("%d %s\n",i,&string[i][0]); */
  free(fn);
  return(1);
}


read_string_(char *keyword, int *l) { read_string(keyword,l);}
read_string(char *keyword, int *l)
/* parse a given keyword returning its position (the last if many) */
/* return 0 if not found */
{
  int i;
  char *k;
  k = (char *) calloc(*l + 1,sizeof(char));
  strncpy(k,keyword,*l);
  k[*l] = '\0';
  i = nstr-1;
  while(strcasecmp(k,&string[i][0]) && i>=0) i--;
  free(k);
  return(i+1);
}

read_int_(char *keyword, int *l, int *ikey) { read_int(keyword,l,ikey);}
read_int(char *keyword, int *l, int *ikey)
/* parse a keyword and a subsequent int (put in ikey) */
/* return 1 if found, 0 if not */
{
  int i;
  long int il;
  char *endptr;
  char *k;
  k = (char *) calloc(*l + 1,sizeof(char));
  strncpy(k,keyword,*l);
  k[*l] = '\0';
  i = nstr-2; /* does not consider string in last position */
  while(i>=0) {
    while(strcasecmp(k,&string[i][0]) && i>=0) i--;
    if(i<0) {free(k); return(0);}
    il = strtol(&string[i+1][0],&endptr,0);
    if(endptr==&string[i+1][0] && il==0) {i--; continue;}
    *ikey = (int) il;
    free(k);
    return(1);
  }
  free(k);
  return(0);
}

read_float_(char *keyword, int *l, float *fkey) {read_float(keyword,l,fkey);}
read_float(char *keyword,int *l, float*fkey)
/* parse a keyword and a subsequent float (put in fkey) */
/* return 1 if found, 0 if not */
{
  int i;
  float f;
  float strtof();
  char *endptr;
  char *k;
  k = (char *) calloc(*l + 1,sizeof(char));
  strncpy(k,keyword,*l);
  k[*l] = '\0';
  i = nstr-2; /* does not consider string in last position */
  while(i>=0) {
    while(strcasecmp(k,&string[i][0]) && i>=0) i--;
    if(i<0) {free(k); return(0);}
    f = strtof(&string[i+1][0],&endptr);
    if(endptr==&string[i+1][0] && f==0.0) {i--; continue;}
    *fkey = f;
    free(k);
    return(1);
  }
  free(k);
  return(0);
}

read_float_array_(char *keyword, int *l, float *fkey, int* pnf) 
  {read_float_array(keyword,l,fkey,pnf);}
read_float_array(char *keyword, int *l, float *fkey, int *pnf)
/* parse a keyword and a subsequent float array of *pnf elements (put in fkey)*/
/* return 1 if found, 0 if not */
{
  int i,j;
  float *f;
  float strtof();
  char *endptr;
  char *k;
  f = (float *) calloc(*pnf,sizeof(float));
  k = (char *) calloc(*l + 1,sizeof(char));
  strncpy(k,keyword,*l);
  k[*l] = '\0';
  i = nstr - *pnf - 1; /* does not consider string in last position */
  while(i>=0) {
    while(strcasecmp(k,&string[i][0]) && i>=0) i--;
    if(i<0) {free(k); return(0);}
    for(j=0;j<*pnf;j++) {
      f[j] = strtof(&string[i+1+j][0],&endptr);
      if(endptr==&string[i+1+j][0] && f[j]==0.0) break;
    }
    if(j<*pnf) {i--; continue;}
    for(j=0;j<*pnf;j++,fkey++) *fkey = f[j];
    free(k);
    return(1);
  }
  free(k);
  free(f);
  return(0);
}
