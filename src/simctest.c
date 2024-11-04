#include<R.h>
#include <Rdefines.h> 

SEXP extendbounds(SEXP n, SEXP level, SEXP U, SEXP L, SEXP porig, SEXP preverr, SEXP maxerr, SEXP returnstopprob){
  SEXP Unew, Lnew, list, list_names;
  int nint;
  double* p0;
  double alphaint, oneminusalphaint;
  char *names[] = {"U", "L", "porig","preverr", "n","Sn","prob"}; // 7 instead of 4 because of stopprob
  int i;

  PROTECT(n = AS_INTEGER(n));
  PROTECT(level = AS_NUMERIC(level));
  nint=*INTEGER_POINTER(n);

  alphaint=*NUMERIC_POINTER(level);
  oneminusalphaint=1-alphaint;
  

  PROTECT(Unew=NEW_INTEGER(nint));
  int* Uptr=INTEGER_POINTER(Unew);

  PROTECT(Lnew=NEW_INTEGER(nint));
  int* Lptr=INTEGER_POINTER(Lnew);

  int p0size=length(porig)+10;//ceil(sqrt(nint*log(nint)));
  p0=  R_Calloc(p0size, double);
  if (p0size==10)  p0[0]=1;
  else
    for (i=0;i<length(porig);i++){
      p0[i]=NUMERIC_POINTER(porig)[i];
    }
  double errL=NUMERIC_POINTER(preverr)[0];
  double errU=NUMERIC_POINTER(preverr)[1];


  PROTECT(returnstopprob = AS_INTEGER(returnstopprob));
  int dostopprob=*INTEGER_POINTER(returnstopprob);
  
  // Begin: for stopprob
  int stopsize=0;
  int nstop=0;
  int* stopn=NULL;
  int* stopSn=NULL;
  double* stopProb=NULL;
  if (dostopprob){
    stopsize = nint+p0size;
    nstop=0;
    stopn=R_Calloc(stopsize,int);
    stopSn=R_Calloc(stopsize,int);
    stopProb=R_Calloc(stopsize,double);
  }
  // End: for stropprob

  PROTECT(L=AS_INTEGER(L));
  int p0start=INTEGER(L)[0]+1;
  PROTECT(U=AS_INTEGER(U));
  int Uakt=INTEGER(U)[0]-p0start;
  int Lakt=INTEGER(L)[0]-p0start;
  double *pos, *pos2;


  for (i=0; i<nint; i++){
    if (Uakt+1>p0size){
      //increase memory
      p0size*=2;
      p0=R_Realloc(p0, p0size,double);
    }
    p0[Uakt]=p0[Uakt-1]*alphaint;
    for (pos=p0+(Uakt-1),pos2=p0+Uakt-2; pos>p0; pos--,pos2--){
      *pos=(*pos)*oneminusalphaint+(*pos2)*alphaint;
    }
    p0[Lakt+1]*=oneminusalphaint;
    double maxerrakt=NUMERIC_POINTER(maxerr)[i];
    while(p0[Uakt]+errU<=maxerrakt){
      errU+=p0[Uakt];

      // Begin: for stopprob
      if (dostopprob){
	if (stopsize<=nstop){
	  //increase memory
	  stopsize*=2;
	  stopn=R_Realloc(stopn,stopsize,int);
	  stopSn=R_Realloc(stopSn,stopsize,int);
	  stopProb=R_Realloc(stopProb,stopsize,double);
	}
	stopn[nstop]=i;
	stopSn[nstop]=Uakt+p0start;
	stopProb[nstop]=p0[Uakt];
	nstop++;
      }
      // End: for stopprob

      Uakt--;      
    }
    while(p0[Lakt+1]+errL<=maxerrakt){
      errL+=p0[Lakt+1];

      // Begin: for stopprob
      if (dostopprob){
	if (stopsize<=nstop){                                       
	  //increase memory
	  stopsize*=2;
	  stopn=R_Realloc(stopn,stopsize,int);
	  stopSn=R_Realloc(stopSn,stopsize,int);
	  stopProb=R_Realloc(stopProb,stopsize,double);
	}
	stopn[nstop]=i;
	stopSn[nstop]=Lakt+p0start+1;
	stopProb[nstop]=p0[Lakt+1];
	nstop++;
      }
      // End: for stopprob

      Lakt++;
    }
    Uakt++;   
    Lptr[i]=Lakt+p0start;
    Uptr[i]=Uakt+p0start;
    if (Lakt>=0){
      //shift everything in the array downward     
      for (pos=p0+Lakt+1,pos2=p0;pos<p0+Uakt;pos++,pos2++){
	*pos2=*pos;
      }
      p0start+=Lakt+1;
      Uakt-=Lakt+1;
      Lakt=-1;
    }

  }
  PROTECT(list = allocVector(VECSXP, 7)); 
  SET_VECTOR_ELT(list, 0, Unew);
  SET_VECTOR_ELT(list, 1, Lnew);
  SEXP porgnew;
  PROTECT(porgnew=NEW_NUMERIC(Uakt));
  for (i=0;i<Uakt;i++)
    NUMERIC_POINTER(porgnew)[i]=p0[i];
  SET_VECTOR_ELT(list, 2, porgnew);

  SEXP preverrnew;
  PROTECT(preverrnew=NEW_NUMERIC(2));
  NUMERIC_POINTER(preverrnew)[0]=errL;
  NUMERIC_POINTER(preverrnew)[1]=errU;
  SET_VECTOR_ELT(list, 3, preverrnew);


  // Begin: for stopprob
  if (dostopprob){
    SEXP res_stopn, res_stopSn, res_stopProb;
    PROTECT(res_stopn=NEW_NUMERIC(nstop));
    PROTECT(res_stopSn=NEW_NUMERIC(nstop));
    PROTECT(res_stopProb=NEW_NUMERIC(nstop));
    for (i=0;i<nstop;i++){
      NUMERIC_POINTER(res_stopn)[i]=stopn[i];
      NUMERIC_POINTER(res_stopSn)[i]=stopSn[i];
      NUMERIC_POINTER(res_stopProb)[i]=stopProb[i];
    }
    SET_VECTOR_ELT(list, 4, res_stopn);
    SET_VECTOR_ELT(list, 5, res_stopSn);
    SET_VECTOR_ELT(list, 6, res_stopProb);
  }
  // End: for stopprob


  PROTECT(list_names = allocVector(STRSXP,4+3*dostopprob));    // 7 instead of 4 because of stopprob
  
  for(i = 0; i < (4+3*dostopprob); i++)                          // 7 instead of 4 because of stopprob
    SET_STRING_ELT(list_names,i,mkChar(names[i]));
  setAttrib(list, R_NamesSymbol, list_names);
  if (dostopprob){
    R_Free(stopn);
    R_Free(stopSn);
    R_Free(stopProb);
  }
  R_Free(p0); /// should be done in a better way - see e.g. pwilcox
  UNPROTECT(11+3*dostopprob);
  return list;
}

