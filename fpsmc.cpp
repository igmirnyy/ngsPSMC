#include <vector>
#include <cassert>
#include <cmath>
#include <pthread.h>
#include "psmcreader.h"
#include "main_psmc.h"
#include "hmm_psmc.h"
#include "bfgs.h"
#include <errno.h>

extern int nThreads;

int nChr = 0;

typedef struct{
  double *tk;
  int tk_l;
  double *epsize;
  double theta;
  double rho;
}shared_forhmm;


typedef struct{
  double **nP;
  double **PP;
  double *tk;
  int tk_l;
  double pix;
  int numWind;
  double theta;
  double rho;
  //  double *epsize;
}oPars;

int remap_l;
double *remap;


shared_forhmm shmm;


fastPSMC **objs = NULL;
oPars *ops = NULL;

/*
  objective function. Function to be optimized, for each chromo
*/

double qFunction_inner(double *tk,int tk_l,const double *epsize,double rho,double pix,int numWind,double **nP,double **PP);


double qFunction(const double *params ,const void *d){
  oPars *data = (oPars*) d;

  return qFunction_inner(data->tk,data->tk_l,params,data->rho,data->pix,data->numWind,data->nP,data->PP);

}
static int ncals=0;
double qFunction_wrapper(const double *pars,const void *){
  ncals++;
  fprintf(stderr,"\t-> calling objective function: remap_l:%d\n\n",remap_l);
  double pars2[ops[0].tk_l];
  int at=0;
  for(int i=0;i<remap_l;i++)
    for(int j=0;j<remap[i];j++)
      pars2[at++] = pars[i]; 
  double ret =0;
  for(int i=0;i<nChr;i++)
    ret += qFunction(pars2,&ops[i]);
  
  fprintf(stderr,"qfun:%f\n",ret);

  return -ret;
}


void make_remapper(psmc_par *pp){
  fprintf(stderr,"makeadsfadsfasfasfadf\n");
  int at=0;
  remap = new double[pp->n_free];
  remap_l=pp->n_free;
  for(int i=0;i<pp->n_free;i++){
    int howMany=0;
    while(pp->par_map[at+howMany]==pp->par_map[at])
      howMany++;
    at+=howMany;
    remap[i] = howMany;
  }

}


void runoptim3(double *tk,int tk_l,double *epsize,double theta,double rho,psmc_par *pp){
#if 0
  fprintf(stderr,"pp->n:%d\n",pp->n);
  fprintf(stderr,"pp->n_free:%d\n",pp->n_free);
  for(int i=0;i<pp->n+1;i++)
    fprintf(stderr,"pars_map[%d]\t%d\n",i,pp->par_map[i]);
#endif
  make_remapper(pp);
  int ndim = pp->n_free;

  double pars[ndim];
  for(int i=0;i<ndim;i++)
    pars[i] =0.5;//1.0;// 5.0*drand48()+1e-6;

  //set bounds
  int nbd[ndim];
  double lbd[ndim];
  double ubd[ndim];
  for(int i=0;i<ndim;i++){
    nbd[i]=1;
    lbd[i]=0.000001;
    ubd[i]=PSMC_T_INF;
  }

  for(int i=0;i<nChr;i++){
    ops[i].nP = objs[i]->nP;
    ops[i].PP = objs[i]->PP;
    ops[i].tk = tk;
    ops[i].tk_l = tk_l;
    ops[i].pix = objs[i]->pix;
    ops[i].numWind=objs[i]->windows.size();
    ops[i].rho= rho;
    ops[i].theta = theta;
  }
  double max_llh = findmax_bfgs(ndim,pars,NULL,qFunction_wrapper,NULL,lbd,ubd,nbd,1);
  fprintf(stderr,"\t-> optim done: after ncalls:%d\n",ncals);
  for(int i=0;i<ndim;i++)
    fprintf(stderr,"newpars: %f\n", pars[i]);
}

void printarray(FILE *fp,double *ary,int l);
void printmatrix(FILE *fp,double **mat,int x,int y);



void printmatrixf(char *fname,double **m,int x,int y){
  //  return ;
  FILE *fp = NULL;
  if(!(fp=fopen(fname,"wb"))){
    fprintf(stderr,"\t-> Problem writing file: \'%s\'\n",fname);
    exit(0);
  }
  printmatrix(fp,m,x,y);
  fclose(fp);
}

void printarrayf(char *fname,double *m,int x){
  //  return;
  FILE *fp = NULL;
  if(!(fp=fopen(fname,"wb"))){
    fprintf(stderr,"\t-> Problem writing file: \'%s\'\n",fname);
    exit(0);
  }
  printarray(fp,m,x);
  fclose(fp);
}

void *run_a_hmm(void *ptr){
  size_t at =(size_t) ptr;
  //  fprintf(stderr,"at:%lu\n",at);
  //  sleep(drand48()*10);
  objs[at]->make_hmm(shmm.tk,shmm.tk_l,shmm.epsize,shmm.theta,shmm.rho);
  pthread_exit(NULL);
}


void main_analysis_make_hmm(double *tk,int tk_l,double *epsize,double theta,double rho){

  fprintf(stderr,"\t-> [%s:%s:%d] nthreads:%d tk_l:%d theta:%f rho:%f\n",__FILE__,__FUNCTION__,__LINE__,nThreads,tk_l,theta,rho);
  shmm.tk=tk;
  shmm.tk_l=tk_l;
  shmm.theta=theta;
  shmm.rho=rho;
  shmm.epsize=epsize;

  pthread_t thread[nThreads];
  //  double qval =0;
  if(nThreads==1)
    for(int i=0;i<nChr;i++){
      objs[i]->make_hmm(shmm.tk,shmm.tk_l,shmm.epsize,shmm.theta,shmm.rho);
  }else {
    int at=0;
    while(at<nChr){
      int thisround = std::min(nChr-at,nThreads);
      for(int t=0;t<thisround;t++){
	size_t index = at+t;
	if(pthread_create( &thread[t], NULL, run_a_hmm, (void*) index)){
	  fprintf(stderr,"[%s] Problem spawning thread\n%s\n",__FUNCTION__,strerror(errno));
	  exit(0);
	}
      }
      for(int t=0;t<thisround;t++){
	if(pthread_join( thread[t], NULL)){
	  fprintf(stderr,"[%s] Problem joining thread\n%s\n",__FUNCTION__,strerror(errno));
	  exit(0);
	}
      }
      at+=thisround;
    }
  }

#if 1
  double fwllh,bwllh,qval;
  fwllh=bwllh=qval=0;
  for(int i=0;i<nChr;i++){
    //    fprintf(stderr,"\t-> hmm.fwllh for chr:%d\n",i);
    fwllh += objs[i]->fwllh;
    bwllh += objs[i]->bwllh;
    qval += objs[i]->qval;
  }
  fprintf(stderr,"\t[total llh]  fwllh:%f\n\t[total llh]  bwllh:%f\n\t[total qval] qval:%f\n",fwllh,bwllh,qval);
#endif
}


void main_analysis_optim(double *tk,int tk_l,double *epsize,double theta,double rho){

  shmm.tk=tk;
  shmm.tk_l=tk_l;
  shmm.theta=theta;
  shmm.rho=rho;
  shmm.epsize=epsize;

  pthread_t thread[nThreads];
  if(nThreads==1)
    for(int i=0;i<nChr;i++)
      objs[i]->make_hmm(shmm.tk,shmm.tk_l,shmm.epsize,shmm.theta,shmm.rho);
  else {
    int at=0;
    while(at<nChr){
      int thisround = std::min(nChr-at,nThreads);
      for(int t=0;t<thisround;t++){
	size_t index = at+t;
	if(pthread_create( &thread[t], NULL, run_a_hmm, (void*) index)){
	  fprintf(stderr,"[%s] Problem spawning thread\n%s\n",__FUNCTION__,strerror(errno));
	  exit(0);
	}
      }
      for(int t=0;t<thisround;t++){
	if(pthread_join( thread[t], NULL)){
	  fprintf(stderr,"[%s] Problem joining thread\n%s\n",__FUNCTION__,strerror(errno));
	  exit(0);
	}
      }
      at+=thisround;
    }
  }
  double fwllh,bwllh,qval;
  fwllh=bwllh=qval=0;
  for(int i=0;i<nChr;i++){
    //    fprintf(stderr,"\t-> hmm.fwllh for chr:%d\n",i);
    fwllh += objs[i]->fwllh;
    bwllh += objs[i]->bwllh;
    qval += objs[i]->qval;
  }
  fprintf(stderr,"\t[total llh]  fwllh:%f\n\t[total llh]  bwllh:%f\n\t[total qval] qval:%f\n",fwllh,bwllh,qval);

}

double calc_all_qvals(double *pars){
  double sum_qval =0;

  return sum_qval;
}


void superoptimizer2000(double *tk,int tk_l,double *epsize,double theta,double rho,psmc_par *pp){
  fprintf(stderr,"pp->n:%d\n",pp->n);
  fprintf(stderr,"pp->n_free:%d\n",pp->n_free);
  for(int i=0;i<pp->n+1;i++)
    fprintf(stderr,"pars_map[%d]\t%d\n",i,pp->par_map[i]);
  
}


void main_analysis(double *tk,int tk_l,double *epsize,double theta,double rho,psmc_par *pp){

  //first make_hmm for all chrs;
  // theta=0.0001;
  //  rho=0.1;
  main_analysis_make_hmm(tk,tk_l,epsize,theta,rho);
  runoptim3(tk,tk_l,epsize,theta,rho,pp);
  if(0){
    double tmptk[tk_l];
    double tmpepsize[tk_l];
    double tmptheta;
    double tmprho;
    void setpars2(char*fname,int which,double *tmptk,double *tmpepsize,double &tmptheta,double &tmprho);
    setpars2((char*)"test8/ms.lh3.fa.gz.psmc",4,tmptk,tmpepsize,tmptheta,tmprho);
    fprintf(stderr,"\\\\\\\\\\NEW PARS\n");
    tmprho=rho;
    tmptheta=theta;
    fprintf(stderr,"rho:\t%f\t->\t%f\n",rho,tmprho);
    fprintf(stderr,"theta:\t%f\t->\t%f\n",theta,tmptheta);
    for(int i=0;i<tk_l;i++)
      fprintf(stderr,"epsize[%d]\t%f\t%f\n",i,epsize[i],tmpepsize[i]);
    //  runoptim2(tk,tk_l,epsize,rho);
    double qvals[objs[0]->tot_index];
    double tot_qval =0;
    for(int i=0;i<objs[0]->tot_index;i++){
      fprintf(stderr,"\t->Running qval extra on: %d ",i);
      qvals[i] = qFunction_inner(tmptk,tk_l,tmpepsize,tmprho,objs[i]->pix,objs[i]->windows.size(),objs[i]->nP,objs[i]->PP);
      tot_qval += qvals[i];
      fprintf(stderr," qval[%d]: %f tot_qval:%f\n",i,qvals[i],tot_qval);
    }

  }
}

int psmc_wrapper(args *pars,int block) {
  fprintf(stderr,"\t-> we are in file: %s function: %s line:%d\n",__FILE__,__FUNCTION__,__LINE__);
  psmc_par *p=pars->par;
#if 1 //print pars
  fprintf(stderr,"\t-> par->n:%d\tpar->n_free:%d\tpar_map:%p\tpar->pattern:%s\tpar->times:%p\tpar->params:%p\n",p->n,p->n_free,p->par_map,p->pattern,p->times,p->params);
  for(int i=0;1&&i<pars->par->n+1;i++)
    fprintf(stderr,"[psmc_wrapper]:%i)\t%f\t%f\n",i,pars->par->times[i],pars->par->params[i]);
  //  exit(0);
#endif
  int tk_l = pars->par->n+1;
  double *tk = new double [tk_l];
  double *epsize = new double [tk_l];
  setEPSize(epsize,tk_l,p->params);
  //(nelems,array,max_t,alpha,array with values from file, can be NULL)
  setTk(tk_l,tk,15,0.01,p->times);//<- last position will be infinity
  //  fprintf(stderr,"[%s] tk=(%f,%f)\n",__FUNCTION__,tk[0],tk[1]);//exit(0);
#if 0
  for(int i=0;i<tk_l;i++)
    fprintf(stderr,"(tk,epsize)[%d]:(%e,%e)\n",i,tk[i],epsize[i]);
#endif
  
  //initialize all hmm (one for each chr), for now just a single
  int nobs = pars->chooseChr?1:pars->perc->mm.size();
  fprintf(stderr,"\t-> nobs/nchr: %d\n",nobs);
  objs = new fastPSMC*[nobs];
  ops = new oPars[nobs];
  for (myMap::const_iterator it = pars->perc->mm.begin() ;it!=pars->perc->mm.end();it++) {
    myMap::const_iterator it2;
    if(pars->chooseChr!=NULL)
      it2 = iter_init(pars->perc,pars->chooseChr,pars->start,pars->stop);
    else
      it2 = iter_init(pars->perc,it->first,pars->start,pars->stop);
    //    fprintf(stderr,"\t-> Parsing chr:%s \n",it2->first);
    fastPSMC *obj=objs[nChr++]=new fastPSMC;
    //    fprintf(stderr,"gls1:%f %f %f %f\n",pars->perc->gls[0],pars->perc->gls[1],pars->perc->gls[2],pars->perc->gls[3]);
    obj->setWindows(pars->perc->gls,pars->perc->pos,pars->perc->last,pars->block);
    //fprintf(stderr,"gls2:%f %f %f %f\n",pars->perc->gls[0],pars->perc->gls[1],pars->perc->gls[2],pars->perc->gls[3]);
    //  obj->printWindows(stdout);exit(0);
    obj->allocate(tk_l);
    if(pars->chooseChr!=NULL)
      break;
  }
  main_analysis(tk,tk_l,epsize,pars->par->TR[0],pars->par->TR[1],pars->par);

  
  
#if 0
  printarrayf("tk",tk,tk_l);
  printmatrixf("fw",objs[0]->fw,tk_l,objs[0]->windows.size()+1);
  printmatrixf("bw",objs[0]->bw,tk_l,objs[0]->windows.size()+1);
  printmatrixf("emis",objs[0]->emis,tk_l,objs[0]->windows.size()+1);
  printmatrixf("P",obj.P,7,tk_l);
#endif
  //printmatrixf("pp",objs[0]->pp,tk_l,objs[0]->windows.size()+1);
  // 
  /*
    printarrayf("stationary",obj.stationary,tk_l);
    
    printarrayf("epsize",epsize,tk_l);

    printmatrixf("emis",obj.emis,tk_l,obj.
windows.size()+1);

    printarrayf("r1",obj.R1,tk_l);
    printarrayf("r2",obj.R2,tk_l);
  */
  
  //calculate window end points
  
  //  obj.printWindows(stdout);
  //allocate internal structures needed
 
  //make an hmm

  // objs[0].make_hmm(tk,tk_l,epsize);
  /*
    printarrayf("stationary",obj.stationary,tk_l);
    printarrayf("tk",tk,tk_l);
    printarrayf("epsize",epsize,tk_l);
    printmatrixf("P",obj.P,7,tk_l);
    printmatrixf("emis",obj.emis,tk_l,obj.windows.size()+1);
    printmatrixf("fw",obj.fw,tk_l,obj.windows.size()+1);
    printarrayf("r1",obj.R1,tk_l);
    printarrayf("r2",obj.R2,tk_l);
  */
  return 1;
}
