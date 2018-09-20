#include <sys/stat.h>
#include <ctime>
#include <htslib/bgzf.h>
#include <htslib/faidx.h>
#include <zlib.h>
#include <cmath>
#include "header.h"
#include "psmcreader.h"



void destroy(myMap &mm){
  for(myMap::iterator it=mm.begin();it!=mm.end();++it)
    free(it->first);
  mm.clear();
}

void perpsmc_destroy(perpsmc *pp){
  if(pp->pf)
    perFasta_destroy(pp->pf);
  bgzf_close(pp->bgzf_gls);
  bgzf_close(pp->bgzf_pos);
  destroy(pp->mm);
  
  if(pp->pos)
    delete [] pp->pos;

  if(pp->tmpgls)
    delete [] pp->tmpgls;

  free(pp->fname);
  
  delete pp;
}

void writepsmc_header(FILE *fp,perpsmc *pp,int onlysubset){
  fprintf(fp,"\t\tInformation from index file: nSites(total):%lu nChr:%lu\n",pp->nSites,pp->mm.size());
  
  int i=0;
  for(myMap::const_iterator it=pp->mm.begin();it!=pp->mm.end();++it){
    datum d = it->second;
    fprintf(fp,"\t\t%d\t%s\t%zu\t%ld\t%ld\n",i++,it->first,d.nSites,(long int)d.pos,(long int)d.saf);
    if(onlysubset && i>9){
      fprintf(stderr,"\t\t Breaking printing of header\n");
      break;
    }
  }

}


int psmcversion(const char *fname){
  gzFile gz=Z_NULL;
  gz = gzopen(fname,"r");
  if(gz==Z_NULL){
    fprintf(stderr,"\t-> Problem opening file: \'%s\'",fname);
    exit(0);
  }
  char buf[8];
  gzread(gz,buf,8*sizeof(char));
  //  fprintf(stderr,"\t-> Magic nr is: \'%s\'\n",buf);
  gzclose(gz);

  if(0==strcmp(buf,"psmcv1"))
    return 1;
  else 
    return 0;
}



perpsmc * perpsmc_init(char *fname,int nChr){
  assert(fname);
  perpsmc *ret = new perpsmc ;
  ret->fname = strdup(fname);
  ret->gls =NULL;
  ret->pos = NULL;
  ret->bgzf_pos=ret->bgzf_gls=NULL;
  ret->pos = NULL;
  ret->pos_l = 0;
  ret->pf = NULL;
  ret->tmpgls = NULL;
  ret->tmpgls_l = 0;
  size_t clen;
  if(!fexists(fname)){
    fprintf(stderr,"\t-> Problem opening file: \'%s\'\n",fname);
    exit(0);
  }
  FILE *fp = NULL;
  fp=fopen(fname,"r");
  if(fp==NULL){
    fprintf(stderr,"\t-> Problem opening file:%s\n",fname);
    exit(0);
  }
  char buf[8];
  assert(fread(buf,1,8,fp)==8);
  ret->version = psmcversion(fname);
  ret->nSites =0;
  fprintf(stderr,"\t-> Version of fname: \'%s\' is:%d\n",fname,ret->version);
  int at=0;//incrementer for breaking out of filereading if -nChr has been supplied

  //loop for fasta
  if(ret->version!=1){
    fprintf(stderr,"\t-> Looks like you are trying to use a version of PSMC that does not exists, assuming its a fastafile\n");
    fclose(fp);
    fp=NULL;
    ret->pf = perFasta_init(fname);
    for(int i=0;i<faidx_nseq(ret->pf->fai);i++){
      if(nChr!=-1&&at++>=nChr)
	break;
      char *chr = strdup(faidx_iseq(ret->pf->fai,i));
      // fprintf(stderr,"\t-> [%s] %d) chr: %s\n",__FUNCTION__,i,chr);
      datum d;
      d.nSites = faidx_seq_len(ret->pf->fai,chr);
      ret->nSites += d.nSites;
      d.pos=d.saf=0;
      myMap::iterator it = ret->mm.find(chr);
      if(it==ret->mm.end())
	ret->mm[chr] =d ;
      else{
	fprintf(stderr,"Problem with chr: %s, key already exists, psmc file needs to be sorted. (sort your -rf that you used for input)\n",chr);
	exit(0);
      }
    }
    return ret;
  }

  //loop for gl
  while(fread(&clen,sizeof(size_t),1,fp)){
    if(nChr!=-1&&at++>=nChr)
      break;
    char *chr = (char*)malloc(clen+1);
    assert(clen==fread(chr,1,clen,fp));
    chr[clen] = '\0';
    
    datum d;
    if(1!=fread(&d.nSites,sizeof(size_t),1,fp)){
      fprintf(stderr,"[%s.%s():%d] Problem reading data: %s \n",__FILE__,__FUNCTION__,__LINE__,fname);
      exit(0);
    }
    ret->nSites += d.nSites;
    if(1!=fread(&d.pos,sizeof(int64_t),1,fp)){
      fprintf(stderr,"[%s->%s():%d] Problem reading data: %s \n",__FILE__,__FUNCTION__,__LINE__,fname);
      exit(0);
    }
    if(1!=fread(&d.saf,sizeof(int64_t),1,fp)){
      fprintf(stderr,"[%s->%s():%d] Problem reading data: %s \n",__FILE__,__FUNCTION__,__LINE__,fname);
      exit(0);
    }
  
    myMap::iterator it = ret->mm.find(chr);
    if(it==ret->mm.end())
      ret->mm[chr] =d ;
    else{
      fprintf(stderr,"Problem with chr: %s, key already exists, psmc file needs to be sorted. (sort your -rf that you used for input)\n",chr);
      exit(0);
    }
  }
  fclose(fp);
  char *tmp =(char*)calloc(strlen(fname)+100,1);//that should do it
  tmp=strncpy(tmp,fname,strlen(fname)-3);
  //  fprintf(stderr,"tmp:%s\n",tmp);
  
  char *tmp2 = (char*)calloc(strlen(fname)+100,1);//that should do it
  snprintf(tmp2,strlen(fname)+100,"%sgz",tmp);
  fprintf(stderr,"\t-> Assuming .psmc.gz file: \'%s\'\n",tmp2);
  ret->bgzf_gls = bgzf_open(tmp2,"r");
  if(ret->bgzf_gls)
    my_bgzf_seek(ret->bgzf_gls,8,SEEK_SET);
  if(ret->bgzf_gls && ret->version!=psmcversion(tmp2)){
    fprintf(stderr,"\t-> Problem with mismatch of version of %s vs %s %d vs %d\n",fname,tmp2,ret->version,psmcversion(tmp2));
    exit(0);
  }

  snprintf(tmp2,strlen(fname)+100,"%spos.gz",tmp);
  fprintf(stderr,"\t-> Assuming .psmc.pos.gz: \'%s\'\n",tmp2);
  ret->bgzf_pos = bgzf_open(tmp2,"r");
  if(ret->pos)
    my_bgzf_seek(ret->bgzf_pos,8,SEEK_SET);
  if(ret->bgzf_pos&& ret->version!=psmcversion(tmp2)){
    fprintf(stderr,"Problem with mismatch of version of %s vs %s\n",fname,tmp2);
    exit(0);
  }
  //assert(ret->pos!=NULL&&ret->saf!=NULL);
  free(tmp);free(tmp2);
  
 return ret;
 }

 //chr start stop is given from commandine
myMap::iterator iter_init(perpsmc *pp,char *chr,int start,int stop,int blockSize){
   assert(chr!=NULL);
   assert(pp->gls==NULL);
   myMap::iterator it = pp->mm.find(chr);
   if(it==pp->mm.end()){
     fprintf(stderr,"\t-> [%s] Problem finding chr: \'%s\'\n",__FUNCTION__,chr);
     exit(0);
     return it;
   }
   //   fprintf(stderr,"%s->%lu,%lu,%lu\n",it->first,it->second.nSites,it->second.saf,it->second.pos);
   if(pp->version==1){
     my_bgzf_seek(pp->bgzf_gls,it->second.saf,SEEK_SET);
     my_bgzf_seek(pp->bgzf_pos,it->second.pos,SEEK_SET);
   }
   //fprintf(stderr,"pp->gls:%p\n",pp->gls);
   if(pp->pos_l<it->second.nSites){
     delete [] pp->pos;
     pp->pos = new int[it->second.nSites];
     pp->pos_l = it->second.nSites;
   }
   pp->gls = new mygltype[it->second.nSites];
   
   if(pp->tmpgls_l<2*it->second.nSites){
     //     fprintf(stderr,"reallocing all the time\n");
     delete [] pp->tmpgls;
     pp->tmpgls = new double[2*it->second.nSites];
     pp->tmpgls_l = 2*it->second.nSites;
   }
   
   if(pp->version==1) {
     my_bgzf_read(pp->bgzf_pos,pp->pos,sizeof(int)*it->second.nSites);
     my_bgzf_read(pp->bgzf_gls,pp->tmpgls,2*sizeof(double)*it->second.nSites);
     for(int i=0;i<it->second.nSites;i++){
       pp->gls[i] = log(0);
       //       fprintf(stderr,"precal res:%f\t0:%f\t1:%f\n",pp->gls[i],tmpgls[2*i],tmpgls[2*i+1]);
       if(pp->tmpgls[2*i]!=pp->tmpgls[2*i+1]){
	 double mmax = std::max(pp->tmpgls[2*i],pp->tmpgls[2*i+1]);
	 pp->tmpgls[2*i] -= mmax;
	 pp->tmpgls[2*i+1] -= mmax;
       }
       // fprintf(stderr,"post scal res:%f\t0:%f\t1:%f\n",pp->gls[i],tmpgls[2*i],tmpgls[2*i+1]);
       if(pp->tmpgls[2*i]>pp->tmpgls[2*i+1])
	 pp->gls[i]=pp->tmpgls[2*i+1];
       else
	 pp->gls[i]=-pp->tmpgls[2*i];
       // fprintf(stderr,"res:%f\t0:%f\t1:%f\n",pp->gls[i],tmpgls[2*i],tmpgls[2*i+1]);
     }
     //   fprintf(stderr," end: %f %f\n",pp->gls[0],pp->gls[1]);
     pp->first=0;
     if(start!=-1)
       while(pp->first<it->second.nSites&&pp->pos[pp->first]<start)
	 pp->first++;
     
     pp->last = it->second.nSites;
     if(stop!=-1&&stop<=pp->pos[pp->last-1]){
       pp->last=pp->first;
       while(pp->pos[pp->last]<stop) 
	 pp->last++;
     }
   }else{
     int asdf = it->second.nSites;
     char *tmp = faidx_fetch_seq(pp->pf->fai, it->first, 0, 0x7fffffff, &asdf);
     for(int i=0;i<it->second.nSites;i++){
       pp->pos[i] = i*blockSize;
       //important relates to problems with divide by zero in compuation of  backward probablity
       //K=het
       if(tmp[i]=='K')
	 pp->gls[i] = 500;// 0;//het 
       else
	 pp->gls[i] = -500;//;//hom

       //ok let me explain. negative means homozygotic and postive means heteroeo. The otherone is always 0.

       //       fprintf(stderr,"%c\n",tmp[i]);
     }
     free(tmp);
    
     pp->first=0;
     pp->last=it->second.nSites;
   }

#ifdef __SHOW_TIME__
   fprintf(stderr, "\t[TIME] cpu-time used =  %.2f sec for reading data\n", (float)(clock() - t) / CLOCKS_PER_SEC);
   fprintf(stderr, "\t[Time] walltime used =  %.2f sec for reading data\n", (float)(time(NULL) - t2));  
#endif
   return it;
 }
