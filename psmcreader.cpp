#include <sys/stat.h>
#include <ctime>
#include <htslib/bgzf.h>
#include <htslib/faidx.h>
#include <htslib/hts.h>
#include <zlib.h>
#include <cmath>
#include "header.h"
#include "psmcreader.h"
#include <htslib/vcf.h>
#include <map>
#include <vector>

void destroy(myMap& mm) {
    for (myMap::iterator it = mm.begin();it != mm.end();++it)
        free(it->first);
    mm.clear();
}


double pl_to_gl(double pl) {
    return pow(10,  -pl / 10);
}



void perpsmc_destroy(perpsmc* pp) {
    if (pp->pf)
        perFasta_destroy(pp->pf);
    if (pp->pb)
        perBcf_destroy(pp->pb);
    free(pp->bgzf_gls);
    free(pp->bgzf_pos);
    destroy(pp->mm);

    free(pp->fname);

    delete pp;
}

void writepsmc_header(FILE* fp, perpsmc* pp, int onlysubset) {
    fprintf(fp, "\t\tInformation from index file: nSites(total):%lu nChr:%lu\n", pp->nSites, pp->mm.size());

    int i = 0;
    for (myMap::const_iterator it = pp->mm.begin();it != pp->mm.end();++it) {
        datum d = it->second;
        fprintf(fp, "\t\t%d\t%s\t%zu\t%ld\t%ld\n", i++, it->first, d.nSites, (long int)d.pos, (long int)d.saf);
        if (onlysubset && i > 9) {
            fprintf(stderr, "\t\t Breaking printing of header\n");
            break;
        }
    }

}


int psmcversion(const char* fname) {//0 - fasta, 1 - psmc, 2 - vcf/bcf
    gzFile gz = Z_NULL;
    gz = gzopen(fname, "r");
    if (gz == Z_NULL) {
        fprintf(stderr, "\t-> Problem opening file: \'%s\'", fname);
        exit(0);
    }
    char buf[8];
    gzread(gz, buf, 8 * sizeof(char));
    //  fprintf(stderr,"\t-> Magic nr is: \'%s\'\n",buf);
    gzclose(gz);

    if (0 == strcmp(buf, "psmcv1"))
        return 1;
    else {
        htsFile* vcf_file = bcf_open(fname, "r");
        bcf_hdr_t* hdr = bcf_hdr_read(vcf_file);
        if (hdr != NULL) {
            bcf_hdr_destroy(hdr);
            return 2;
        }
        bcf_hdr_destroy(hdr);
        return 0;
    }
}


void init_for_fasta(const char* fname, perpsmc* ret, int nChr, int& at){
    ret->pf = perFasta_init(fname);
    for (int i = 0; i < faidx_nseq(ret->pf->fai); i++) {
        if (nChr != -1 && at++ >= nChr)
            break;
        char* chr = strdup(faidx_iseq(ret->pf->fai, i));
        // fprintf(stderr,"\t-> [%s] %d) chr: %s\n",__FUNCTION__,i,chr);
        datum d;
        d.nSites = faidx_seq_len(ret->pf->fai, chr);
        ret->nSites += d.nSites;
        d.pos = d.saf = 0;
        myMap::iterator it = ret->mm.find(chr);
        if (it == ret->mm.end())
            ret->mm[chr] = d;
        else {
            fprintf(stderr,
                "Problem with chr: %s, key already exists, try to reindex fasta file\n",
                chr);
            exit(1);
        }
    }
}
void init_for_bcf(const char* fname, perpsmc* ret, int nChr, int& at){
ret->pb = perBcf_init(fname);
int n_seqs, i, rec_id;
uint64_t mapped, unmapped;
char* chr;
const char ** seqnames = bcf_index_seqnames(ret->pb->idx, ret->pb->hdr, &n_seqs);
for (i =0; i<n_seqs; i++){
    chr = strdup(seqnames[i]);
    datum d;
    rec_id = bcf_hdr_name2id(ret->pb->hdr, seqnames[i]);
    if (rec_id >= 0){
        hts_idx_get_stat(ret->pb->idx, rec_id, &mapped, &unmapped);
        d.nSites += mapped;
    }
    else {
        fprintf(stderr,
            "Problem with chr: %s, try to reindex bcf file\n",
            chr);
        exit(1);
    }
    ret->nSites += d.nSites;
    d.pos = d.saf = 0;
    myMap::iterator it = ret->mm.find(chr);
    if (it == ret->mm.end())
        ret->mm[chr] = d;
    else {
        fprintf(stderr,
            "Problem with chr: %s, key already exists, try to reindex bcf file\n",
            chr);
        exit(0);
    }
 }
}
void init_for_psmc(const char * index_fname, FILE* fp, perpsmc* ret, int nChr, int& at){
size_t clen;
while (fread(&clen, sizeof(size_t), 1, fp)) {
    if (nChr != -1 && at++ >= nChr)
        break;
    char* chr = (char*)malloc(clen + 1);
    assert(clen == fread(chr, 1, clen, fp));
    chr[clen] = '\0';

    datum d;
    if (1 != fread(&d.nSites, sizeof(size_t), 1, fp)) {
        fprintf(stderr, "[%s.%s():%d] Problem reading data: %s \n", __FILE__, __FUNCTION__, __LINE__, index_fname);
        exit(0);
    }
    ret->nSites += d.nSites;
    if (1 != fread(&d.pos, sizeof(int64_t), 1, fp)) {
        fprintf(stderr, "[%s->%s():%d] Problem reading data: %s \n", __FILE__, __FUNCTION__, __LINE__, index_fname);
        exit(0);
    }
    if (1 != fread(&d.saf, sizeof(int64_t), 1, fp)) {
        fprintf(stderr, "[%s->%s():%d] Problem reading data: %s \n", __FILE__, __FUNCTION__, __LINE__, index_fname);
        exit(0);
    }

    myMap::iterator it = ret->mm.find(chr);
    if (it == ret->mm.end())
        ret->mm[chr] = d;
    else {
        fprintf(stderr,
            "Problem with chr: %s, key already exists, psmc file needs to be sorted. (sort your -rf that you used for input)\n",
            chr);
        exit(1);
    }
}
fclose(fp);
char* base_fname = (char*)calloc(strlen(index_fname) + 100, 1);//that should do it
base_fname = strncpy(base_fname, index_fname, strlen(index_fname) - 3);
//  fprintf(stderr,"tmp:%s\n",tmp);

char* gl_fname = (char*)calloc(strlen(base_fname) + 100, 1);//that should do it
snprintf(gl_fname, strlen(index_fname) + 100, "%sgz", base_fname);
fprintf(stderr, "\t-> Assuming .psmc.gz file: \'%s\'\n", gl_fname);
ret->bgzf_gls = strdup(gl_fname);
BGZF* tmpfp = NULL;
tmpfp = bgzf_open(ret->bgzf_gls, "r");
if (tmpfp)
    my_bgzf_seek(tmpfp, 8, SEEK_SET);
if (tmpfp && ret->version != psmcversion(gl_fname)) {
    fprintf(stderr, "\t-> Problem with mismatch of version of %s vs %s %d vs %d\n", index_fname, gl_fname, ret->version,
        psmcversion(gl_fname));
    exit(0);
}
bgzf_close(tmpfp);
tmpfp = NULL;

char* pos_fname = gl_fname;
snprintf(pos_fname, strlen(index_fname) + 100, "%spos.gz", base_fname);
fprintf(stderr, "\t-> Assuming .psmc.pos.gz: \'%s\'\n", pos_fname);
ret->bgzf_pos = strdup(pos_fname);
tmpfp = bgzf_open(ret->bgzf_pos, "r");
if (tmpfp)
    my_bgzf_seek(tmpfp, 8, SEEK_SET);
if (tmpfp && ret->version != psmcversion(pos_fname)) {
    fprintf(stderr, "Problem with mismatch of version of %s vs %s\n", index_fname, pos_fname);
    exit(0);
}
bgzf_close(tmpfp);

free(base_fname);
free(gl_fname);
}



perpsmc* perpsmc_init(char* fname, int nChr) {
    assert(fname);
    perpsmc* ret = new perpsmc;
    ret->fname = strdup(fname);
    ret->bgzf_pos = ret->bgzf_gls = NULL;
    ret->pf = NULL;


    if (!fexists(fname)) {
        fprintf(stderr, "\t-> Problem opening file: \'%s\'\n", fname);
        exit(0);
    }
    FILE* fp = NULL;
    fp = fopen(fname, "r");
    if (fp == NULL) {
        fprintf(stderr, "\t-> Problem opening file:%s\n", fname);
        exit(0);
    }
    char buf[8];
    assert(fread(buf, 1, 8, fp) == 8);
    ret->version = psmcversion(fname);
    ret->nSites = 0;
    fprintf(stderr, "\t-> Version of fname: \'%s\' is:%d\n", fname, ret->version);
    int at = 0;//incrementer for breaking out of filereading if -nChr has been supplied

    //loop for fasta
    switch (ret->version)
    {
    case 0:
        fprintf(stderr,
        "\t-> Looks like you are trying to use a version of PSMC that does not exists, assuming its a fastafile\n");
        fclose(fp);
        init_for_fasta(fname, ret, nChr, at);
        break;
    case 1:
        init_for_psmc(fname, fp, ret, nChr, at);
        break;
    case 2:
        fprintf(stderr,
        "\t-> Reading a bcf file\n");
        fclose(fp);
        init_for_bcf(fname, ret, nChr, at);
    default:
        break;
    }
    return ret;
}


BGZF* bgzf_open_seek(char* fname, int64_t offs) {
    BGZF* ret = NULL;
    ret = bgzf_open(fname, "r");
    my_bgzf_seek(ret, offs, SEEK_SET);
    return ret;
}


void readstuff_from_psmc(perpsmc* pp, myMap::iterator it, rawdata ret){
    double* tmpgls = new double[2 * it->second.nSites];
    BGZF* bgzf_gls = bgzf_open_seek(pp->bgzf_gls, it->second.saf);
    BGZF* bgzf_pos = bgzf_open_seek(pp->bgzf_pos, it->second.pos);

    my_bgzf_read(bgzf_pos, ret.pos, sizeof(int) * it->second.nSites);
    my_bgzf_read(bgzf_gls, tmpgls, 2 * sizeof(double) * it->second.nSites);
    for (int i = 0;i < it->second.nSites;i++) {
        ret.gls[i] = log(0);
        if (tmpgls[2 * i] != tmpgls[2 * i + 1]) {
            double mmax = std::max(tmpgls[2 * i], tmpgls[2 * i + 1]);
            double val = std::min(tmpgls[2 * i], tmpgls[2 * i + 1]) - mmax;

            ret.gls[i] = val;
            if (tmpgls[2 * i] < tmpgls[2 * i + 1])
                ret.gls[i] = -ret.gls[i];

            //code here should be implemented for using phredstyle gls //if(sizeof(mygltype))

        }
    }
    delete[] tmpgls;
    bgzf_close(bgzf_gls);
    bgzf_close(bgzf_pos);
}

void readstuff_from_fasta(perpsmc* pp, myMap::iterator it, rawdata ret, int blockSize){
    int asdf = it->second.nSites;
    char* tmp = faidx_fetch_seq(pp->pf->fai, it->first, 0, 0x7fffffff, &asdf);
    for (int i = 0;i < it->second.nSites;i++) {
        ret.pos[i] = i * blockSize;
        //important relates to problems with divide by zero in compuation of  backward probablity
        //K=het
        if (tmp[i] == 'K')
            ret.gls[i] = 500.0;// 0;//het
        else
            ret.gls[i] = -500.0;//;//hom

        //ok let me explain. negative means homozygotic and postive means heteroeo. The otherone is always 0.

        //       fprintf(stderr,"%c\n",tmp[i]);
    }
    free(tmp);
}

long readstuff_from_bcf(perpsmc* pp, myMap::iterator it, rawdata ret){
    long i = 0;
    int* ploidy;
    int pl_arr_len;
    hts_itr_t* iter = bcf_itr_querys(pp->pb->idx, pp->pb->hdr, it->first);
    bcf1_t* rec = bcf_init();
    double homo_pl, hetero_pl;
    while  (bcf_itr_next(pp->pb->bcf_file, iter, rec) == 0) {
        bcf_unpack(rec, BCF_UN_STR);
        if (bcf_get_info_flag(pp->pb->hdr, rec, "INDEL", NULL, NULL) == 1 || rec->d.als[0] == 'N') continue;
        ret.pos[i] = rec->pos + 1;
        bcf_get_format_int32(pp->pb->hdr, rec, "PL", &ploidy, &pl_arr_len);
        //Counting PL scores
        switch (rec->n_allele) {
        case 2:
            homo_pl = log(pl_to_gl(ploidy[0]) + pl_to_gl(ploidy[2]));
            hetero_pl = log(pl_to_gl(ploidy[1]));
            break;
        case 3:
            homo_pl =log(pl_to_gl(ploidy[0]) + pl_to_gl(ploidy[3]) + pl_to_gl(ploidy[5]));
            hetero_pl = log(pl_to_gl(ploidy[1]) + pl_to_gl(ploidy[2]) + pl_to_gl(ploidy[4]));
            break;
        case 4:
            homo_pl = log(pl_to_gl(ploidy[0]) + pl_to_gl(ploidy[4]) + pl_to_gl(ploidy[7]) + pl_to_gl(ploidy[9]));
            hetero_pl = log(pl_to_gl(ploidy[1]) + pl_to_gl(ploidy[2]) + pl_to_gl(ploidy[3]) + pl_to_gl(ploidy[5]) + pl_to_gl(ploidy[6]) + pl_to_gl(ploidy[8]));
        default:
            break;
        }
        ret.gls[i] = hetero_pl - homo_pl;
        if (i % 4000 == 0){
            fprintf(stderr, "%ld pos: %d gl %lf\n", i, ret.pos[i], ret.gls[i]);
        }
        i++;
        if (i > 100000){
            break;
        }
    }
    fprintf(stderr, "Read chromosome %s expected %ld got %ld", it->first, it->second.nSites , i);
    ret.len = i;
    exit(0);
}

//this functions returns the emissions
rawdata readstuff(perpsmc* pp, char* chr, int blockSize, int start, int stop) {
    rawdata ret;
    assert(chr != NULL);

    myMap::iterator it = pp->mm.find(chr);
    if (it == pp->mm.end()) {
        fprintf(stderr, "\t-> [%s] Problem finding chr: \'%s\'\n", __FUNCTION__, chr);
        exit(0);
    }

    ret.pos = new int[it->second.nSites];
    ret.len = it->second.nSites;
    ret.gls = new double[it->second.nSites];

    switch (pp->version)
    {
    case 0:
        readstuff_from_fasta(pp, it, ret, blockSize);
        break;
    case 1:
        readstuff_from_psmc(pp, it, ret);
        break;
    case 2:
        ret.len = readstuff_from_bcf(pp, it, ret);
    default:
        break;
    } 
    ret.firstp = 0;
    ret.lastp = it->second.nSites;

#if 1 //the code below should be readded if we ever want to run on specific specified regions
    ret.firstp = 0;
    if (start != -1)
        while (ret.firstp < it->second.nSites && ret.pos[ret.firstp] < start)
            ret.firstp++;

    ret.lastp = it->second.nSites;
    if (stop != -1 && stop <= ret.pos[ret.lastp - 1]) {
        ret.lastp = ret.firstp;
        while (ret.pos[ret.lastp] < stop)
            ret.lastp++;
    }
#endif
    return ret;
}

