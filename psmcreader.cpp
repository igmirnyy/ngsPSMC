#include <sys/stat.h>
#include <ctime>
#include <htslib/bgzf.h>
#include <htslib/faidx.h>
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

void perpsmc_destroy(perpsmc* pp) {
    if (pp->pf)
        perFasta_destroy(pp->pf);
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





perpsmc* perpsmc_init(char* fname, int nChr) {
    assert(fname);
    perpsmc* ret = new perpsmc;
    ret->fname = strdup(fname);
    ret->bgzf_pos = ret->bgzf_gls = NULL;
    ret->pf = NULL;

    size_t clen;
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
    if (ret->version != 2) {
        if (ret->version != 1) {
            fprintf(stderr,
                "\t-> Looks like you are trying to use a version of PSMC that does not exists, assuming its a fastafile\n");
            fclose(fp);
            fp = NULL;
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
                        "Problem with chr: %s, key already exists, psmc file needs to be sorted. (sort your -rf that you used for input)\n",
                        chr);
                    exit(0);
                }
            }
            return ret;
        }

        //loop for gl
        while (fread(&clen, sizeof(size_t), 1, fp)) {
            if (nChr != -1 && at++ >= nChr)
                break;
            char* chr = (char*)malloc(clen + 1);
            assert(clen == fread(chr, 1, clen, fp));
            chr[clen] = '\0';

            datum d;
            if (1 != fread(&d.nSites, sizeof(size_t), 1, fp)) {
                fprintf(stderr, "[%s.%s():%d] Problem reading data: %s \n", __FILE__, __FUNCTION__, __LINE__, fname);
                exit(0);
            }
            ret->nSites += d.nSites;
            if (1 != fread(&d.pos, sizeof(int64_t), 1, fp)) {
                fprintf(stderr, "[%s->%s():%d] Problem reading data: %s \n", __FILE__, __FUNCTION__, __LINE__, fname);
                exit(0);
            }
            if (1 != fread(&d.saf, sizeof(int64_t), 1, fp)) {
                fprintf(stderr, "[%s->%s():%d] Problem reading data: %s \n", __FILE__, __FUNCTION__, __LINE__, fname);
                exit(0);
            }

            myMap::iterator it = ret->mm.find(chr);
            if (it == ret->mm.end())
                ret->mm[chr] = d;
            else {
                fprintf(stderr,
                    "Problem with chr: %s, key already exists, psmc file needs to be sorted. (sort your -rf that you used for input)\n",
                    chr);
                exit(0);
            }
        }
        fclose(fp);
        char* tmp = (char*)calloc(strlen(fname) + 100, 1);//that should do it
        tmp = strncpy(tmp, fname, strlen(fname) - 3);
        //  fprintf(stderr,"tmp:%s\n",tmp);

        char* tmp2 = (char*)calloc(strlen(fname) + 100, 1);//that should do it
        snprintf(tmp2, strlen(fname) + 100, "%sgz", tmp);
        fprintf(stderr, "\t-> Assuming .psmc.gz file: \'%s\'\n", tmp2);
        ret->bgzf_gls = strdup(tmp2);
        BGZF* tmpfp = NULL;
        tmpfp = bgzf_open(ret->bgzf_gls, "r");
        if (tmpfp)
            my_bgzf_seek(tmpfp, 8, SEEK_SET);
        if (tmpfp && ret->version != psmcversion(tmp2)) {
            fprintf(stderr, "\t-> Problem with mismatch of version of %s vs %s %d vs %d\n", fname, tmp2, ret->version,
                psmcversion(tmp2));
            exit(0);
        }
        bgzf_close(tmpfp);
        tmpfp = NULL;

        snprintf(tmp2, strlen(fname) + 100, "%spos.gz", tmp);
        fprintf(stderr, "\t-> Assuming .psmc.pos.gz: \'%s\'\n", tmp2);
        ret->bgzf_pos = strdup(tmp2);
        tmpfp = bgzf_open(ret->bgzf_pos, "r");
        if (tmpfp)
            my_bgzf_seek(tmpfp, 8, SEEK_SET);
        if (tmpfp && ret->version != psmcversion(tmp2)) {
            fprintf(stderr, "Problem with mismatch of version of %s vs %s\n", fname, tmp2);
            exit(0);
        }
        bgzf_close(tmpfp);

        free(tmp);
        free(tmp2);
    }
    return ret;
}

BGZF* bgzf_open_seek(char* fname, int64_t offs) {
    BGZF* ret = NULL;
    ret = bgzf_open(fname, "r");
    my_bgzf_seek(ret, offs, SEEK_SET);
    return ret;
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
    double* tmpgls = new double[2 * it->second.nSites];

    if (pp->version == 1) {
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
    else {
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

double pl_to_gl(double pl) {
    return exp(log(10) * -pl / 10);
}

std::map<const char*, rawdata> get_vcf_data(perpsmc* pp, int start, int stop) {
    std::map<const char*, rawdata> vcf_data;
    std::map<const char*, std::vector<int > > positions;
    std::map<const char*, std::vector<double > > likelihoods;
    int* ploidy = NULL;
    int pl_arr_len = 0;
    htsFile* input_file = bcf_open(pp->fname, "r");
    bcf_hdr_t* header = bcf_hdr_read(input_file);
    bcf1_t* record = bcf_init();
    int i = 0;
    //Reading data from vcf file
    while (bcf_read(input_file, header, record) == 0) {
        bcf_unpack(record, BCF_UN_STR);
        //Skipping INDELS and N in REF
        if (bcf_get_info_flag(header, record, "INDEL", NULL, NULL) == 1 || record->d.als[0] == 'N') continue;
        i++;
        double homo_pl = 0;
        double hetero_pl = 0;
        double likelihood;
        bcf_get_format_int32(header, record, "PL", &ploidy, &pl_arr_len);
        //Counting PL scores
        switch (record->n_allele) {
        case 2:
            homo_pl = pl_to_gl(ploidy[0]) + pl_to_gl(ploidy[2]);
            hetero_pl = pl_to_gl(ploidy[1]);

            break;
        case 3:
            homo_pl = pl_to_gl(ploidy[0]) + pl_to_gl(ploidy[3]) + pl_to_gl(ploidy[5]);
            hetero_pl = pl_to_gl(ploidy[1]) + pl_to_gl(ploidy[2]) + pl_to_gl(ploidy[4]);
            break;
        case 4:
            homo_pl = pl_to_gl(ploidy[0]) + pl_to_gl(ploidy[4]) + pl_to_gl(ploidy[7]) + pl_to_gl(ploidy[9]);
            hetero_pl = pl_to_gl(ploidy[1]) + pl_to_gl(ploidy[2]) + pl_to_gl(ploidy[3]) + pl_to_gl(ploidy[5]) + pl_to_gl(ploidy[6]) + pl_to_gl(ploidy[8]);
        default:
            break;
        }

        homo_pl = homo_pl;
        hetero_pl = hetero_pl;
        likelihood =   hetero_pl > homo_pl? -5: 5;
        //Storing positions and likelihoods
        std::vector<int>& positions_vector = positions[bcf_hdr_id2name(header, record->rid)];
        std::vector<double>& likelihoods_vector = likelihoods[bcf_hdr_id2name(header, record->rid)];
        positions_vector.push_back(record->pos + 1);
        likelihoods_vector.push_back(likelihood);
    }
    bcf_hdr_destroy(header);
    bcf_destroy(record);
    bcf_close(input_file);
    rawdata output_rawdata;
    //Creating rawdata objects from vectors
    for (std::map<const char*, std::vector<int > >::iterator it = positions.begin();it != positions.end();it++) {
        output_rawdata = rawdata();
        output_rawdata.pos = new int[positions[it->first].size()];
        memcpy(output_rawdata.pos, positions[it->first].data(), positions[it->first].size() * sizeof(int));
        output_rawdata.len = positions[it->first].size();
        positions[it->first].clear();
        output_rawdata.gls = new double[likelihoods[it->first].size()];
        memcpy(output_rawdata.gls, likelihoods[it->first].data(), likelihoods[it->first].size() * sizeof(double));
        likelihoods[it->first].clear();
        output_rawdata.firstp = 0;
        output_rawdata.lastp = output_rawdata.len;
#if 0 //the code below should be read if we ever want to run on specific specified regions
        if (start != -1)
            while (output_rawdata.firstp < output_rawdata.len && output_rawdata.pos[output_rawdata.firstp] < start)
                output_rawdata.firstp++;


        if (stop != -1 && stop <= output_rawdata.pos[output_rawdata.lastp - 1]) {
            output_rawdata.lastp = output_rawdata.firstp;
            while (output_rawdata.pos[output_rawdata.lastp] < stop)
                output_rawdata.lastp++;
        }
#endif
        vcf_data[it->first] = std::move(output_rawdata);
    }
    fprintf(stderr, "\t-> VCF file successfully read\n");
    return vcf_data;
}