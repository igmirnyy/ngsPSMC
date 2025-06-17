#include <cstdio>
#include <htslib/vcf.h>
#include <htslib/tbx.h>
#include "header.h"
#include "bcf.h"

perBcf* perBcf_init(const char* fname){
    perBcf* pb = new perBcf;
    pb->bcf_file = bcf_open(fname, "r");
    if (!pb->bcf_file) {
        fprintf(stderr, "Failed to open BCF file: %s\n", fname);
        exit(1);
    }
    
    // Read the header
    pb->hdr = bcf_hdr_read(pb->bcf_file);
    if (!pb->hdr) {
        fprintf(stderr, "Failed to read BCF header from file: %s\n", fname);
        bcf_close(pb->bcf_file);
        exit(1);
    }
    
    // Load the index
    pb->idx = bcf_index_load(fname);
    if (!pb->idx) {
        fprintf(stderr, "Failed to read BCF index to file: %s\n", fname);
        bcf_close(pb->bcf_file);
        bcf_hdr_destroy(pb->hdr);
        exit(1);
    }
    pb->tbx = tbx_index_load(fname);
    if (!pb->tbx) {
        fprintf(stderr, "Failed to read BCF index to file: %s\n", fname);
        bcf_close(pb->bcf_file);
        bcf_hdr_destroy(pb->hdr);
        exit(1);
    }
    pb->rec = bcf_init1();
    if (!pb->rec){
        fprintf(stderr, "Failed to allocate BCF reader\n");
        bcf_close(pb->bcf_file);
        bcf_hdr_destroy(pb->hdr);
        hts_idx_destroy(pb->idx);
        exit(1);
    }
    pb->bcf_name = strdup(fname);
    return pb;

}
void perBcf_destroy(perBcf* pb){
    free(pb->bcf_name);
    bcf_hdr_destroy(pb->hdr);
    hts_idx_destroy(pb->idx);
    tbx_destroy(pb->tbx);
    bcf_destroy(pb->rec);
    bcf_close(pb->bcf_file);
}