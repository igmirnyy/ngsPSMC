#include <cstdio>
#include <htslib/vcf.h>
#include "header.h"
#include "bcf.h"

perBcf* perBcf_init(const char* fname){
    perBcf* pb = new perBcf;
    
    pb->sr = bcf_sr_init();
    if (!pb->sr) {
        fprintf(stderr, "Failed to create BCF reader file: %s\n", fname);
        exit(1);
    }
    
    if (!bcf_sr_add_reader(pb->sr, fname)) {
        fprintf(stderr, "Failed to open BCF file: %s\n", fname);
        exit(1);
    }
    
    // Read the header
    pb->hdr = bcf_sr_get_header(pb->sr, 0);
    if (!pb->hdr) {
        fprintf(stderr, "Failed to read BCF header from file: %s\n", fname);
        bcf_sr_destroy(pb->sr);
        exit(1);
    }
    
    // Load the index
    pb->idx = bcf_index_load(fname);
    if (!pb->idx) {
        fprintf(stderr, "Failed to read BCF index to file: %s\n", fname);
        bcf_sr_destroy(pb->sr);
        bcf_hdr_destroy(pb->hdr);
        exit(1);
    }
    pb->bcf_name = strdup(fname);
    return pb;

}
void perBcf_destroy(perBcf* pb){
    free(pb->bcf_name);
    bcf_hdr_destroy(pb->hdr);
    hts_idx_destroy(pb->idx);
    bcf_sr_destroy(pb->sr);
}