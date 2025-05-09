#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>


typedef struct {
    char* bcf_name;
    bcf_srs_t* sr;
    bcf_hdr_t* hdr; //bcf header
    hts_idx_t* idx; //bcf index
  }perBcf;

perBcf* perBcf_init(const char* fname);
void perBcf_destroy(perBcf* pb);