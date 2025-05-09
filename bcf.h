#include <htslib/vcf.h>


typedef struct {
    char* bcf_name;
    htsFile* bcf_file;
    bcf_hdr_t* hdr; //bcf header
    hts_idx_t* idx; //bcf index
  }perBcf;

perBcf* perBcf_init(const char* fname);
void perBcf_destroy(perBcf* pb);