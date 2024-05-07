#include <stdio.h>
#include <htslib/vcf.h>

int main(int argc, char *argv[]) {
	
	if (argc != 4) {
        fprintf(stderr, "Usage: %s input.vcf output.vcf mnv_radius\n", argv[0]);
        return 1;
    }

    const char *input_file = argv[1];
    const char *output_file = argv[2];
	int mnv_radius = atoi(argv[3]);

    htsFile *fp_in = hts_open(input_file, "r");
    if (fp_in == NULL) {
	fprintf(stderr, "Error opening input vcf file: %s\n", input_file);
	return 1;
    } 

    bcf_hdr_t *header = bcf_hdr_read(fp_in);

    const char *format_id = "MNV";
    int format_id_index = bcf_hdr_id2int(header, BCF_DT_ID, format_id);
    if (format_id_index < 0) {
	char format_annotation[1024];
	snprintf(format_annotation, sizeof(format_annotation), "##FORMAT=<ID=%s,Number=1,Type=Integer,Description=\"Flag indicating presence of %s\">", format_id, format_id);

	int ret = bcf_hdr_append(header, format_annotation);
	if (ret < 0) {
	    fprintf(stderr, "Error adding new FORMAT annotation %s\n", format_id);
	    bcf_hdr_destroy(header);
	    hts_close(fp_in);
	    return 1;
	}


	format_id_index = bcf_hdr_id2int(header, BCF_DT_ID, format_id);
	if (format_id_index < 0) {
	    fprintf(stderr, "Error getting index for %s\n", format_id);
	    bcf_hdr_destroy(header);
	    hts_close(fp_in);
	    return 1;
	}
    }

    htsFile *fp_out = hts_open(output_file, "w");
    if (fp_out == NULL) {
	fprintf(stderr, "Error opening output VCF file: %s\n", output_file);
	return 1;
    }

    if (bcf_hdr_write(fp_out, header) != 0) {
	fprintf(stderr, "Error writing header to output VCF file: %s\n", output_file);
	bcf_hdr_destroy(header);
	hts_close(fp_in);
	hts_close(fp_out);
	return 1;
    }


    int nsamples = bcf_hdr_nsamples(header); 
    //int mnv_radius = 5;
    int max_stored = mnv_radius + 1;
    //int32_t *mnv_cr = calloc(nsamples, sizeof(int32_t));
    //int32_t **mnv_array = calloc(mnv_radius, sizeof(int32_t *));
    int32_t mnv_cr[nsamples];
    memset(mnv_cr, bcf_int32_missing, nsamples * sizeof(int32_t));
    int32_t mnv_array[max_stored][nsamples];


    bcf1_t *cr = bcf_init();
    bcf1_t **prs = (bcf1_t **) malloc(max_stored * sizeof(bcf1_t *));
    for (int i = 0; i < max_stored; i++) {
	prs[i] = bcf_init();
    }
    int n_prs = 0;
 
    // TODO: Check if using bitfield variables decreases runtime? It might have a substantial impact on the combined duration of memset calls.
    int i, j, missing;
    int ngt_arr = 0, *gt_arr = NULL, ngt;

    while(bcf_read(fp_in, header, cr) >= 0) {
	// Clear stored records if they get out of the mnv_radius
	for (i = n_prs-1; i >= 0; i--) {
	    if ((cr->rid != prs[i]->rid) || ((cr->pos - mnv_radius) > prs[i]->pos)) {
		ngt = bcf_get_genotypes(header, prs[i], &gt_arr, &ngt_arr);
		for (j = 0; j < nsamples; j++) {
		    mnv_array[i][j] = bcf_gt_allele(gt_arr[j * 2]) < 0 ? bcf_int32_missing : mnv_array[i][j];
		}

		bcf_update_format_int32(header, prs[i], format_id, &mnv_array[i], nsamples);
		bcf_write(fp_out, header, prs[i]);
		bcf_clear(prs[i]);
		n_prs--;
	    } else {
		ngt = bcf_get_genotypes(header, cr, &gt_arr, &ngt_arr);

		for (j = 0; j < nsamples; j++) {
		    missing = bcf_gt_allele(gt_arr[j * 2]) + bcf_gt_allele(gt_arr[j * 2 + 1]);
		    if (missing > 0) {
			mnv_cr[j] = 1;
			mnv_array[i][j] = 1;
		    }
		}
	    }
	}

	for (i = n_prs; i > 0; i--) {
	    bcf_copy(prs[i], prs[i-1]);
	    memcpy(mnv_array[i], mnv_array[i-1], nsamples * sizeof(int32_t));
	}

	bcf_copy(prs[0], cr);
	memcpy(mnv_array[0], mnv_cr, nsamples * sizeof(int32_t));
	memset(mnv_cr, bcf_int32_missing, nsamples * sizeof(int32_t));
	
	n_prs++;
	bcf_clear(cr);
    }
    bcf_destroy(cr);

    // Close out trailing records leftover from the loop
    for (i = n_prs-1; i>= 0; i--) {
	bcf_update_format_int32(header, prs[i], format_id, &mnv_array[i], nsamples);
	bcf_write(fp_out, header, prs[i]);
	bcf_clear(prs[i]);
	n_prs--;
    }

    for (i = 0; i < mnv_radius; i++) {
	bcf_destroy(prs[i]);
    } 
    free(prs); 
    bcf_hdr_destroy(header);
    hts_close(fp_in);
    hts_close(fp_out);

    return 0;
}
