## Running annotate_mnv

```annotate_mnv``` is written in C, and uses HTSLib functions. HTSLib is required to compile and execute.

If an HTSLib module is already available on your system, it should be loaded before compiling or running the executable:

```
module spider htslib
module load bioinf/htslib/1.9
```

If HTSLib is unavailable on your system, you can [download HTSLib here](https://www.htslib.org/download/) and [install from source](https://github.com/samtools/htslib/blob/develop/INSTALL). If ```annotate_af``` can't find HTSLib at runtime, execution will fail with ```./annotate_mnv: error while loading shared libraries: libhts.so.2: cannot open shared object file: No such file or directory```. I compile and run with HTSLib 1.9; other versions of HTSLib have not been tested and may cause unexpected behaviour.

```annotate_mnv``` expects three arguments:
* An input VCF file (- for STDIN); these can be compressed or uncompressed, indexed or unindexed
* An output VCF path (- for STDOUT)
* An mnv_radius; the distance between two variants for them to be considered MNVs (in base pairs)

To test whether ```annotate_mnv``` is working correctly, you can run the following commands:

```
./annotate_mnv test/test.vcf.gz test/test_mnv_annotated.vcf
./annotate_mnv test/test_merged.vcf.gz test/test_merged_mnv_annotated.vcf
```

If you want to receive input from STDIN or write to STDOUT as part of a pipe, you can use the following syntax:

```
cat test/test.vcf.gz | ./annotate_mnv - - > test/test_mnv_annotated.vcf
cat test/test_merged.vcf.gz | ./annotate_mnv - - > test/test_merged_mnv_annotated.vcf
```

## Compiling annotate_af

If you need to compile ```annotate_mnv``` from source, first load HTSLib. You can then compile the executable using the following command:

```gcc -o annotate_mnv vcf_annotate_mnv.c -lhts -lz```

Compilation was done using gcc 10.2.0 on Linux.
