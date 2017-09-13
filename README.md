# vcfToFasta
Convert SNP data in a somatypus VCF file to an alignment in Fasta format

### Compilation
This is a D program, and can be built using dub. The binary in the repo was built on ```Ubuntu 12.04.5 LTS``` using the command
    
    dub build --compiler=ldc2 --build=release

### Basic usage
    Usage: vcfToFasta [--help] [--excludeInvariant] [--useGenotypeInfo]
                    [--min_total_cov=<uint>] [--min_alt_cov=<uint>] vcffile

    Positional arguments:
    vcffile         VCF file (SNPs) to convert to alignment

    Optional arguments:
    --help, -h      This help information
    --excludeInvariant, -e
                    Filter out invariant sites from the alignment
    --useGenotypeInfo, -g
                    Use the genotype information in the VCF file to confirm variant
                    sites, rather than numerical filters (default = off)
    --min_total_cov, -c <uint>
                    Minimum number of reads (total) needed to consider as a
                    potentially variant site (default = 10)
    --min_alt_cov, -a <uint>
                    Minimum number of variant-containing reads needed to confirm as
                    a variant site (default = 5)
