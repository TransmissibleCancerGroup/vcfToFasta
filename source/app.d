import std.algorithm;
import std.algorithm: sort;
import std.conv;
import std.file;
import std.path: baseName;
import std.range;
import std.stdio;
import std.string;
import darg;

/*
Script to extract SNP information from a VCF and make an alignment
*/


Sequence[] getSequences(string header) {
    auto sequences = appender!(Sequence[]);
    foreach (name; header.chompPrefix("#").splitter("\t").drop(9)) {
        sequences.put(Sequence(name));
    }
    return sequences.data;
}

enum VariantResult {
    IsRef,
    IsAlt,
    IsAmbiguous
}

enum string[string] ambiguityCodes = [
    "AC": "M", "AG": "R", "AT": "W",
    "CA": "M", "CG": "S", "CT": "Y",
    "GA": "R", "GC": "S", "GT": "K",
    "TA": "W", "TC": "Y", "TG": "K"
];

string getAmbiguityCode(string chars)
in {
    assert(chars.length == 2);
}
body {
    return ambiguityCodes[chars];
}

unittest {
    static assert(getAmbiguityCode("AC") == "M");
    string refBase = "A";
    string altBase = "T";
    string ambig = getAmbiguityCode(refBase ~ altBase);
    assert(ambig == "W");
    assert(getAmbiguityCode("TA") == "W");
}

// Reduce 2 Variant results to 1 - if they differ, result is ambiguous
// if they are the same, output result is same type
auto composeResults(VariantResult a, VariantResult b) {
    if (a != b) return VariantResult.IsAmbiguous;
    else return a;
}

unittest {
    auto refnce = VariantResult.IsRef;
    auto alt = VariantResult.IsAlt;
    auto ambig = VariantResult.IsAmbiguous;

    assert(composeResults(refnce, refnce) == refnce);
    assert(composeResults(refnce, alt) == ambig);
    assert(composeResults(refnce, ambig) == ambig);
    assert(composeResults(alt, ambig) == ambig);
    assert(composeResults(ambig, ambig) == ambig);
    assert(composeResults(alt, alt) == alt);
    assert(composeResults(alt, refnce) == ambig);
}

auto examineGenotype(string genotype) {
    if (genotype.startsWith("0/0")) {
        return VariantResult.IsRef;
    }
    else if (genotype.startsWith("./.")) {
        return VariantResult.IsAmbiguous;
    }
    else {
        return VariantResult.IsAlt;
    }
}

unittest {
    assert(examineGenotype("1/1:::::") == VariantResult.IsAlt);
    assert(examineGenotype("0/1:::::") == VariantResult.IsAlt);
    assert(examineGenotype("1/0:::::") == VariantResult.IsAlt);
    assert(examineGenotype("0/0:::::") == VariantResult.IsRef);
    assert(examineGenotype("./.:::::") == VariantResult.IsAmbiguous);
}

// Consume the front item of a range
auto pop(Range)(ref Range r) if (isInputRange!Range) {
    if (r.empty) {
        ElementType!Range val;
        return val;
    }
    auto f = r.front;
    r.popFront();
    return f;
}

unittest {
    string s = "a\tb\tc\td\te";
    auto spl = s.splitter('\t');
    string f = pop(spl);
    assert(f == "a");
    assert(spl == s.splitter('\t').drop(1));
    string g = pop(spl);
    assert(g == "b");
    assert(spl == s.splitter('\t').drop(2));

    s = "";
    spl = s.splitter('\t');
    f = pop(spl);
    assert(spl.empty);
    assert(f == "");
}

auto examineCoverage(string sampleData, int min_total_cov, int min_alt_cov) {
    // Assuming sample data format is GT:GL:GOF:GQ:NR:NV
    auto splitLine = sampleData.splitter(':').drop(4);

    int nr = to!int(pop(splitLine)); // parse out 2nd-last field
    int nv = to!int(pop(splitLine)); // parse out last field

    if (nr < min_total_cov) {
        return VariantResult.IsAmbiguous;
    }
    else if (nv < min_alt_cov) {
        return VariantResult.IsRef;
    }
    else {
        return VariantResult.IsAlt;
    }
}

unittest {
    string data = "1/1:0.0,-2,03:45:20:8:0";
    assert(examineCoverage(data, 10, 5) == VariantResult.IsAmbiguous);
    assert(examineCoverage("::::8:0", 10, 5) == VariantResult.IsAmbiguous);
    assert(examineCoverage("::::8:8", 10, 5) == VariantResult.IsAmbiguous);
    assert(examineCoverage("::::10:0", 10, 5) == VariantResult.IsRef);
    assert(examineCoverage("::::10:8", 10, 5) == VariantResult.IsAlt);
}


enum int[dchar] baseToBits = [
    'A': 0x1, 'C': 0x2, 'G': 0x4, 'T': 0x8,
    'M': 0x1|0x2, 'R': 0x1|0x4, 'W': 0x1|0x8,
    'S': 0x2|0x4, 'Y': 0x2|0x8, 'K': 0x4|0x8
];

// Checks if a site (i.e. an alignment column) is invariant (i.e. has no phylogenetic information)
// -NO LONGER USED-
bool isSiteInvariant(string site)
in {
    assert(!site.empty);
}
body {
    int front = baseToBits[site[0]];
    foreach (c; site) {
        if (!(baseToBits[c] & front)) return false;
    }
    return true;
}

unittest {
    static assert(isSiteInvariant("TKKKKKKKKKKTKKKKKKKKTKKTTTTKTTTTK"));
    assert(isSiteInvariant("TYTYYYYYYTTTTYYYTYYTTYYTTTTYTTYTY"));
    assert(isSiteInvariant("TWWTTTTWTTTWTWWWWWWTTWWWTTTWWWTTT"));
    assert(isSiteInvariant("AAAAAAAAAAAMAMAMAAAAAMAAAAAAAAAAA"));
    assert(isSiteInvariant("MCMCCCMAAAACAACACCCM"));
    assert(isSiteInvariant("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"));
    assert(!isSiteInvariant("ACAAAAAAAA"));
}

// Parse a list of key-value pairs into a hash map.
// Input format is "K1=V1;K2=V2;..."
auto parseKeyValue(string input) {
    string[string] hashmap;
    foreach (pair; input.splitter(";")) {
        auto keyVal = pair.splitter("=").array;
        hashmap[keyVal[0]] = keyVal[1];
    }
    return hashmap;
}

unittest {
    string[string] test = [
        "QD": "20",
        "SC": "TGTGTAAGCATTGGGCCTGAA",
        "SbPval": "0"
    ];
    auto info = "QD=20;SC=TGTGTAAGCATTGGGCCTGAA;SbPval=0";
    assert(parseKeyValue(info) == test);
    assert(parseKeyValue(info)["SC"] == test["SC"]);
}

auto makeTripletList() {
    auto alphabet = "ACGT";
	return cartesianProduct(alphabet, alphabet, alphabet)
		.map!(a => format("%s%s%s", a[0], a[1], a[2])).array;
}

string compressTriplet(string triplet) {
    static tripletList = makeTripletList;
    return to!string(to!char(countUntil(tripletList, triplet) + 63));
}

unittest {
    static ctTripletList = makeTripletList;
    auto rtTripletList = makeTripletList;
    assert(ctTripletList == rtTripletList);
    assert(compressTriplet("TCG") == "u");
    assert(compressTriplet("ACG") == "E");
}

// TODO: reduce number of params by passing cmdline options struct directly
void processVcfLine(ref Sequence[] sequences, ref Sequence reference, string line,
        bool excludeInvariant, int mintot, int minalt, bool gInfo, bool ambiguityIsRef,
        bool fullContext)
{
    // Grab what we want from the tab-separated fields
    auto splitLine = line.splitter("\t");
    auto chrom = pop(splitLine);
    auto pos = pop(splitLine);
    splitLine = splitLine.drop(1); // ID
    string refBase = pop(splitLine);
    auto altBase = pop(splitLine);

    // Examine sequence context
    string pre, post;
    if (fullContext) {
        splitLine = splitLine.drop(2); // QUAL, FILTER
        auto info = pop(splitLine);
        auto infoHashMap = parseKeyValue(info);
        auto context = infoHashMap["SC"];
        auto centre = (context.length - 1) / 2;
        pre = to!string(context[centre-1]);
        post = to!string(context[centre+1]);
    }

    // add a collector here to check if this site is invariant
    auto site = appender!(string)();

    // This bool monitors whether we ever see an ALT residue
    bool isVariantSite = false;

    // For each sample, examine its genotype and decide if it is
    // REF or ALT. Add the choice to the site
    foreach (genotypeData; line.splitter("\t").drop(9)) {
        VariantResult result;
        if (gInfo) result = examineGenotype(genotypeData);
        else result = examineCoverage(genotypeData, mintot, minalt);
        string nextSite;
        switch (result) {
            default:
                nextSite = refBase;
                break;
            case VariantResult.IsRef:
                nextSite = refBase;
                break;
            case VariantResult.IsAlt:
                isVariantSite = true;
                nextSite = altBase;
                break;
            case VariantResult.IsAmbiguous:
                if (ambiguityIsRef) nextSite = refBase;
                else nextSite = getAmbiguityCode(refBase ~ altBase);
                break;
        }

        // If using base triplet compression then we need to convert
        // our current refBase and altBase. If not, the compressedRefBase
        // is the same as refBase
        if (fullContext) {
            nextSite = compressTriplet(pre ~ nextSite ~ post);
        }
        site ~= nextSite;
    }

    // Compress the refBase
    if (fullContext) {
        refBase = compressTriplet(pre ~ refBase ~ post);
    }

    // If excludeInvariant is true, and site is invariant, then don't add
    // (return early)
    if (excludeInvariant && !isVariantSite) {
        return;
    }

    // Otherwise, add the site to all sequences including the reference
    reference.seq ~= refBase;
    foreach (ref sequence, character; lockstep(sequences, site.data)) {
        sequence.seq ~= character;
    }
}


// This struct represents a sequence. It holds the sequence name (`name`)
// and an appender (`seq`) to hold the growing sequence
struct Sequence {
    string name;
    Appender!(string) seq;
}

unittest {
    string line = "\t\t\tT\tC\t\t\tSC=TGTTTCG\tGT:GL:GOF:GQ:NR:NV\t1/1::::8:0\t1/1::::15:6";
    Sequence[] seqs = [Sequence("a"), Sequence("b")];
    Sequence refseq = Sequence("ref");

    // Use genotype strings - both sites are alt
    processVcfLine(seqs, refseq, line, false, 10, 5, true, true, false);
    assert(seqs[0].seq.data[0] == 'C');
    assert(seqs[1].seq.data[0] == 'C');
    assert(refseq.seq.data[0] == 'T');

    // Use coverage counts - site 1 is ref, 2 is alt
    processVcfLine(seqs, refseq, line, false, 8, 5, false, true, false);
    assert(seqs[0].seq.data[1] == 'T');
    assert(seqs[1].seq.data[1] == 'C');
    assert(refseq.seq.data[1] == 'T');

    // Use coverage counts - both sites are ref at these thresholds
    // this site is invariant, but not filtered out
    processVcfLine(seqs, refseq, line, false, 8, 8, false, true, false);
    assert(seqs[0].seq.data[2] == 'T' && seqs[0].seq.data.length == 3);
    assert(seqs[1].seq.data[2] == 'T' && seqs[1].seq.data.length == 3);
    assert(refseq.seq.data[2] == 'T' && refseq.seq.data.length == 3);

    // As above, but filtered out (array is same size)
    processVcfLine(seqs, refseq, line, true, 8, 8, false, true, false);
    assert(seqs[0].seq.data.length == 3);
    assert(seqs[1].seq.data.length == 3);
    assert(refseq.seq.data.length == 3);

    // Compressed sequence context
    processVcfLine(seqs, refseq, line, false, 10, 5, true, true, true);
    assert(seqs[0].seq.data[3] == 'v');
    assert(seqs[1].seq.data[3] == 'v');
    assert(refseq.seq.data[3] == '~');
}


// This struct holds command line arguments
struct Options {
    @Option("help", "h")
    @Help("This help information")
    OptionFlag help;

    @Option("excludeInvariant", "e")
    @Help("Filter out invariant sites from the alignment")
    OptionFlag excludeInvariant;

    @Option("useGenotypeInfo", "g")
    @Help("Use the genotype information in the VCF file to confirm "~
            "variant sites, rather than numerical filters "~
            "(default = off)")
    OptionFlag useGenotypeInfo;

    @Option("ambiguityIsRef", "r")
    @Help("Ambiguous calls are conservatively called as the reference base. "~
            "If switched off, sites are called as the IUPAC ambiguity code "~
            "standing for Ref/Base (default = off)")
    OptionFlag ambiguityIsRef;

    @Option("fullContext", "f")
    @Help("EXPERIMENTAL: output in a format that carries triplet sequence context information. "~
            "Each triplet is encoded as its position in the alphabetical list [AAA, AAC, AAG, ... TTT] "~
            "+63, converted to ascii. Not compatible with ambiguity codes, so enforces ambiguityIsRef "~
            "to be switched on when fullContext is active")
    OptionFlag fullContext;

    @Option("min_total_cov", "c")
    @Help("Minimum number of reads (total) needed to consider "~
            "as a potentially variant site (default = 10)")
    uint min_total_cov = 10;

    @Option("min_alt_cov", "a")
    @Help("Minimum number of variant-containing reads needed to "~
            "confirm as a variant site (default = 5)")
    uint min_alt_cov = 5;

    @Argument("vcffile")
    @Help("VCF file (SNPs) to convert to alignment")
    string infilename;
}

// Generate the usage and help string at compile time.
immutable usage = usageString!Options("vcfToFasta ");
immutable help = helpString!Options;

// main handles the commandline and dispatches all work to `run`
int main(string[] args)
{
    Options options;

    try {
        options = parseArgs!Options(args[1 .. $]);
    }
    catch (ArgParseError e) {
        writeln(e.msg);
        writeln(usage);
        return 1;
    }
    catch (ArgParseHelp e) {
        writeln(usage);
        write(help);
        return 0;
    }
    if (options.fullContext == OptionFlag.yes) options.ambiguityIsRef = OptionFlag.yes;
    return run(options);
}

// Does all the work
int run(Options options) {
    // Won't get too far if the infile doesn't exist
    if (!exists(options.infilename)) {
        stderr.writefln("Couldn't open %s", options.infilename);
        return 1;
    }

    auto file = File(options.infilename, "r");
    scope(exit) file.close();

    // This array will hold the growing sequence info for each sample
    Sequence[] sequences;

    // Also keep track of the reference sequence
    auto refseq = Sequence("Reference");

    // Need to initialise the sequences only once - when we reach
    // the VCF column headers line, which starts with a single '#'
    bool initialised = false;

    ulong loopcounter = 0; // Just for monitoring progress

    // Loop over the VCF file and process each line
    foreach (line; file.byLine()) {
        if (line.startsWith("##")) {
            // Ignore these lines
            continue;
        }
        try {
            if (line.startsWith("#") && !initialised) {
                // Reached the headers - initialise data structure
                sequences = getSequences(to!string(line));
                initialised = true;
            }
            else {
                // Process a line and extract next sequence residue
                processVcfLine(sequences, refseq, to!string(line),
                    options.excludeInvariant, options.min_total_cov,
                    options.min_alt_cov, options.useGenotypeInfo,
                    options.ambiguityIsRef, options.fullContext);
            }
        }
        catch (Throwable e) {
            writefln("Error parsing %s as VCF. Check input and try again.",
                     options.infilename);
            // throw(e);
            return 1;
        }

        // Rudimentary progress monitor
        loopcounter++;
        if (loopcounter % 100000 == 0) {
            stderr.writefln("Processed %d SNPs", loopcounter);
            stderr.flush();
        }
    }
    stderr.writefln("Processed %d SNPs in total.", loopcounter);
    stderr.flush();

    // Write Fasta to standard out
    foreach (ref sequence; chain([refseq], sequences)) {
        writefln(">%s", sequence.name);
        foreach (ref chunk; sequence.seq.data.chunks(80)) {
            stdout.writeln(chunk);
        }
        stdout.flush();
    }
    return 0;
}

