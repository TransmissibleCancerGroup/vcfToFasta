import std.algorithm;
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

auto examineCoverage(string sampleData, int min_total_cov, int min_alt_cov) {
    // Assuming sample data format is GT:GL:GOF:GQ:NR:NV
    auto splitLine = sampleData.splitter(':').drop(4);

    int nr = to!int(splitLine.front); splitLine.popFront(); // parse out 2nd-last field
    int nv = to!int(splitLine.front); splitLine.popFront(); // parse out last field

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

void processVcfLine(ref Sequence[] sequences, ref Sequence reference, string line,
        bool excludeInvariant, int mintot, int minalt, bool gInfo) {
    auto splitLine = line.splitter("\t");
    auto chrom = splitLine.front; splitLine.popFront();
    auto pos = splitLine.front; splitLine.popFront();
    splitLine = splitLine.drop(1);
    string refBase = splitLine.front; splitLine.popFront();
    auto altBase = splitLine.front;

    // add a collector here to check if this site is invariant
    auto site = appender!(string)();

    // For each sample, examine its genotype and decide if it is
    // REF or ALT. Add the choice to the site
    foreach (genotypeData; line.splitter("\t").drop(9)) {
        VariantResult result;
        if (gInfo) result = examineGenotype(genotypeData);
        else result = examineCoverage(genotypeData, mintot, minalt);
        switch (result) {
            default:
                site ~= refBase;
                break;
            case VariantResult.IsRef:
                site ~= refBase;
                break;
            case VariantResult.IsAlt:
                site ~= altBase;
                break;
            case VariantResult.IsAmbiguous:
                site ~= refBase;
                break;
        }
    }

    // Have a look at the site, check if it's invariant. Remember to check the reference, too.
    bool siteIsInvariant = array(chain(refBase, site.data)).sort().uniq.count == 1;

    // If excludeInvariant is true, and siteIsInvariant, then don't add
    // (return early)
    if (excludeInvariant && siteIsInvariant) {
        return;
    }

    // Otherwise, add the site to all sequences including the reference
    reference.seq ~= refBase;
    foreach (ref sequence, character; lockstep(sequences, site.data)) {
        sequence.seq ~= character;
    }
}

struct Sequence {
    string name;
    Appender!(string) seq;
}

unittest {
    string line = "\t\t\tT\tC\t\t\t\tGT:GL:GOF:GQ:NR:NV\t1/1::::8:0\t1/1::::15:6";
    Sequence[] seqs = [Sequence("a"), Sequence("b")];
    Sequence refseq = Sequence("ref");
    processVcfLine(seqs, refseq, line, false, 10, 5, true);
    assert(seqs[0].seq.data[0] == 'C');
    assert(seqs[1].seq.data[0] == 'C');
    assert(refseq.seq.data[0] == 'T');

    processVcfLine(seqs, refseq, line, false, 10, 5, false);
    assert(seqs[0].seq.data[1] == 'T');
    assert(seqs[1].seq.data[1] == 'C');
    assert(refseq.seq.data[1] == 'T');

    processVcfLine(seqs, refseq, line, false, 100, 50, false);
    assert(seqs[0].seq.data[2] == 'T');
    assert(seqs[1].seq.data[2] == 'T');
    assert(refseq.seq.data[2] == 'T');

    processVcfLine(seqs, refseq, line, true, 100, 50, false);
    assert(seqs[0].seq.data.length == 3);
    assert(seqs[1].seq.data.length == 3);
    assert(refseq.seq.data.length == 3);
}

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
    return run(options);
}

int run(Options options) {

    if (!exists(options.infilename)) {
        stderr.writefln("Couldn't open %s", options.infilename);
        return 1;
    }

    auto file = File(options.infilename, "r");
    scope(exit) file.close();

    Sequence[] sequences;

    // add reference seq
    auto refseq = Sequence("Reference");

    bool initialised = false;

    ulong counter = 0;

    foreach (line; file.byLine()) {
        if (line.startsWith("##")) continue;
        try {
            if (line.startsWith("#") && !initialised) {
                sequences = getSequences(to!string(line));
                initialised = true;
            }
            else {
                processVcfLine(sequences, refseq, to!string(line),
                    options.excludeInvariant, options.min_total_cov,
                    options.min_alt_cov, options.useGenotypeInfo);
            }
        }
        catch (Throwable) {
            writefln("Error parsing %s as VCF. Check input and try again.",
                     options.infilename);
            return 1;
        }
        counter++;
        if (counter % 100000 == 0) {
            stderr.writefln("Processed %d SNPs", counter);
            stderr.flush();
        }
        // if (counter > 1000000) break;
    }

    foreach (ref sequence; chain([refseq], sequences)) {
        writefln(">%s", sequence.name);
        foreach (ref chunk; sequence.seq.data.chunks(80)) {
            stdout.writeln(chunk);
        }
        stdout.flush();
    }
    return 0;
}

