use std::vec::Vec;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::io::{stdout, sink};
use anyhow::{Result, anyhow};

use bam2bedgraph::*;
use bam2bedgraph::indexed_annotation::*;

use structopt::StructOpt;

#[derive(StructOpt, Debug)]
#[structopt(name = "gff2bb", about = "Convert GFF/GTF to bigBed format")]
struct Options {
    // input files
    #[structopt(long="gff", help = "A genome annotation file in gff3 format", name="ANNOT_GFF_FILE")]
    annotfile_gff: Option<String>,
    #[structopt(long="gtf", help = "A genome annotation file in gtf format", name="ANNOT_GTF_FILE")]
    annotfile_gtf: Option<String>,
    #[structopt(long="chrmap", help = "Optional tab-delimited chr name mapping file", name="CHRMAP_FILE")]
    chrmap_file: Option<String>,
    #[structopt(long="vizchrmap", help = "Optional tab-delimited chr name mapping file for bigwig/bigbed exports", name="VIZCHRMAP_FILE")]
    vizchrmap_file: Option<String>,
    #[structopt(long="sizes", help = "Optional chr sizes file", name="SIZES_FILE")]
    sizes_file: Option<String>,
    // output file
    #[structopt(long="out", short="o", help = "Output file", name="OUT_FILE", default_value="-")]
    outfile: String,
    // feature types filter
    #[structopt(long="exon_type", help = "The exon type(s) to search for", name="EXON_TYPE")]
    exon_type: Vec<String>,
    #[structopt(long="transcript_type", help = "The transcript type(s) to search for", name="TRANSCRIPT_TYPE")]
    transcript_type: Vec<String>,
    #[structopt(long="gene_type", help = "The gene type(s) to search for", name="GENE_TYPE")]
    gene_type: Vec<String>,
    #[structopt(long="cds_type", help = "The CDS type(s) to search for", name="CDS_TYPE")]
    cds_type: Vec<String>,
    #[structopt(long="trackdb", help = "Write a UCSC trackDb.txt file with all the bigwigs/bigbeds", name="TRACKDB_FILE")]
    trackdb: Option<String>,
}

fn run() -> Result<()> {
    let mut options = Options::from_args();
    // set defaults for feature types
    options.gene_type =
        (if options.gene_type.is_empty() { vec!["gene".to_string()] } 
         else { options.gene_type.clone() }).into_iter().collect();
    //options.transcript_type =
    //    (if options.transcript_type.is_empty() { vec!["mRNA".to_string()] }
    //     else { options.transcript_type.clone() }).into_iter().collect();
    options.cds_type =
        (if options.cds_type.is_empty() { vec!["CDS".to_string()] }
         else { options.cds_type.clone() }).into_iter().collect();
    options.exon_type =
        (if options.exon_type.is_empty() { vec!["exon".to_string()] }
         else { options.exon_type.clone() }).into_iter().collect();
    // set debug options if --debug flag is set

    let transcript_type = String::from("transcript");
    let mut annot = if let Some(annotfile_gff) = options.annotfile_gff.clone() {
        eprintln!("Reading annotation file {:?}", &annotfile_gff);
        IndexedAnnotation::from_gff(
            &annotfile_gff, 
            &options.chrmap_file,
            &options.vizchrmap_file)?
    } else if let Some(annotfile_gtf) = options.annotfile_gtf.clone() {
        eprintln!("Reading annotation file {:?}", &annotfile_gtf);
        IndexedAnnotation::from_gtf(&annotfile_gtf, 
            options.gene_type.get(0).ok_or(anyhow!("NoneError"))?,
            options.transcript_type.get(0).unwrap_or(&transcript_type),
            &options.chrmap_file,
            &options.vizchrmap_file)?
    } else {
        Options::clap().print_help()?;
        eprintln!("\n\nNo annotation file was given!");
        std::process::exit(1);
    };
    if let Some(sizes_file) = options.sizes_file {
        annot.refs = read_sizes_file(&sizes_file, &annot.chrmap)?
    }

    // set up the trackdb writer
    let mut trackdb: BufWriter<Box<dyn Write>> = BufWriter::new(
        match options.trackdb.as_ref().map(String::as_ref) {
            Some("-") => Box::new(stdout()),
            Some(f) => Box::new(File::create(f)?),
            None => Box::new(sink()),
        });

    annot.to_bigbed(
        &options.outfile,
        &options.exon_type,
        &options.cds_type,
        &options.transcript_type,
        &options.gene_type,
        &mut trackdb)?;

    Ok(())
}

fn main() -> Result<()> {
    // enable stack traces
    std::env::set_var("RUST_BACKTRACE", "full");
    run()
}
