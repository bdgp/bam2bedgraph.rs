#![recursion_limit = "1024"]
#![cfg_attr(feature = "cargo-clippy", allow(cyclomatic_complexity, trivial_regex))]
use std::ascii::AsciiExt;
use std::collections::BTreeMap;
use std::io::Write;
use std::io::BufWriter;
use std::io::stderr;
use std::fs::File;
use std::str;
use std::path::{PathBuf, Path};
use std::process::Command;
use std::vec::Vec;

extern crate regex;
use regex::Regex;

extern crate theban_interval_tree;
extern crate memrange;
use memrange::Range;

extern crate rust_htslib;
use rust_htslib::bam::Read;
use rust_htslib::bam::Reader;

extern crate structopt;
#[macro_use]
extern crate structopt_derive;

use structopt::StructOpt;

extern crate bam2bedgraph;
use bam2bedgraph::errors::*;
use bam2bedgraph::errors::ToResult;
use bam2bedgraph::cigar2exons;

fn open_file(options: &Options,
             read_number: i32,
             strand: &str,
             split_strand: &str,
             fhs: &mut BTreeMap<String, Option<File>>)
             -> Result<String> {
    let mut prefix = PathBuf::new();
    prefix.set_file_name(&options.bamfile);
    prefix.set_extension("");
    let track_name = vec![if !options.trackname.is_empty() {
                              options.trackname.clone()
                          } else {
                              let p = prefix.as_path().to_str().ok_or("Failed to parse path!")?;
                              p.to_string()
                          },
                          if options.split_read && read_number > 0 {
                              format!(".r{}", read_number)
                          } else {
                              "".to_string()
                          },
                          if split_strand != "uu" && !strand.is_empty() {
                              format!(".{}", strand)
                          } else {
                              "".to_string()
                          }]
        .join("");

    let filename = vec![if !options.out.is_empty() {
                            options.out.clone()
                        } else {
                            let p = prefix.as_path().to_str().ok_or("Failed to parse path!")?;
                            p.to_string()
                        },
                        if options.split_read && read_number > 0 {
                            format!(".r{}", read_number)
                        } else {
                            "".to_string()
                        },
                        if split_strand != "uu" && !strand.is_empty() {
                            format!(".{}", strand)
                        } else {
                            "".to_string()
                        },
                        ".bedgraph".to_string()]
        .join("");

    // initialize the file if needed
    if !fhs.contains_key(&filename) {
        let mut f = File::create(&filename)?;
        if options.trackline && !options.bigwig {
            writeln!(f,
                     "track type=bedGraph name=\"{}\" description=\"{}\" visibility=full",
                     track_name,
                     track_name)?;
        }
        fhs.insert(filename.clone(), Some(f));
    }
    Ok(filename)
}

fn write_chr(options: &Options,
             chr: &(u32, String),
             histogram: &BTreeMap<(i32, String), Vec<i32>>,
             fhs: &mut BTreeMap<String, Option<File>>,
             split_strand: &str)
             -> Result<()> {
    for (key, histo) in histogram {
        let read_number = key.0;
        let strand = &key.1;
        let filename = open_file(options, read_number, strand, split_strand, fhs)?;
        let mut f = fhs.get_mut(&filename).r()?;
        let file = f.as_mut().r()?;
        let mut writer = BufWriter::new(file);

        // scan the histogram to produce the bedgraph data
        let mut start: usize = 0;
        let mut end: usize = 0;
        let ref_length: usize = chr.0 as usize;
        while start < ref_length {
            while (if end < histo.len() { histo[end] } else { 0 }) ==
                  (if start < histo.len() { histo[start] } else { 0 }) &&
                  end < ref_length {
                end += 1
            }
            if options.zero || (if start < histo.len() { histo[start] } else { 0 }) > 0 {
                writeln!(writer,
                         "{}\t{}\t{}\t{}",
                         chr.1,
                         start,
                         end,
                         if strand == "-" {
                             -histo[start]
                         } else {
                             histo[start]
                         })?;
            }
            start = end;
        }
    }
    Ok(())
}

fn analyze_bam(options: &Options,
               split_strand: &str,
               autostrand_pass: bool,
               intervals: &Option<BTreeMap<String, theban_interval_tree::IntervalTree<u8>>>)
               -> Result<()> {
    if !Path::new(&options.bamfile).exists() {
        return Err(format!("Bam file {} could not be found!", &options.bamfile).into());
    }
    let bam = Reader::from_path(&options.bamfile)?;
    let header = bam.header();

    let mut refs: Vec<(u32, String)> = Vec::new();
    refs.resize(header.target_count() as usize, (0, "".to_string()));
    let target_names = header.target_names();
    for target_name in target_names {
        let tid = header.tid(target_name).r()?;
        let target_len = header.target_len(tid).r()?;
        let target_name = std::str::from_utf8(target_name)?;
        refs[tid as usize] = (target_len, target_name.to_string());
    }

    if options.fixchr {
        for r in &mut refs {
            let regex = Regex::new(r"^(chr|Zv9_)")?;
            if regex.is_match(&r.1) {
                let refname = r.1.to_string();
                r.1.clear();
                r.1.push_str(&format!("chr{}", refname));
            }
        }
    }
    if autostrand_pass {
        writeln!(stderr(),
                 "Running strand detection phase on {}",
                 options.bamfile)?;
    } else {
        writeln!(stderr(), "Building histograms for {}", options.bamfile)?;
    }

    // build a lookup map for the refseqs
    let mut refmap: BTreeMap<String, usize> = BTreeMap::new();
    for (i, _) in refs.iter().enumerate() {
        refmap.insert(refs[i].1.to_string(), i);
    }

    let mut lastchr: i32 = -1;
    let mut fhs: BTreeMap<String, Option<File>> = BTreeMap::new();
    let mut histogram: BTreeMap<(i32, String), Vec<i32>> = BTreeMap::new();

    let mut autostrand_totals: BTreeMap<char, i64> = BTreeMap::new();
    autostrand_totals.insert('s', 0);
    autostrand_totals.insert('r', 0);
    let mut autostrand_totals2: BTreeMap<char, i64> = BTreeMap::new();
    autostrand_totals2.insert('s', 0);
    autostrand_totals2.insert('r', 0);

    let mut read = rust_htslib::bam::record::Record::new();
    while bam.read(&mut read).is_ok() {
        // if we've hit a new chr, write out the bedgraph data and clear the histogram
        if lastchr == -1 || read.tid() != lastchr {
            if !autostrand_pass && !histogram.is_empty() && lastchr != -1 {
                write_chr(options,
                          &refs[lastchr as usize],
                          &histogram,
                          &mut fhs,
                          split_strand)?;
            }
            histogram.clear();
            lastchr = read.tid();
        }

        // skip this read if it's no good
        let paired = read.is_paired();
        let proper = read.is_proper_pair();
        let primary = !read.is_secondary();
        if (options.paired_only && !paired) || (options.primary_only && !primary) ||
           (options.proper_only && !proper) {
            continue;
        }

        // skip if it's not unique and we want unique alignments
        if options.uniq {
            let hits = read.aux("NH".to_string().as_bytes());
            if hits == None || hits.r()?.integer() != 1 {
                continue;
            }
        }

        let mut exons: Vec<(u64, u64)> = Vec::new();
        let mut get_exons: Vec<(u64, u64)> = Vec::new();
        cigar2exons(&mut get_exons, &read.cigar(), read.pos() as u64)?;
        if !options.split_exons && !get_exons.is_empty() {
            let first = get_exons.get(0).r()?;
            let last = get_exons.get(get_exons.len() - 1).r()?;
            exons = vec![(first.0, last.1)];
        } else {
            exons.append(&mut get_exons);
        }
        let read_number = if read.is_last_in_template() { 2 } else { 1 };

        // attempt to determine the strandedness of the transcript
        // read numbers match, is not reverse, is not flipped
        // let xs = read.aux("XS".as_bytes());
        let strand = 
        //if xs.is_some() {
         //   str::from_utf8(xs.r()?.string())?
        //} else 
        if read_number == 1 {
                if split_strand.chars().nth(0).r()? == 'r' {
                    if read.is_reverse() { "+" } else { "-" }
                } else if split_strand.chars().nth(0).r()? == 's' {
                    if read.is_reverse() { "-" } else { "+" }
                } else {
                    ""
                }
            } else if read_number == 2 {
                if split_strand.chars().nth(1).r()? == 's' {
                    if read.is_reverse() { "-" } else { "+" }
                } else if split_strand.chars().nth(1).r()? == 'r' {
                    if read.is_reverse() { "+" } else { "-" }
                } else {
                    ""
                }
            } else {
                ""
            };

        let read_num = if options.split_read { read_number } else { 0 };

        let ref_length = refs[read.tid() as usize].0;

        // add the read to the histogram
        for exon in exons {
            // try to determine the strandedness of the data
            if autostrand_pass {
                if intervals.is_some() {
                    let intervals = intervals.as_ref().r()?;
                    if intervals.contains_key(&refs[lastchr as usize].1) {
                        for (_, interval) in intervals[&refs[lastchr as usize].1]
                            .range(exon.0 + 1u64, exon.1)
                            .enumerate() {
                            let overlap_length = std::cmp::min(exon.1, interval.0.max) -
                                                 std::cmp::max(exon.0, interval.0.min - 1);

                            let strandtype = if read.is_reverse() == (*interval.1 == b'-') {
                                's'
                            } else {
                                'r'
                            };
                            if read_number == 1 {
                                let at = autostrand_totals.get_mut(&strandtype).r()?;
                                *at += overlap_length as i64
                            } else if read_number == 2 {
                                let at2 = autostrand_totals2.get_mut(&strandtype).r()?;
                                *at2 += overlap_length as i64
                            }
                        }
                    }
                }
            } else {
                let tuple = (read_num, strand.to_string());
                if !histogram.contains_key(&tuple) {
                    histogram.insert(tuple.clone(), Vec::new());
                }
                // keep track of chromosome sizes
                if ref_length < exon.1 as u32 {
                    refs[read.tid() as usize].0 = exon.1 as u32;
                }
                if histogram[&tuple].len() < ref_length as usize {
                    let h = histogram.get_mut(&tuple).r()?;
                    h.resize(ref_length as usize, 0);
                }

                for pos in exon.0..exon.1 {
                    let h = histogram.get_mut(&tuple).r()?;
                    (*h)[pos as usize] += 1;
                }
            }
        }
    }

    if !autostrand_pass && !histogram.is_empty() && lastchr != -1 {
        write_chr(options,
                  &refs[lastchr as usize],
                  &histogram,
                  &mut fhs,
                  split_strand)?;
    }

    // make sure empty files were created
    if histogram.is_empty() && !autostrand_pass {
        for read_number in if options.split_read {
            vec![1, 2]
        } else {
            vec![0]
        } {
            for s in if split_strand != "uu" {
                vec!["+", "-"]
            } else {
                vec![""]
            } {
                open_file(options, read_number, s, split_strand, &mut fhs)?;
            }
        }
    }

    // close the filehandles
    for fh in &mut fhs {
        *fh.1 = None;
    }

    if autostrand_pass {
        // get the read 1 and read2 totals
        let mut total1: i64 = 0;
        let mut total2: i64 = 0;
        for i in &autostrand_totals {
            total1 += *i.1;
        }
        for i in &autostrand_totals2 {
            total2 += *i.1;
        }
        // figure out the best and second-best strand types for reads 1 and 2
        let mut best1: (char, i64) = ('\0', 0);
        let mut second_best1: (char, i64) = ('\0', 0);
        let mut best2: (char, i64) = ('\0', 0);
        let mut second_best2: (char, i64) = ('\0', 0);
        for i in autostrand_totals {
            if i.1 > 0 {
                writeln!(stderr(),
                         "Total evidence for read 1 strand type {}: {}",
                         i.0,
                         i.1)?;
            }
            if best1.1 < i.1 {
                second_best1 = best1;
                best1 = i;
            } else if second_best1.1 < i.1 {
                second_best1 = i;
            }
        }
        for i in autostrand_totals2 {
            if i.1 > 0 {
                writeln!(stderr(),
                         "Total evidence for read 2 strand type {}: {}",
                         i.0,
                         i.1)?;
            }
            if best2.1 < i.1 {
                second_best2 = best2;
                best2 = i;
            } else if second_best2.1 < i.1 {
                second_best2 = i;
            }
        }
        let threshold: f64 = 0.0; //threshold isn't working, set to zero
        let strand1 = if total1 > 0 {
            if threshold < (best1.1 - second_best1.1) as f64 / (total1) as f64 {
                best1.0
            } else {
                'u'
            }
        } else {
            'u'
        };
        let strand2 = if total2 > 0 {
            if threshold < (best2.1 - second_best2.1) as f64 / (total2) as f64 {
                best2.0
            } else {
                'u'
            }
        } else {
            'u'
        };

        let best_strand = format!("{}{}", strand1, strand2);
        writeln!(stderr(),
                 "autostrand_pass found best strand type: {}",
                 best_strand)?;

        // re-run analyzeBam with the strand type indicated
        analyze_bam(options, &best_strand, false, intervals)?;
    }

    if !autostrand_pass && options.bigwig {
        // Convert bedgraph file to bigwig file
        for fh in &fhs {
            // write the genome file for bigwigs
            let genome_filename = format!("{}.genome", fh.0);
            {
                let mut genome_fh = File::create(&genome_filename)?;
                for r in &refs {
                    writeln!(genome_fh, "{}\t{}", r.1, r.0)?;
                }
            }

            // run bedGraphToBigWig
            let regex = Regex::new(r"\.bedgraph$")?;
            let bigwig_file = String::from(regex.replace(fh.0, ".bw"));
            let sorted_bedgraph = String::from(regex.replace(fh.0, ".sorted.bedgraph"));
            for command in vec![Command::new("sort")
                                    .args(&["-k1,1", "-k2,2n", "-o", &sorted_bedgraph, fh.0])
                                    .env("LC_COLLATE", "C"),
                                Command::new("bedGraphToBigWig")
                                    .args(&[&sorted_bedgraph, &genome_filename, &bigwig_file])] {
                let mut child = command.spawn()?;
                let exit_code = child.wait()?;
                if !exit_code.success() {
                    if exit_code.code().is_some() {
                        let code = exit_code.code().r()?;
                        return Err(format!("Nonzero exit code {} returned from \
                                                       command: {:?}",
                                           code,
                                           command)
                            .into());
                    } else {
                        return Err(format!("Command was interrupted: {:?}", command).into());
                    }
                }
            }
            // remove the bedgraph file
            std::fs::remove_file(fh.0)?;
            // remove the sorted bedgraph file
            std::fs::remove_file(&sorted_bedgraph)?;
            // remove the genome file
            std::fs::remove_file(&genome_filename)?;
        }
    };
    Ok(())
}

#[derive(StructOpt, Debug)]
#[structopt(name = "bam2bedgraph", about = "Convert bam files to bedgraph/bigWig format")]
struct Options {
    #[structopt(long = "split", help = "Use CIGAR string to split alignment into separate exons (default)", default_value="true")]
    split_exons: bool,
    #[structopt(long = "read", help = "Split output bedgraph by read number")]
    split_read: bool,
    #[structopt(long = "zero", help = "Pad output bedgraph with zeroes")]
    zero: bool,
    #[structopt(long = "paired", help = "Only output paired read alignments")]
    paired_only: bool,
    #[structopt(long = "proper", help = "Only output proper-paired read alignments")]
    proper_only: bool,
    #[structopt(long = "primary", help = "Only output primary read alignments")]
    primary_only: bool,
    #[structopt(long = "trackline", help = "Output a UCSC track line (default)", default_value="true")]
    trackline: bool,
    #[structopt(long = "bigwig", help = "Output bigwig files (requires bedGraphToBigWig in $PATH)")]
    bigwig: bool,
    #[structopt(long = "uniq", help = "Keep only unique alignments (NH:i:1)")]
    uniq: bool,
    #[structopt(long = "fixchr", help = "Transform chromosome names to be UCSC-compatible")]
    fixchr: bool,
    #[structopt(help = "Convert a bam file into a bedgraph/bigwig file.", name="BAMFILE")]
    bamfile: String,
    #[structopt(long = "trackname", help = "Name of track for the track line", name="TRACKNAME", default_value="")]
    trackname: String,
    #[structopt(long = "out", help = "Output file prefix", name="PREFIX", default_value="")]
    out: String,
    #[structopt(long = "autostrand", help = 
        "Attempt to determine the strandedness of the input data using an \
        annotation file. Must be a .bam file.", name="AUTOSTRAND_FILE", default_value="")]
    autostrand: String,
    #[structopt(long = "split_strand", help =
        "Split output bedgraph by strand: Possible values: u s r uu us ur su ss \
        sr ru rs rr, first char is read1, second is read2, u=unstranded, \
        s=stranded, r=reverse", default_value = "uu", name="DESCRIPTION")]
    split_strand: String,
}

fn run() -> Result<()> {
    let mut options = Options::from_args();
    options.split_strand = options.split_strand.to_ascii_lowercase();
    if options.split_strand.len() == 1 {
        options.split_strand = options.split_strand + "u";
    }
    let regex = Regex::new(r"^[usr][usr]$")?;
    if !regex.is_match(&options.split_strand) {
        return Err(format!("Invalid value for split_strand: \"{}\": values must be \
                                       one of: u s r uu us ur su ss sr ru rs rr",
                           options.split_strand)
            .into());
    }

    // read in the annotation file
    let mut intervals: Option<BTreeMap<String, theban_interval_tree::IntervalTree<u8>>> = None;
    if !options.autostrand.is_empty() {
        if !Path::new(&options.autostrand).exists() {
            return Err(format!("Autostrand Bam file {} could not be found!",
                               &options.autostrand)
                .into());
        }
        let bam = rust_htslib::bam::Reader::from_path(&options.autostrand)?;
        let header = bam.header();

        let mut refs: Vec<(u32, String)> = Vec::new();
        refs.resize(header.target_count() as usize, (0, "".to_string()));
        let target_names = header.target_names();
        for target_name in target_names {
            let tid = header.tid(target_name).r()?;
            let target_len = header.target_len(tid).r()?;
            let target_name = std::str::from_utf8(target_name)?;
            refs[tid as usize] = (target_len, target_name.to_string());
        }
        if options.fixchr {
            for r in &mut refs {
                let regex = Regex::new(r"^(chr|Zv9_)")?;
                if regex.is_match(&r.1) {
                    let refname = r.1.to_string();
                    r.1.clear();
                    r.1.push_str(&format!("chr{}", refname));
                }
            }
        }

        let mut interval_lists: BTreeMap<String, Vec<(Range, u8)>> = BTreeMap::new();
        let mut read = rust_htslib::bam::record::Record::new();
        while bam.read(&mut read).is_ok() {
            let chr = refs[read.tid() as usize].1.clone();
            if !interval_lists.contains_key(&chr) {
                interval_lists.insert(chr.clone(), Vec::new());
            }

            let mut exons: Vec<(u64, u64)> = Vec::new();
            cigar2exons(&mut exons, &read.cigar(), read.pos() as u64)?;
            let interval_list = interval_lists.get_mut(&chr).r()?;
            for exon in exons {
                interval_list.push((Range::new(exon.0 + 1, exon.1),
                                    if read.is_reverse() { b'-' } else { b'+' }));
            }
        }
        for (chr, list) in &interval_lists {
            if intervals.is_none() {
                intervals = Some(BTreeMap::new());
            }
            let interval = intervals.as_mut().r()?;
            let mut tree = theban_interval_tree::IntervalTree::<u8>::new();
            for l in list {
                tree.insert(l.0, l.1);
            }
            interval.insert(chr.clone(), tree);
        }
    }

    // // analyze the bam file and produce histograms
    if !options.autostrand.is_empty() {
        // make both stranded and unstranded files
        analyze_bam(&options,
                    &options.split_strand,
                    !options.autostrand.is_empty(),
                    &intervals)?;
        analyze_bam(&options, "uu", false, &intervals)?;
    } else {
        analyze_bam(&options,
                    &options.split_strand,
                    !options.autostrand.is_empty(),
                    &intervals)?;
    }
    Ok(())
}

fn main() {
    // enable stack traces
    std::env::set_var("RUST_BACKTRACE", "1");

    if let Err(ref e) = run() {
        println!("error: {}", e);
        for e in e.iter().skip(1) {
            println!("caused by: {}", e);
        }
        if let Some(backtrace) = e.backtrace() {
            println!("backtrace: {:?}", backtrace);
        }
        ::std::process::exit(1);
    }
}