//! OxoG Context Filter
//!
//! Filters a GVCF/VCF for OxoG (8-oxoguanine) artifacts with sequence-context
//! awareness (CpG islands, homopolymer runs).
//!
//! OxoG artifacts manifest as:
//! - C>A on the forward strand (REF=C, ALT=A)
//! - G>T on the reverse strand (REF=G, ALT=T)
//! with characteristic strand bias (high SOR, negative ReadPosRankSum).
//!
//! Sequence-context modules:
//! 1. CpG islands  — Stricter filtering; OxoG C>A in a CpG context is almost
//!                   certainly artifactual (real CpG changes are C>T transitions).
//! 2. Homopolymers — Variants in G-runs / C-runs are flagged as likely artifacts
//!                   because polymerase slippage + oxidative damage compound.

use anyhow::{bail, Context, Result};
use clap::Parser;
use log::info;
use rust_htslib::bcf::{header::HeaderRecord, Format, Read, Reader, Writer};
use rust_htslib::faidx;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::PathBuf;

// ---------------------------------------------------------------------------
// CLI Arguments
// ---------------------------------------------------------------------------
#[derive(Parser, Debug)]
#[command(
    name = "oxoq_filt",
    about = "OxoG artifact filter with CpG island and homopolymer context awareness",
    long_about = r#"
Filters a GVCF/VCF for OxoG (8-oxoguanine) artifacts with sequence-context
awareness (CpG islands, homopolymer runs).

OxoG artifacts manifest as:
  - C>A on the forward strand (REF=C, ALT=A)
  - G>T on the reverse strand (REF=G, ALT=T)
  with characteristic strand bias (high SOR, negative ReadPosRankSum).

Examples:
  # Soft filter (annotate FILTER column, keep all records):
  oxoq_filt -i input.g.vcf.gz -o filtered.g.vcf.gz -r ref.fa

  # Hard filter (remove OxoG artifacts from output):
  oxoq_filt -i input.g.vcf.gz -o filtered.g.vcf.gz -r ref.fa --hard-filter

  # With CpG island annotations + custom thresholds:
  oxoq_filt -i input.g.vcf.gz -o filtered.g.vcf.gz -r ref.fa \
      --cpg-bed cpgIslandExt.bed --sor-threshold 3.0 --cpg-sor-threshold 1.5
"#
)]
struct Args {
    /// Input GVCF/VCF (plain or .gz)
    #[arg(short, long)]
    input: PathBuf,

    /// Output VCF (.vcf or .vcf.gz)
    #[arg(short, long)]
    output: PathBuf,

    /// Reference FASTA (must have .fai index)
    #[arg(short, long)]
    reference: PathBuf,

    /// BED file of CpG islands (e.g., UCSC cpgIslandExt)
    #[arg(long)]
    cpg_bed: Option<PathBuf>,

    /// Remove filtered variants from output (default: soft-filter / annotate only)
    #[arg(long, default_value = "false")]
    hard_filter: bool,

    /// SOR threshold for standard OxoG filter
    #[arg(long, default_value = "2.5")]
    sor_threshold: f64,

    /// ReadPosRankSum threshold for standard OxoG filter
    #[arg(long, default_value = "-2.0")]
    rprs_threshold: f64,

    /// SOR threshold in CpG context (stricter)
    #[arg(long, default_value = "2.0")]
    cpg_sor_threshold: f64,

    /// ReadPosRankSum threshold in CpG context (stricter)
    #[arg(long, default_value = "-1.0")]
    cpg_rprs_threshold: f64,

    /// SOR threshold in homopolymer context
    #[arg(long, default_value = "2.0")]
    homo_sor_threshold: f64,

    /// ReadPosRankSum threshold in homopolymer context
    #[arg(long, default_value = "-1.0")]
    homo_rprs_threshold: f64,

    /// SOR-only threshold when ReadPosRankSum is missing (stricter)
    #[arg(long, default_value = "4.0")]
    sor_only_threshold: f64,

    /// SOR-only threshold in CpG context when ReadPosRankSum is missing
    #[arg(long, default_value = "3.0")]
    cpg_sor_only_threshold: f64,

    /// SOR-only threshold in homopolymer context when ReadPosRankSum is missing
    #[arg(long, default_value = "3.0")]
    homo_sor_only_threshold: f64,

    /// Minimum run length to classify as homopolymer
    #[arg(long, default_value = "4")]
    homopolymer_min_length: usize,

    /// Flanking bp to inspect for homopolymer runs
    #[arg(long, default_value = "5")]
    homopolymer_flank: usize,

    /// Write filtering summary to this file
    #[arg(long)]
    summary_out: Option<PathBuf>,

    /// Verbose logging
    #[arg(short, long, default_value = "false")]
    verbose: bool,
}

// ---------------------------------------------------------------------------
// Configuration
// ---------------------------------------------------------------------------
#[derive(Debug, Clone)]
struct FilterConfig {
    sor_threshold: f64,
    rprs_threshold: f64,
    cpg_sor_threshold: f64,
    cpg_rprs_threshold: f64,
    homo_sor_threshold: f64,
    homo_rprs_threshold: f64,
    sor_only_threshold: f64,
    cpg_sor_only_threshold: f64,
    homo_sor_only_threshold: f64,
    homopolymer_min_length: usize,
    homopolymer_flank: usize,
    hard_filter: bool,
}

impl From<&Args> for FilterConfig {
    fn from(args: &Args) -> Self {
        FilterConfig {
            sor_threshold: args.sor_threshold,
            rprs_threshold: args.rprs_threshold,
            cpg_sor_threshold: args.cpg_sor_threshold,
            cpg_rprs_threshold: args.cpg_rprs_threshold,
            homo_sor_threshold: args.homo_sor_threshold,
            homo_rprs_threshold: args.homo_rprs_threshold,
            sor_only_threshold: args.sor_only_threshold,
            cpg_sor_only_threshold: args.cpg_sor_only_threshold,
            homo_sor_only_threshold: args.homo_sor_only_threshold,
            homopolymer_min_length: args.homopolymer_min_length,
            homopolymer_flank: args.homopolymer_flank,
            hard_filter: args.hard_filter,
        }
    }
}

// ---------------------------------------------------------------------------
// CpG Island Index
// ---------------------------------------------------------------------------
/// Fast interval lookup for CpG islands.
/// Loads a BED file into a hashmap of sorted interval vectors, then does
/// binary search for overlap queries.
struct CpGIndex {
    intervals: HashMap<String, Vec<(u64, u64)>>,
}

impl CpGIndex {
    fn new(bed_path: Option<&PathBuf>) -> Result<Self> {
        let mut idx = CpGIndex {
            intervals: HashMap::new(),
        };
        if let Some(path) = bed_path {
            idx.load(path)?;
        }
        Ok(idx)
    }

    fn load(&mut self, bed_path: &PathBuf) -> Result<()> {
        info!("Loading CpG island BED: {:?}", bed_path);
        let file = File::open(bed_path).context("Failed to open CpG BED file")?;
        let reader = BufReader::new(file);
        let mut count = 0;

        for line in reader.lines() {
            let line = line?;
            if line.starts_with('#') || line.starts_with("track") || line.starts_with("browser") {
                continue;
            }
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() < 3 {
                continue;
            }
            let chrom = parts[0].to_string();
            let start: u64 = parts[1].parse().context("Invalid start position in BED")?;
            let end: u64 = parts[2].parse().context("Invalid end position in BED")?;  
            self.intervals.entry(chrom).or_default().push((start, end));
            count += 1;
        }

        // Sort each chromosome's intervals for binary search
        for ivs in self.intervals.values_mut() {
            ivs.sort_by_key(|&(s, _)| s);
        }

        info!(
            "Loaded {} CpG island intervals across {} contigs",
            count,
            self.intervals.len()
        );
        Ok(())
    }

    /// Check if a 0-based position overlaps any CpG island.
    fn overlaps(&self, chrom: &str, pos: u64) -> bool {
        let ivs = match self.intervals.get(chrom) {
            Some(v) => v,
            None => return false,
        };

        // Binary search for first interval whose end > pos
        let idx = ivs.partition_point(|&(_, end)| end <= pos);

        // Check if 'idx' interval contains pos
        if idx < ivs.len() && ivs[idx].0 <= pos && pos < ivs[idx].1 {
            return true;
        }
        false
    }
}

// ---------------------------------------------------------------------------
// Homopolymer Detector
// ---------------------------------------------------------------------------
struct HomopolymerDetector {
    fasta: faidx::Reader,
    min_length: usize,
    flank: usize,
}

impl HomopolymerDetector {
    fn new(fasta_path: &PathBuf, min_length: usize, flank: usize) -> Result<Self> {
        let fasta = faidx::Reader::from_path(fasta_path)
            .context("Failed to open reference FASTA")?;
        Ok(HomopolymerDetector {
            fasta,
            min_length,
            flank,
        })
    }

    /// Check if pos (0-based) sits inside or adjacent to a homopolymer run.
    /// Returns (is_homopolymer, run_length, run_base)
    fn is_in_homopolymer(&self, chrom: &str, pos: u64) -> (bool, usize, char) {
        let chrom_len = self.fasta.fetch_seq_len(chrom);
        if chrom_len == 0 {
            return (false, 0, ' ');
        }

        let start = pos.saturating_sub(self.flank as u64);
        let end = (pos + self.flank as u64 + 1).min(chrom_len);

        let seq = match self.fasta.fetch_seq_string(chrom, start as usize, end as usize - 1) {
            Ok(s) => s.to_uppercase(),
            Err(_) => return (false, 0, ' '),
        };

        if seq.is_empty() {
            return (false, 0, ' ');
        }

        let seq_bytes = seq.as_bytes();
        let var_idx = (pos - start) as usize;

        let mut best_len = 0;
        let mut best_base = ' ';
        let mut i = 0;

        while i < seq_bytes.len() {
            let base = seq_bytes[i] as char;
            let run_start = i;
            while i < seq_bytes.len() && seq_bytes[i] as char == base {
                i += 1;
            }
            let run_end = i;
            let run_len = run_end - run_start;

            // Check if this run covers the variant position
            if run_start <= var_idx && var_idx < run_end && run_len >= self.min_length {
                if run_len > best_len {
                    best_len = run_len;
                    best_base = base;
                }
            }
        }

        if best_len >= self.min_length {
            (true, best_len, best_base)
        } else {
            (false, 0, ' ')
        }
    }
}

// ---------------------------------------------------------------------------
// Reference FASTA helper
// ---------------------------------------------------------------------------
struct ReferenceHelper {
    fasta: faidx::Reader,
}

impl ReferenceHelper {
    fn new(fasta_path: &PathBuf) -> Result<Self> {
        let fasta = faidx::Reader::from_path(fasta_path)
            .context("Failed to open reference FASTA")?;
        Ok(ReferenceHelper { fasta })
    }

    /// Check if the trinucleotide context (upstream base + ref + downstream base)
    /// has GC content > 0.5 (i.e., at least 2 of 3 bases are G or C).
    /// pos is 0-based.
    fn is_cpg_dinucleotide(&self, chrom: &str, pos: u64, ref_base: char) -> bool {
        let chrom_len = self.fasta.fetch_seq_len(chrom);
        if chrom_len == 0 || pos == 0 || pos + 1 >= chrom_len {
            return false;
        }

        let trinuc = match self.fasta.fetch_seq_string(chrom, (pos - 1) as usize, (pos + 1) as usize) {
            Ok(s) => s.to_uppercase(),
            Err(_) => return false,
        };

        if trinuc.len() != 3 {
            return false;
        }

        let gc_count = trinuc.chars().filter(|&c| c == 'G' || c == 'C').count();
        gc_count > 1
    }
}

// ---------------------------------------------------------------------------
// Core OxoG filter logic
// ---------------------------------------------------------------------------

/// True if this is an OxoG-signature substitution: C>A or G>T.
fn is_oxog_candidate(ref_base: char, alt_base: char) -> bool {
    (ref_base == 'C' && alt_base == 'A') || (ref_base == 'G' && alt_base == 'T')
}

/// Convert a strand bias (SB) table to a Strand Odds Ratio (SOR) value.
/// SB is a 4-element array: [refFwd, refRev, altFwd, altRev]
fn sb_to_sor(sb: &[i32]) -> Option<f64> {
    if sb.len() < 4 {
        return None;
    }

    let ref_fwd = sb[0] as f64;
    let ref_rev = sb[1] as f64;
    let alt_fwd = sb[2] as f64;
    let alt_rev = sb[3] as f64;

    let pseudocount = 1.0;

    // Symmetric ratio (cross-product ratio)
    let ratio = (ref_fwd * alt_rev + pseudocount) / (ref_rev * alt_fwd + pseudocount);

    // Ref strand balance: closer to 1 means balanced
    let ref_ratio = (ref_fwd.min(ref_rev) + pseudocount) / (ref_fwd.max(ref_rev) + pseudocount);

    // Alt strand balance
    let alt_ratio = (alt_fwd.min(alt_rev) + pseudocount) / (alt_fwd.max(alt_rev) + pseudocount);

    let sor = ratio.ln() + ref_ratio.ln() - alt_ratio.ln();
    Some(sor)
}

// ---------------------------------------------------------------------------
// Filter Result
// ---------------------------------------------------------------------------
#[derive(Debug, Default)]
struct FilterResult {
    filters_applied: Vec<String>,
    annotations: HashMap<String, String>,
}

impl FilterResult {
    fn is_filtered(&self) -> bool {
        !self.filters_applied.is_empty()
    }
}

// ---------------------------------------------------------------------------
// Filter Statistics
// ---------------------------------------------------------------------------
#[derive(Debug, Default)]
struct FilterStats {
    total_variants: u64,
    total_snps: u64,
    oxog_candidates: u64,
    filtered_standard: u64,
    filtered_sor_only: u64,
    filtered_cpg: u64,
    filtered_homopolymer: u64,
    filtered_homopoly_strandbias: u64,
    hard_removed: u64,
    passed: u64,
    filter_counts: HashMap<String, u64>,
}

impl FilterStats {
    fn summary(&self) -> String {
        let mut lines = vec![
            "=".repeat(60),
            "OxoG + Context Filter Summary".to_string(),
            "=".repeat(60),
            format!("Total variant records processed:   {:>10}", self.total_variants),
            format!("  SNPs examined:                   {:>10}", self.total_snps),
            format!("  OxoG candidates (C>A / G>T):     {:>10}", self.oxog_candidates),
            "-".repeat(60),
            format!("Filtered (OxoG standard):          {:>10}", self.filtered_standard),
            format!("Filtered (OxoG SOR-only):           {:>10}", self.filtered_sor_only),
            format!("Filtered (OxoG in CpG context):    {:>10}", self.filtered_cpg),
            format!("Filtered (OxoG in homopolymer):    {:>10}", self.filtered_homopolymer),
            format!("Flagged  (HomoPoly strand bias):   {:>10}", self.filtered_homopoly_strandbias),
            "-".repeat(60),
        ];

        if self.hard_removed > 0 {
            lines.push(format!("Hard-removed from output:          {:>10}", self.hard_removed));
        }

        lines.push(format!("PASS (after filtering):            {:>10}", self.passed));
        lines.push("=".repeat(60));
        lines.push(String::new());
        lines.push("Per-filter breakdown:".to_string());

        let mut sorted_counts: Vec<_> = self.filter_counts.iter().collect();
        sorted_counts.sort_by(|a, b| b.1.cmp(a.1));
        for (fname, count) in sorted_counts {
            lines.push(format!("  {:<35} {:>8}", fname, count));
        }

        lines.join("\n")
    }
}

// ---------------------------------------------------------------------------
// Main processing
// ---------------------------------------------------------------------------
fn get_info_float(record: &rust_htslib::bcf::Record, key: &[u8]) -> Option<f64> {
    match record.info(key).float() {
        Ok(Some(vals)) => vals.first().copied().map(|v| v as f64),
        _ => None,
    }
}

fn get_info_int_array(record: &rust_htslib::bcf::Record, key: &[u8]) -> Option<Vec<i32>> {
    match record.info(key).integer() {
        Ok(Some(vals)) => Some(vals.to_vec()),
        _ => None,
    }
}

fn get_format_int_array(record: &rust_htslib::bcf::Record, key: &[u8], sample_idx: usize) -> Option<Vec<i32>> {
    match record.format(key).integer() {
        Ok(vals) => {
            if sample_idx < vals.len() {
                Some(vals[sample_idx].to_vec())
            } else {
                None
            }
        }
        _ => None,
    }
}

fn get_sor_value(record: &rust_htslib::bcf::Record) -> (Option<f64>, &'static str) {
    // Try INFO/SOR first
    if let Some(sor) = get_info_float(record, b"SOR") {
        return (Some(sor), "SOR");
    }

    // Try FORMAT/SB (per-sample)
    if let Some(sb) = get_format_int_array(record, b"SB", 0) {
        if let Some(sor) = sb_to_sor(&sb) {
            return (Some(sor), "SB_format");
        }
    }

    // Try INFO/SB
    if let Some(sb) = get_info_int_array(record, b"SB") {
        if let Some(sor) = sb_to_sor(&sb) {
            return (Some(sor), "SB_info");
        }
    }

    (None, "missing")
}

/// Compute SOR from SB if SB data is available, regardless of other fields.
fn compute_sor_from_sb(record: &rust_htslib::bcf::Record) -> Option<f64> {
    if let Some(sb) = get_format_int_array(record, b"SB", 0) {
        if let Some(sor) = sb_to_sor(&sb) {
            return Some(sor);
        }
    }
    if let Some(sb) = get_info_int_array(record, b"SB") {
        if let Some(sor) = sb_to_sor(&sb) {
            return Some(sor);
        }
    }
    None
}

fn evaluate_variant(
    record: &rust_htslib::bcf::Record,
    config: &FilterConfig,
    cpg_index: &CpGIndex,
    homo_detector: &HomopolymerDetector,
    ref_helper: &ReferenceHelper,
    header: &rust_htslib::bcf::header::HeaderView,
) -> FilterResult {
    let mut result = FilterResult::default();

    let alleles = record.alleles();
    // In GVCF, variants carry an extra <NON_REF> allele — filter it out
    let real_alleles: Vec<&[u8]> = alleles.iter().filter(|a| *a != b"<NON_REF>").copied().collect();
    if real_alleles.len() != 2 {
        return result;
    }

    let ref_allele = real_alleles[0];
    let alt_allele = real_alleles[1];

    // Only process biallelic SNPs
    if ref_allele.len() != 1 || alt_allele.len() != 1 {
        return result;
    }

    let ref_base = ref_allele[0].to_ascii_uppercase() as char;
    let alt_base = alt_allele[0].to_ascii_uppercase() as char;

    if !is_oxog_candidate(ref_base, alt_base) {
        return result;
    }

    // Retrieve metrics
    let (sor, _sor_source) = get_sor_value(record);
    let rprs = get_info_float(record, b"ReadPosRankSum");


    // If SOR is missing entirely, we can't evaluate
    if sor.is_none() {
        result.annotations.insert("OxoG_Note".to_string(), "Missing_SOR/SB".to_string());
        return result;
    }

    let sor = sor.unwrap();

    // Get chromosome and position (0-based)
    let chrom = String::from_utf8_lossy(header.rid2name(record.rid().unwrap()).unwrap()).to_string();
    let pos = record.pos() as u64;

    // --- CpG context ---
    let in_cpg_island = cpg_index.overlaps(&chrom, pos);
    let in_cpg_dinuc = ref_helper.is_cpg_dinucleotide(&chrom, pos, ref_base);
    let in_cpg_context = in_cpg_island || in_cpg_dinuc;

    // --- Homopolymer context ---
    let (in_homo, homo_len, homo_base) = homo_detector.is_in_homopolymer(&chrom, pos);

    // Annotate context
    if in_cpg_context {
        let mut parts = Vec::new();
        if in_cpg_island {
            parts.push("CpG_island");
        }
        if in_cpg_dinuc {
            parts.push("CpG_dinuc");
        }
        result.annotations.insert("CpG_Context".to_string(), parts.join("|"));
    }

    if in_homo {
        result.annotations.insert("Homopolymer".to_string(), format!("{}x{}", homo_base, homo_len));
    }

    // --- Apply filters with context-dependent thresholds ---
    if let Some(rprs) = rprs {
        // Both SOR and RPRS available — use standard dual thresholds
        let (sor_thresh, rprs_thresh, filter_prefix) = if in_cpg_context {
            (config.cpg_sor_threshold, config.cpg_rprs_threshold, "OxoG_CpG")
        } else if in_homo {
            (config.homo_sor_threshold, config.homo_rprs_threshold, "OxoG_HomoPoly")
        } else {
            (config.sor_threshold, config.rprs_threshold, "OxoG")
        };

        if sor > sor_thresh && rprs < rprs_thresh {
            let direction = if ref_base == 'C' && alt_base == 'A' { "fwd" } else { "rev" };
            let filter_name = format!("{}_{}", filter_prefix, direction);
            result.filters_applied.push(filter_name);
        }
    } else {
        // RPRS missing — fall back to SOR-only with stricter threshold
        let (sor_thresh, filter_prefix) = if in_cpg_context {
            (config.cpg_sor_only_threshold, "OxoG_CpG_SORonly")
        } else if in_homo {
            (config.homo_sor_only_threshold, "OxoG_HomoPoly_SORonly")
        } else {
            (config.sor_only_threshold, "OxoG_SORonly")
        };

        if sor > sor_thresh {
            let direction = if ref_base == 'C' && alt_base == 'A' { "fwd" } else { "rev" };
            let filter_name = format!("{}_{}", filter_prefix, direction);
            result.filters_applied.push(filter_name);
        }
    }

    // Additionally flag homopolymer variants even if OxoG thresholds aren't met
    if in_homo && !result.is_filtered() && sor > config.sor_threshold {
        result.filters_applied.push("HomoPoly_StrandBias".to_string());
    }

    result
}

fn process_vcf(args: &Args, config: &FilterConfig) -> Result<()> {
    info!("Opening input: {:?}", args.input);
    let mut vcf_in = Reader::from_path(&args.input).context("Failed to open input VCF")?;

    // Prepare reference
    info!("Opening reference: {:?}", args.reference);
    let ref_helper = ReferenceHelper::new(&args.reference)?;
    let homo_detector = HomopolymerDetector::new(
        &args.reference,
        config.homopolymer_min_length,
        config.homopolymer_flank,
    )?;

    // Prepare CpG index
    let cpg_index = CpGIndex::new(args.cpg_bed.as_ref())?;

    // Set up output header with new filters and annotations
    let mut header = rust_htslib::bcf::Header::from_template(vcf_in.header());

    // Add filter headers
    let filter_descriptions = [
        ("OxoG_fwd", format!("OxoG artifact: C>A with SOR>{} and ReadPosRankSum<{}", config.sor_threshold, config.rprs_threshold)),
        ("OxoG_rev", format!("OxoG artifact: G>T with SOR>{} and ReadPosRankSum<{}", config.sor_threshold, config.rprs_threshold)),
        ("OxoG_CpG_fwd", format!("OxoG artifact in CpG context: C>A with SOR>{} and ReadPosRankSum<{}", config.cpg_sor_threshold, config.cpg_rprs_threshold)),
        ("OxoG_CpG_rev", format!("OxoG artifact in CpG context: G>T with SOR>{} and ReadPosRankSum<{}", config.cpg_sor_threshold, config.cpg_rprs_threshold)),
        ("OxoG_HomoPoly_fwd", format!("OxoG artifact in homopolymer: C>A with SOR>{} and ReadPosRankSum<{}", config.homo_sor_threshold, config.homo_rprs_threshold)),
        ("OxoG_HomoPoly_rev", format!("OxoG artifact in homopolymer: G>T with SOR>{} and ReadPosRankSum<{}", config.homo_sor_threshold, config.homo_rprs_threshold)),
        ("HomoPoly_StrandBias", format!("Variant in homopolymer run with elevated strand bias (SOR>{})", config.sor_threshold)),
        ("OxoG_SORonly_fwd", format!("OxoG artifact: C>A with SOR>{} (ReadPosRankSum missing)", config.sor_only_threshold)),
        ("OxoG_SORonly_rev", format!("OxoG artifact: G>T with SOR>{} (ReadPosRankSum missing)", config.sor_only_threshold)),
        ("OxoG_CpG_SORonly_fwd", format!("OxoG artifact in CpG context: C>A with SOR>{} (ReadPosRankSum missing)", config.cpg_sor_only_threshold)),
        ("OxoG_CpG_SORonly_rev", format!("OxoG artifact in CpG context: G>T with SOR>{} (ReadPosRankSum missing)", config.cpg_sor_only_threshold)),
        ("OxoG_HomoPoly_SORonly_fwd", format!("OxoG artifact in homopolymer: C>A with SOR>{} (ReadPosRankSum missing)", config.homo_sor_only_threshold)),
        ("OxoG_HomoPoly_SORonly_rev", format!("OxoG artifact in homopolymer: G>T with SOR>{} (ReadPosRankSum missing)", config.homo_sor_only_threshold)),
    ];

    for (filt_id, desc) in &filter_descriptions {
        header.push_record(format!("##FILTER=<ID={},Description=\"{}\">", filt_id, desc).as_bytes());
    }

    // Add INFO annotations
    let info_annotations = [
        ("CpG_Context", "String", "CpG context: CpG_island, CpG_dinuc, or both"),
        ("Homopolymer", "String", "Homopolymer run at variant (e.g., Gx6)"),
        ("OxoG_Note", "String", "Note about OxoG evaluation"),
    ];

    for (info_id, typ, desc) in &info_annotations {
        header.push_record(format!("##INFO=<ID={},Number=1,Type={},Description=\"{}\">", info_id, typ, desc).as_bytes());
    }

    // Ensure INFO/SOR header exists
    let has_sor = vcf_in.header().header_records().iter().any(|r| {
        matches!(r, HeaderRecord::Info { values, .. } if values.get("ID") == Some(&"SOR".to_string()))
    });
    if !has_sor {
        header.push_record(b"##INFO=<ID=SOR,Number=1,Type=Float,Description=\"Strand Odds Ratio (computed from FORMAT/SB when not present in input)\">");
    }

    // Open output
    info!("Writing output: {:?}", args.output);
    let output_path_str = args.output.to_string_lossy();
    let is_compressed = output_path_str.ends_with(".gz") || output_path_str.ends_with(".bgz");
    
    let mut vcf_out = if is_compressed {
        Writer::from_path(&args.output, &header, false, Format::Vcf)?
    } else {
        Writer::from_path(&args.output, &header, true, Format::Vcf)?
    };

    let mut stats = FilterStats::default();
    let header_view = vcf_in.header().clone();

    for result in vcf_in.records() {
        let mut record = result.context("Failed to read VCF record")?;
        stats.total_variants += 1;

        if stats.total_variants % 100_000 == 0 {
            info!("Processed {} records...", stats.total_variants);
        }

        // Check if it's a SNP (ignore <NON_REF> alleles from GVCF)
        let alleles = record.alleles();
        let real_alleles: Vec<&[u8]> = alleles.iter().filter(|a| *a != b"<NON_REF>").copied().collect();
        let is_snp = real_alleles.len() == 2
            && real_alleles[0].len() == 1
            && real_alleles[1].len() == 1
            && real_alleles[1] != b"*";

        if is_snp {
            stats.total_snps += 1;

            let ref_base = real_alleles[0][0].to_ascii_uppercase() as char;
            let alt_base = real_alleles[1][0].to_ascii_uppercase() as char;
            if is_oxog_candidate(ref_base, alt_base) {
                stats.oxog_candidates += 1;
            }
        }

        // Evaluate filters
        let fr = evaluate_variant(&record, config, &cpg_index, &homo_detector, &ref_helper, &header_view);

        // Update stats
        for f in &fr.filters_applied {
            *stats.filter_counts.entry(f.clone()).or_insert(0) += 1;
            if f.contains("SORonly") {
                stats.filtered_sor_only += 1;
            } else if f.contains("CpG") {
                stats.filtered_cpg += 1;
            } else if f.contains("HomoPoly_StrandBias") {
                stats.filtered_homopoly_strandbias += 1;
            } else if f.contains("HomoPoly") {
                stats.filtered_homopolymer += 1;
            } else {
                stats.filtered_standard += 1;
            }
        }

        // Hard filter: skip writing this record entirely
        if config.hard_filter && fr.is_filtered() {
            stats.hard_removed += 1;
            continue;
        }

        // Translate record to the output header (preserves all fields including FORMAT/sample data)
        vcf_out.translate(&mut record);

        // Add new annotations
        for (key, val) in &fr.annotations {
            record.push_info_string(key.as_bytes(), &[val.as_bytes()])?;
        }

        // Always compute and write SOR from SB when SB data is available
        if let Some(sor) = compute_sor_from_sb(&record) {
            record.push_info_float(b"SOR", &[sor as f32])?;
        }

        // Set filters
        if !fr.filters_applied.is_empty() {
            for f in &fr.filters_applied {
                record.push_filter(f.as_bytes())?;
            }
        } else {
            stats.passed += 1;
        }

        vcf_out.write(&record)?;
    }

    // Print summary
    let summary_text = stats.summary();
    println!("{}", summary_text);

    if let Some(summary_path) = &args.summary_out {
        let mut file = File::create(summary_path)?;
        writeln!(file, "{}", summary_text)?;
        info!("Summary written to: {:?}", summary_path);
    }

    info!("Done.");
    Ok(())
}

fn main() -> Result<()> {
    let args = Args::parse();

    // Initialize logging
    if args.verbose {
        env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("debug")).init();
    } else {
        env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();
    }

    // Validate inputs
    if !args.input.exists() {
        bail!("Input file not found: {:?}", args.input);
    }
    if !args.reference.exists() {
        bail!("Reference FASTA not found: {:?}", args.reference);
    }
    let fai_path = args.reference.with_extension("fa.fai");
    let fai_path2 = PathBuf::from(format!("{}.fai", args.reference.display()));
    if !fai_path.exists() && !fai_path2.exists() {
        bail!(
            "Reference index not found: {:?}\nRun: samtools faidx {:?}",
            fai_path2,
            args.reference
        );
    }
    if let Some(cpg_bed) = &args.cpg_bed {
        if !cpg_bed.exists() {
            bail!("CpG BED file not found: {:?}", cpg_bed);
        }
    }

    let config = FilterConfig::from(&args);

    info!("Filter configuration:");
    info!(
        "  Standard OxoG:   SOR > {:.1}  AND  ReadPosRankSum < {:.1}",
        config.sor_threshold, config.rprs_threshold
    );
    info!(
        "  CpG context:     SOR > {:.1}  AND  ReadPosRankSum < {:.1}  (stricter)",
        config.cpg_sor_threshold, config.cpg_rprs_threshold
    );
    info!(
        "  Homopolymer:     SOR > {:.1}  AND  ReadPosRankSum < {:.1}  (stricter)",
        config.homo_sor_threshold, config.homo_rprs_threshold
    );
    info!("  Homopolymer min length: {} bp", config.homopolymer_min_length);
    info!(
        "  SOR-only fallback: standard>{:.1}  CpG>{:.1}  HomoPoly>{:.1}",
        config.sor_only_threshold, config.cpg_sor_only_threshold, config.homo_sor_only_threshold
    );
    info!(
        "  Mode: {}",
        if config.hard_filter {
            "HARD filter (remove)"
        } else {
            "SOFT filter (annotate)"
        }
    );

    process_vcf(&args, &config)?;

    Ok(())
}
