use crate::plot::PlotData;
use anyhow::Result;
use bio::{
    bio_types::strand::Strand,
    io::gff::{GffType, Reader},
};
use std::{collections::BTreeMap, path::PathBuf};

const SOURCES: [&str; 4] = ["ORFfinder", "cmscan", "tRNAscan-SE", "barrnap:0.9"];

pub struct Row {
    pub feature_name: String,
    pub source: String,
    pub feature_type: String,
    pub start: u64,
    pub end: u64,
    pub strand: Strand,
}

pub fn parse_gff(file: PathBuf) -> Result<PlotData> {
    let mut reader = Reader::from_file(file, GffType::GFF3)?;

    let mut map = BTreeMap::new();
    for record in reader.records() {
        let r = record?;
        let seqname = r.seqname();
        let source = r.source();
        let feature_type = r.feature_type();
        // skip trnascan-se's exons
        if feature_type == "exon" {
            continue;
        }
        let start = r.start();
        let end = r.end();
        let strand = r.strand();
        // don't bother with phase
        let attributes = r.attributes();

        let feature_name = match source {
            "ORFfinder" => {
                let pfam_accession = attributes.get("pfam_accession").unwrap();
                let target_name = attributes.get("target_name").unwrap();
                let description = attributes.get("description").unwrap();

                format!("Pfam accession: {pfam_accession}\nTarget name: {target_name}\nDescription: {description}")
            }
            "cmscan" => {
                let description = attributes.get("description").unwrap();
                format!("Description: {description}")
            }
            "tRNAscan-SE" => {
                let gene_biotype = attributes.get("gene_biotype").unwrap();
                let anticodon = attributes.get("anticodon").unwrap();

                format!("{gene_biotype}: {anticodon}")
            }
            "barrnap:0.9" => {
                let product = attributes.get("product").unwrap();
                let note = attributes.get("note").unwrap();

                format!("{product}: {note}")
            }
            u => panic!("{u} not found"),
        };

        let row = Row {
            feature_name,
            source: source.to_string(),
            feature_type: feature_type.to_string(),
            start: *start,
            end: *end,
            strand: strand.unwrap(),
        };

        // add to BTreeMap
        let entry = map.entry(seqname.to_string()).or_insert(Vec::new());
        entry.push(row);
    }

    Ok(PlotData { data: map })
}
