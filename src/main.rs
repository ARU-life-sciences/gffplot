use anyhow::Result;
use clap::{arg, command, value_parser};
use gffplot::{parse::parse_gff, plot::PlotData};
use std::path::PathBuf;

fn main() -> Result<()> {
    let matches = command!()
        .arg(arg!(<GFF> "The input GFF file").value_parser(value_parser!(PathBuf)))
        .get_matches();

    let gff = matches.get_one::<PathBuf>("GFF").unwrap();

    let parsed_gff = parse_gff(gff.clone())?;

    parsed_gff.plot(PathBuf::from("./test.html"))?;

    Ok(())
}
