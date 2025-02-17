use anyhow::Result;
use clap::{arg, command, value_parser};
use gffplot::parse::parse_gff;
use std::path::PathBuf;

fn main() -> Result<()> {
    let matches = command!()
        .help_template("gffplot {version}\nUsage: gffplot in.gff > out.html")
        .arg_required_else_help(true)
        .arg(arg!(<GFF> "The input GFF file").value_parser(value_parser!(PathBuf)))
        .get_matches();

    let gff = matches.get_one::<PathBuf>("GFF").unwrap();

    let parsed_gff = parse_gff(gff.clone())?;

    parsed_gff.plot()?;

    Ok(())
}
