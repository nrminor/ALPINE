use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader, LineWriter, Write};

fn main() {

    let args: Vec<String> = std::env::args().collect();

    if args.len() != 3 {
        eprintln!("Usage: {} <input.fasta> <accessions.txt>", args[0]);
        std::process::exit(1);
    }

    let fasta_path = &args[1];
    let accessions_path = &args[2];

    // Load accession numbers from the text file into a HashSet
    let accessions_of_interest = load_accessions(accessions_path).expect("Failed to load accessions. Please double check that a text file called accessions.txt is in the current working directory.");

    // Open the input FASTA file
    let file = File::open(fasta_path).expect("Failed to open input file");
    let mut reader = BufReader::new(file);

    // Open the output FASTA file for writing
    let output_file = File::create("filtered-to-geography.fasta").expect("Failed to create output file");
    let mut writer = LineWriter::new(output_file);

    let mut current_accession = String::new();
    let mut current_sequence = String::new();
    let mut is_sequence_of_interest = false;

    let mut line = String::new();
    while reader.read_line(&mut line).unwrap() > 0 {
    
        if line.starts_with('>') {

            // Process the previous record if it matches an accession of interest and is not empty
            if is_sequence_of_interest && !current_sequence.is_empty() {
                write_record(&mut writer, &current_accession, &current_sequence);
            }

            // clear accession and sequence contents now that we are onto the next record
            current_sequence.clear();
            current_accession.clear();
    
            // Update the current accession
            let parts: Vec<&str> = line.trim().split_whitespace().collect();
            let potential_accession = &parts[0][1..];
            is_sequence_of_interest = accessions_of_interest.contains(potential_accession);
            if is_sequence_of_interest {
                current_accession = potential_accession.to_string();
            }

        } else if is_sequence_of_interest {

            // Append the sequence line to the current sequence
            current_sequence.push_str(&line.trim());

        }
        
    line.clear();

    }

    // Process the last record in the file if it matches an accession of interest
    if is_sequence_of_interest && !current_sequence.is_empty() {
        write_record(&mut writer, &current_accession, &current_sequence);
    }


    println!("FASTA subsetting completed!");
}

fn load_accessions(filename: &str) -> std::io::Result<HashSet<String>> {
    let file = File::open(filename)?;
    let reader = BufReader::new(file);

    let mut accessions = HashSet::new();
    for line in reader.lines() {
        let line = line?;
        accessions.insert(line);
    }

    Ok(accessions)
}

fn write_record(writer: &mut LineWriter<File>, accession: &str, sequence: &str) {
    writeln!(writer, ">{}", accession).expect("Failed to write record name");

    for chunk in sequence.as_bytes().chunks(80) {
        writeln!(writer, "{}", String::from_utf8_lossy(chunk)).expect("Failed to write sequence");
    }
}
