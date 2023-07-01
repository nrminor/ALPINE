use std::collections::HashSet;
use std::fs::File;
use std::io::{Read, BufReader, LineWriter, Write, BufRead};
use std::str;

const CHUNK_SIZE: usize = 10_485_760; // 10 MiB

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
    let mut reader = BufReader::with_capacity(CHUNK_SIZE, file);

    // Open the output FASTA file for writing
    let output_file = File::create("filtered-to-geography.fasta").expect("Failed to create output file");
    let mut writer = LineWriter::new(output_file);

    let mut current_accession = String::new();
    let mut current_sequence = String::new();
    let mut is_sequence_of_interest = false;
    let mut buffer = Vec::new();
    let buffer_limit = 5000000;

    let mut buf = vec![0; CHUNK_SIZE];
    let mut leftover_line = String::new();
    while let Ok(bytes_read) = reader.read(&mut buf) {
        if bytes_read == 0 {
            break;
        }
        let mut content = str::from_utf8(&buf[..bytes_read]).unwrap().to_string();
        content.insert_str(0, &leftover_line);
        let mut lines = content.lines().collect::<Vec<_>>();
        if !content.ends_with('\n') {
            leftover_line = lines.pop().unwrap().to_string();
        } else {
            leftover_line.clear();
        }

        for line in lines {
    
            if line.starts_with('>') {

                // Process the previous record if it matches an accession of interest and is not empty
                if is_sequence_of_interest && !current_sequence.is_empty() {
                    write_record(&current_accession, &current_sequence, &mut buffer, buffer_limit, &mut writer);
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
            
        }

    }

    // Process the last record in the file if it matches an accession of interest
    if is_sequence_of_interest && !current_sequence.is_empty() {
        write_record(&current_accession, &current_sequence, &mut buffer, buffer_limit, &mut writer);
    }

    // Flush the buffer at the end
    flush_buffer(&mut writer, &mut buffer);

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

fn write_record(accession: &str, sequence: &str, buffer: &mut Vec<String>, buffer_limit: usize, writer: &mut LineWriter<File>) {
    buffer.push(format!(">{}", accession));

    for chunk in sequence.as_bytes().chunks(80) {
        buffer.push(String::from_utf8_lossy(chunk).into_owned());
    }

    // Check if the buffer size has reached the limit
    if buffer.len() >= buffer_limit {
        flush_buffer(writer, buffer);  // Error occurs here because `writer` is not defined in this scope
    }
}

fn flush_buffer(writer: &mut LineWriter<File>, buffer: &mut Vec<String>) {
    for line in buffer.drain(..) {
        writeln!(writer, "{}", line).expect("Failed to write sequence");
    }
    writer.flush().expect("Failed to flush writer");
}

