use std::{fs::File, io::Read};

pub fn get_utau_frq_sub(fname: &str) -> f32 {

    // Try with original filename + ".frq"
    let frq_fname: String = format!("{}.frq", fname);
    let mut fp: File = match File::open(&frq_fname) {
        Ok(file) => file,
        Err(_) => {
            // Try with extension replaced by underscore + ".frq"
            let mut modified_fname: String = fname.to_string();
            if let Some(dot_pos) = fname.rfind('.') {
                modified_fname.replace_range(dot_pos..dot_pos+1, "_");
            }
            let alt_frq_fname: String = format!("{}.frq", modified_fname);
            
            match File::open(&alt_frq_fname) {
                Ok(file) => file,
                Err(_) => {
                    println!("There is no UTAU frequency file.");
                    println!("Please generate .frq files using resampler.exe or fresamp.");
                    return 0.0;
                }
            }
        }
    };

    // Check header
    let mut temp_str: [u8; 4] = [0u8; 4];
    if fp.read_exact(&mut temp_str).is_err() || &temp_str != b"FREQ" {
        println!("Header Error \"FREQ\".");
        return 0.0;
    }

    // Check version
    if fp.read_exact(&mut temp_str).is_err() || &temp_str != b"0003" {
        println!("Header Error Version Problem.");
        return 0.0;
    }

    // Read pitch step (we don't actually use this value)
    let mut pitch_step: [u8; 4] = [0u8; 4];
    if fp.read_exact(&mut pitch_step).is_err() {
        println!("Error reading pitch step.");
        return 0.0;
    }

    // Read average frequency
    let mut avr_fo_bytes: [u8; 8] = [0u8; 8];
    if fp.read_exact(&mut avr_fo_bytes).is_err() {
        println!("Error reading average frequency.");
        return 0.0;
    }

    f64::from_le_bytes(avr_fo_bytes) as f32
    
}

