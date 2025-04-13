#![allow(unused_assignments, unused_mut)]
use tn_fnds_PBR::{constants::{DEFAULT_F0, FLOOR_F0, FRAMEPERIOD, LN_2, TAU}, frq::*, utau::*, platinum::pt101, wham77::wham77};
use std::{io::{BufReader, BufWriter, Read}, time::Instant, path::Path, fs::File, process, env};
use hound::{SampleFormat, WavReader, WavWriter, WavSpec, Error};
use ndarray::{Array1, Array2, s};

fn main() {

    let args: Vec<String> = env::args().collect();

    if args.len() < 3 {
        eprintln!("error: insufficient arguments");
        return;
    }   

    let mut flag_b: usize = 0;
    let mut flag_t: f32 = 0.0;
    let mut flag_g: f32 = 0.0;
    let mut flag_w: f32 = -1.0;
    let mut flag_a: usize = 0;
    let mut g_ratio: f32 = 0.0;

    if args.len() > 5 {

        let arg5: &String = &args[5];

        // Flag b: parameter with range 0-100
        if let Some(pos) = arg5.find('b') {
            if let Ok(num) = arg5[pos+1..].parse::<usize>() {
                flag_b = num.clamp(0, 100);
            }
        }
        
        // Flag t: simple parameter (no range restriction)
        if let Some(pos) = arg5.find('t') {
            if let Ok(num) = arg5[pos+1..].parse::<usize>() {
                flag_t = num as f32;
            }
        }
        
        // Flag g: parameter between -100 and 100, converted to ratio
        if let Some(pos) = arg5.find('g') {
            if let Ok(num) = arg5[pos+1..].parse::<f32>() {
                flag_g = num.clamp(-100.0, 100.0);
            }
        }

        // Using natural exponent (e^x)
        g_ratio = libm::powf(10.0, -flag_g / 200.0);
        
        // Flag W: special parameter with complex rules
        if let Some(pos) = arg5.find('W') {
            if let Ok(num) = arg5[pos+1..].parse::<f32>() {
                flag_w = if num > 1000.0 {
                    1000.0
                } else if num < 50.0 {
                    0.0
                } else {
                    num
                };
            }
        }
        
        // Flag A: parameter with range 0-100
        if let Some(pos) = arg5.find('A') {
            if let Ok(num) = arg5[pos+1..].parse::<usize>() {
                flag_a = num.clamp(0, 100);
            }
        }

    }

    // Parse numeric arguments safely
    let mut offset: f32 = match args[6].parse::<f32>() {
        Ok(num) => num as f32,
        Err(_) => {
            eprintln!("Error: offset must be a valid integer");
            process::exit(1);
        }
    };

    let mut ed_length_msec: f32 = match args[9].parse::<f32>() {
        Ok(num) => num,
        Err(_) => {
            eprintln!("Error: edLengthMsec must be a valid integer");
            process::exit(1);
        }
    };

    let sample: WavReader<BufReader<File>> = WavReader::open(Path::new(&args[1])).unwrap();
    let spec: WavSpec = sample.spec();
    let mut x: Array1<f32> = sample.into_samples::<i16>().into_iter().map(|x: Result<i16, Error>| x.unwrap() as f32 ).collect::<Array1<f32>>();

    let fs: f32 = spec.sample_rate as f32;
    let n_bit: f32 = spec.bits_per_sample as f32;
    
    x = {

        let x_len: i32 = x.len() as i32;
        let (mut cut_offset, mut cut_ed_length_msec) = (offset as i32, ed_length_msec as i32);

        if cut_ed_length_msec < 0 {
            cut_ed_length_msec = (x_len * 1000 / fs as i32) - (cut_offset - cut_ed_length_msec);
        }

        let (st, ed): (i32, i32) = (0.max((x_len - 1).min((cut_offset - 100) * fs as i32 / 1000)), 0.max((x_len - 1).min(x_len - 0.max(cut_ed_length_msec - 100) * fs as i32 / 1000))); 

        offset = (cut_offset - (st * 1000 / fs as i32)) as f32;
        ed_length_msec = ((ed * 1000 / fs as i32) - ((x_len * 1000 / fs as i32) - cut_ed_length_msec)) as f32;

        x.slice(s![st .. ed + 1]).mapv(|x: f32| x / fs)
      
    };

    let x_len: usize = x.len();
    let t_len: usize = (x_len as f32 / fs / (FRAMEPERIOD / 1000.0) ) as usize + 1;

    let mut t: Array1<f32> = Array1::<f32>::range(0.0, t_len as f32, 1.0) * FRAMEPERIOD;
    let mut f0: Array1<f32> = Array1::<f32>::zeros(t_len);
    let mut start_time: Instant = Instant::now();

    println!("File information");
    println!("Sampling : {} Hz {} Bit", fs, n_bit);
    println!("Length {} [sample]", x.len());
    println!("Length {:.6} [sec]", x.len() as f32 / fs);

    if flag_w == -1.0 {

        println!("\nAnalysis");

        //Limits to determine how much trust to be placed on DIO
	    //Advisable for lowerlimit to not be below 0.6, upperlimit not to be above 1.7
	    //I chose slightly tighter limits
        let lower_limit: f32 = 0.7;
        let upper_limit: f32 = 1.6;

        start_time = Instant::now();

        wham77(&x, &mut t, &mut f0, fs, x_len, t_len);

        println!("WHAM: {} [msec]", start_time.elapsed().as_millis());
        
        let utau_avg_f0: f32 = get_utau_frq_sub(&args[1]);
        let mut avg_f0: f32 = get_freq_avg(&f0, t_len);
        
        if utau_avg_f0 != 0.0 {
            
            if 0.95 * utau_avg_f0 > avg_f0 || avg_f0 > 1.05 * utau_avg_f0 {
                avg_f0 = utau_avg_f0;
            }

        }

        f0 = f0.mapv(|f0: f32| if lower_limit * avg_f0 > f0 || f0 > upper_limit * avg_f0 { avg_f0 } else { f0 });

        f0 = f0.mapv(|f0: f32| if f0 == 0.0 { avg_f0 } else { f0 });
        
    } else {
        
        if flag_w == 0.0 {
            flag_w = get_utau_frq_sub(&args[1]);
        }

        f0.fill(flag_w);

    }

    t = Array1::<f32>::range(0.0, t_len as f32, 1.0) * FRAMEPERIOD / 1000.0;

    let fftl: f32 = libm::powf(2.0, 1.0 + libm::roundf(libm::logf(3.0 * fs / FLOOR_F0 + 1.0) / LN_2));

    start_time = Instant::now();

    let (p_count, mut residual_spectrum, mut residual_spectrum_index, mut residual_spectrum_length): (usize, Array2<f32>, Array1<usize>, Array1<usize>) = pt101(&x, &mut f0, &mut t, fs, x_len, t_len);

    println!("PLATINUM: {} [msec]", start_time.elapsed().as_millis());

    if flag_g != 0.0 {

        // g Factor
        for i in 0 .. p_count - 1 {

            if residual_spectrum_length[i] == 0 { continue; }

            let new_length: usize = libm::fmaxf(libm::fminf(libm::roundf(residual_spectrum_length[i] as f32 / g_ratio + 0.5), fftl - 1.0), 0.0) as usize;
 
            if g_ratio > 1.0 {

                for (j, k) in (0 .. new_length).zip(Array1::<f32>::range(0.0, new_length as f32, 1.0)) {

                    let position: f32 = libm::fminf(k * g_ratio, residual_spectrum_length[i] as f32 - 1.0001);
                    let sindex: usize = position as usize;
                    let eindex: usize = sindex + 1;
                    
                    // Linear interpolation
                    residual_spectrum[[i, j]] = residual_spectrum[[i, sindex]] + (residual_spectrum[[i, eindex]] - residual_spectrum[[i, sindex]]) * (position - sindex as f32);

                }

            } else {

                for (j, k) in (0 .. new_length).rev().zip(Array1::<f32>::range(new_length as f32, 0.0, -1.0)) {

                    let position: f32 = libm::fminf(k * g_ratio, residual_spectrum_length[i] as f32 - 1.0001);
                    let sindex: usize = position as usize;
                    let eindex: usize = sindex + 1;
                    
                    // Linear interpolation
                    residual_spectrum[[i, j]] = residual_spectrum[[i, sindex]] + (residual_spectrum[[i, eindex]] - residual_spectrum[[i, sindex]]) * (position - sindex as f32);

                }

            }
            
            residual_spectrum_length[i] = new_length;

        }

    }

    // Pulse Residual Window
    let (mut length, mut coefficient): (f32, f32);

    for i in 0 .. p_count {
        length = residual_spectrum_length[i] as f32 + 1.0;
        coefficient = TAU / length;
        for (j, k) in (0 .. (length - 1.0) as usize).zip(Array1::<f32>::range(0.0, length - 1.0, 1.0)) {
            residual_spectrum[[i, j]] *= 0.5 - 0.5 * libm::cosf((k + 1.0) * coefficient);
        }
    }

    let length_msec: f32 = match args[7].parse::<f32>() {
        Ok(num) => num,
        Err(_) => {
            eprintln!("Invalid number for milliseconds");
            std::process::exit(1);
        }
    };

    let st_length_msec: f32 = match args[8].parse::<f32>() {
        Ok(num) => num,
        Err(_) => {
            eprintln!("Invalid number for milliseconds");
            std::process::exit(1);
        }
    };

    let velocity: f32 = match args[4].parse::<f32>() {
        Ok(num) => num / 100.0,
        Err(_) => {
            eprintln!("Invalid number for milliseconds");
            std::process::exit(1);
        }
    };
    
    let (input_length_msec, x_len2, mut t_len2): (usize, usize, usize) = ((t_len as f32 * FRAMEPERIOD) as usize, (length_msec / 1000.0 * fs) as usize, (length_msec / FRAMEPERIOD) as usize);
    let (mut fixed_f0, mut fixed_volume, mut fixed_residual_specgram_index): (Array1<f32>, Array1<f32>, Array1<usize>) = (Array1::<f32>::zeros(t_len2), Array1::<f32>::ones(t_len2), Array1::<usize>::zeros(t_len2));
    let mut y: Array1<f32> = Array1::<f32>::zeros(x_len2);

    let modulation: f32 = match args[11].parse::<f32>() {
        Ok(num) => num,
        Err(_) => {
            eprintln!("Invalid number for milliseconds");
            std::process::exit(1);
        }
    };
    
    // Equalizing Picth
    equalizing_pitch(&mut f0, t_len, &args[3], modulation, flag_t);

    // Stretch Time
    let v_ratio: f32 = libm::powf(2.0, 1.0 - velocity); 
    let (st, ed, os): (usize, usize, usize) = (((st_length_msec + offset) / FRAMEPERIOD) as usize, (((input_length_msec as f32 - ed_length_msec) / FRAMEPERIOD) as usize).min(t_len - 1), (offset / FRAMEPERIOD) as usize);
    let (st2, ed2): (usize, usize) = (t_len2.min(((st - os) as f32 * v_ratio + 0.5) as usize), t_len2.min((length_msec + 0.5) as usize));
    
    let mut prestretch_position: usize;

    for (i, j) in (0 .. st2).zip(Array1::<f32>::range(0.0, st2 as f32, 1.0)) {
        prestretch_position = 0.max((t_len - 1).min((j / v_ratio) as usize + os));
        (fixed_f0[i], fixed_residual_specgram_index[i]) = (f0[prestretch_position], residual_spectrum_index[prestretch_position]);
    }

    let mut stretching_position: usize = st2;

    while stretching_position < ed2 {

        for forward in st .. ed - 2 {
            if stretching_position > ed2 - 1 { break; }
            (fixed_f0[stretching_position], fixed_residual_specgram_index[stretching_position]) = (f0[forward], residual_spectrum_index[forward]);
            stretching_position += 1;
        }

        for go_back in (st ..= ed - 1).rev() {
            if stretching_position > ed2 - 1 { break; }
            (fixed_f0[stretching_position], fixed_residual_specgram_index[stretching_position]) = (f0[go_back], residual_spectrum_index[go_back]);
            stretching_position += 1;
        }
        
    }


    let mut tempo: f32 = 120.0;
    let mut p_len: usize = t_len2;
    let mut p_step: usize = 256;

    let mut pit: Array1<f32> = Array1::<f32>::zeros(p_len + 1);
    
    if args.len() > 13 {

        let cp: &String = &args[12];
        let tempo_str: &str = &cp[1..]; // Ignore the first character
        tempo = match tempo_str.parse::<f32>() {
            Ok(val) => val,
            Err(e) => {
                eprintln!("Error parsing tempo: {}", e);
                std::process::exit(1);
            }
        };
 
        p_step = libm::roundf(60.0 / 96.0 / tempo * fs) as usize;
        p_len = x_len2 / p_step + 1;

        decpit(&args[13], &mut pit);
            
    }

    let (mut u, mut m): (f32, usize);

    for (i, j) in (0 .. t_len2).zip(Array1::<f32>::range(0.0, t_len2 as f32, 1.0)) {
        u = FRAMEPERIOD * j * 0.001 * fs / p_step as f32;
		m = u as usize; u -= m as f32;
        if m >= p_len { m = p_len - 1; }
        fixed_f0[i] *= libm::powf(2.0, (pit[m] * (1.0 - u) + pit[m + 1] * u) / 1200.0)
    }

    // Auto Volume
    if flag_a != 0 {

        for i in 0..t_len - 1 {

            if fixed_f0[i] == 0.0 {
                fixed_volume[i] = 1.0;
                continue;
            }
    
            if fixed_f0[i + 1] != 0.0 {
                let auto_pow: f32 = (frq_to_pit(fixed_f0[i + 1]) - frq_to_pit(fixed_f0[i])) * (441.0 / (fs * FRAMEPERIOD)) * flag_a as f32;
                fixed_volume[i] = libm::fminf(1.2, libm::powf(2.0, auto_pow));
                continue;
            }
    
            if i > 0 && fixed_f0[i - 1] != 0.0 {
                fixed_volume[i] = fixed_volume[i - 1];
                continue;
            }
    
            fixed_volume[i] = 1.0;

        }
    
        if t_len > 1 && fixed_f0[t_len - 1] != 0.0 && fixed_f0[t_len - 2] != 0.0 {
            fixed_volume[t_len - 1] = fixed_volume[t_len - 2];
        }
        
    }

    // Consonant 2 Amp
    if flag_b != 0 {

        let frame_len: usize = 5; // fixed window size
        let mut add_count: f32 = 0.0;
        let mut add_volume: f32 = 0.0;
        let ratio: f32 = flag_b as f32 / 20.0; // When b=100, ratio=5

        // Initialize the sliding window
        for i in 0 .. t_len.min(frame_len + 1) {
            add_count += 1.0;
            if fixed_f0[i] == 0.0 {
                add_volume += ratio;
            }
        }

        // Process each frame
        for i in 0 .. t_len - 1 {
            fixed_volume[i] *= if add_count != 0.0 {
                add_volume / add_count + 1.0
            } else {
                1.0
            };

            // Update the sliding window
            if i >= frame_len {
                add_count -= 1.0;
                if fixed_f0[i - frame_len] == 0.0 {
                    add_volume -= ratio;
                }
            }
            if i <= t_len - 1 - frame_len - 1 {
                add_count += 1.0;
                if fixed_f0[i + frame_len + 1] == 0.0 {
                    add_volume += ratio;
                }
            }
        }

    }

    let fixed_default_f0: f32 = DEFAULT_F0 * g_ratio;

    println!("Synthesis\n");

    start_time = Instant::now();

    // Synthesis Pt101
    let (mut current_time, mut current_position, mut current_frame, mut final_frame): (f32, usize, usize, usize) = (0.0, 0, 0, residual_spectrum_length[fixed_residual_specgram_index[0]]);
    
    loop {
        
        for j in 0 .. final_frame {
            if j + current_position >= x_len2 { break; }
            y[j + current_position] += residual_spectrum[[fixed_residual_specgram_index[current_frame], j]] * fixed_volume[current_frame];
        }
        
        current_time += 1.0 / if fixed_f0[current_frame] == 0.0 { fixed_default_f0 } else { fixed_f0[current_frame] };
        (current_frame, current_position) = ((current_time / (FRAMEPERIOD / 1000.0)) as usize, (current_time * fs) as usize);
        
        if final_frame + current_position >= x_len2 || current_frame >= t_len2 { break; }

    }

    println!("WORLD: {} [msec]", start_time.elapsed().as_millis());

    // Volume adjustment
    let volume: f32 = args[10].parse::<f32>().unwrap_or(100.0) / 100.0;
    let max_amp: f32 = y.iter().fold(0.0, |max, &val| libm::fmaxf(max, libm::fabsf(val)));
    
    let output: Vec<i16> = y.iter().map(|&val| (32768.0 * (val * 0.5 * volume / max_amp) as f32) as i16).collect::<Vec<i16>>();
    
    // Read header from input file (first 22 bytes)
    let mut header: Vec<u8> = vec![0u8; 44];

    if args.len() > 1 {
        let mut fp = File::open(&args[1]).expect("Failed to open input file");
        fp.read_exact(&mut header[0..22]).expect("Failed to read header");
    }
    
    // Update header values (little-endian)
    // Channels = 1
    header[22..24].copy_from_slice(&1i16.to_le_bytes());
    // Sample rate
    header[24..28].copy_from_slice(&(fs as i32).to_le_bytes());
    // Bytes per second
    header[28..32].copy_from_slice(&((fs * n_bit / 8.0) as i32).to_le_bytes());
    // Block align
    header[32..34].copy_from_slice(&((n_bit / 8.0) as i16).to_le_bytes());
    // Bits per sample
    header[34..36].copy_from_slice(&(n_bit as i16).to_le_bytes());
    // "data" marker
    header[36..40].copy_from_slice(b"data");
    
    // Write WAV file using hound
    let spec: WavSpec = WavSpec {
        channels: 1,
        sample_rate: fs as u32,
        bits_per_sample: n_bit as u16,
        sample_format: SampleFormat::Int,
    };
    
    let mut writer: WavWriter<BufWriter<File>> = WavWriter::create(&args[2], spec).expect("Failed to create WAV writer");

    for sample in output {
        writer.write_sample(sample).expect("Failed to write sample");
    }

    writer.finalize().expect("Failed to finalize WAV file");

    println!("Complete.\n");

}
