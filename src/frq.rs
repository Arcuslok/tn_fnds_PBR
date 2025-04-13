use crate::constants::LOG2_E;
use ndarray::Array1;

pub fn equalizing_pitch(f0: &mut Array1<f32>, t_len: usize, scale_param: &str, modulation: f32, flag_t: f32) {

    // Calculate modulation factor
    let modulation: f32 = modulation / 100.0;

    // Get average F0
    let average_f0: f32 = get_freq_avg(f0, t_len);

    let mut bias = 0;
    let scale_param_chars: Vec<char> = scale_param.chars().collect();
    
    // Determine the bias for sharp notes
    if scale_param_chars.len() > 1 && scale_param_chars[1] == '#' {
        bias = 1;
    }
    // Determine the base scale
    let scale: i32 = match scale_param_chars[0] {
        'C' => -9 + bias,
        'D' => -7 + bias,
        'E' => -5,
        'F' => -4 + bias,
        'G' => -2 + bias,
        'A' => bias,
        'B' => 2,
        _ => 0, // default case (shouldn't happen with valid input)
    };

    // Calculate octave and target frequency
    let octave_char: char = if bias == 1 { scale_param_chars[2] } else { scale_param_chars[1] };
    
    let mut target_f0: f32 = 440.0 * libm::powf(2.0, (octave_char.to_digit(10).unwrap() as i32 - 4) as f32) * libm::powf(2.0, scale as f32 / 12.0);
    target_f0 *= libm::powf(2.0, flag_t / 120.0);

    if average_f0 != 0.0 {

        for i in 0 .. t_len {
            if f0[i] != 0.0 {
                let tmp = (f0[i] - average_f0) * modulation + average_f0;
                f0[i] = tmp * target_f0 / average_f0;
            }
        }

    } else {

        for i in 0 .. t_len {
            if f0[i] != 0.0 {
                f0[i] = target_f0;
            }
        }

    }

}

pub fn get_freq_avg(f0: &Array1<f32>, t_len: usize) -> f32 {

    let mut value: f32;
    let (mut freq_avg, mut base_value): (f32, f32) = (0.0, 0.0);

    for i in 0 .. t_len {

        value = f0[i];

        if value < 1000.0 && value > 55.0 {

            let mut r: f32 = 1.0;

            // Weight calculation for nearby values
            for j in 0..=5 {

                let p_j: f32 = if i > j {
                    let q: f32 = f0[i - j - 1] - value;
                    value / (value + q * q)
                } else {
                    1.0 / (1.0 + value)
                };

                r *= p_j;

            }

            freq_avg += value * r;
            base_value += r;

        }

    }

    if base_value > 0.0 {
        freq_avg / base_value
    } else {
        0.0 // or return freq_avg as is (which would be 0.0)
    }

}

pub fn frq_to_pit(frq: f32) -> f32 {
	libm::logf(frq / 220.0) * LOG2_E * 1200.0 + 5700.0
}

pub fn decpit(s: &str, dst: &mut Array1<f32>) -> usize {

    fn get64(c: char) -> i32 {
    
        match c {
            '0'..='9' => (c as i32) - ('0' as i32) + 52,
            'A'..='Z' => (c as i32) - ('A' as i32),
            'a'..='z' => (c as i32) - ('a' as i32) + 26,
            '+' => 62,
            '/' => 63,
            _ => 0,
        }
    
    }

    let cnt: usize = dst.len();
    let chars: Vec<char> = s.chars().collect();
    let len: usize = chars.len();

    let mut k: usize = 0;
    let mut n: i32 = 0;
    let mut i: usize = 0;
    
    while i < len {
        if chars[i] == '#' {
            i += 1;
            let num_start = i;
            
            // Extract all digits after the #
            while i < len && chars[i].is_ascii_digit() {
                i += 1;
            }
            
            let num_str: String = chars[num_start..i].iter().collect();
            let num: i32 = num_str.parse().unwrap_or(0);
            
            // Repeat the previous value 'num' times
            for _ in 0..num {
                if k >= cnt {
                    break;
                }
                dst[k] = n as f32;
                k += 1;
            }
            
            // Find the next #
            while i < len && chars[i] != '#' {
                i += 1;
            }

        } else {

            // Process normal character pair
            if i + 1 < len {
                let c1 = chars[i];
                let c2 = chars[i + 1];
                n = get64(c1) * 64 + get64(c2);
                
                // Adjust values ​​greater than 2047
                if n > 2047 {
                    n -= 4096;
                }
                
                if k < cnt {
                    dst[k] = n as f32;
                    k += 1;
                }
                i += 1; // We advance an extra character
            }

        }

        i += 1;

    }
    
    len

}
