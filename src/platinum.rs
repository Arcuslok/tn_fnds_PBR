#![allow(unused_assignments)]
use super::constants::{DEFAULT_F0, LOWEST_F0, FLOOR_F0, PI_HALF, TAU};

use ndarray::{Array1, Array2, s};

// Aperiodicity estimation based on PLATINUM
pub fn pt101(x: &Array1<f32>, f0: &mut Array1<f32>, t: &mut Array1<f32>, fs: f32, x_length: usize, t_length: usize) -> (usize, Array2<f32>, Array1<usize>, Array1<usize>) {

    let (frame_period, mut vuv_num): (f32, usize) = (t[1] - t[0], 1);

    f0.iter().zip(f0.iter().skip(1)).for_each(|(f0_p, f0_s)| {
        if *f0_s != 0.0 && *f0_p == 0.0 { vuv_num += 1; } if *f0_s == 0.0 && *f0_p != 0.0 { vuv_num += 1; }
    });
    
    let (mut st_list, mut ed_list): (Array1<usize>, Array1<usize>) = (Array1::<usize>::zeros(vuv_num), Array1::<usize>::zeros(vuv_num));
    let (mut st_count, mut ed_count, mut index): (usize, usize, usize) = (1, 0, 1);
    
    if f0[0] != 0.0 {
        for i in 1 .. t_length {
            if f0[i] == 0.0 && f0[i-1] != 0.0 { 
                (ed_list[0], st_list[1], index) = (i-1, i, i); ed_count += 1; st_count += 1; break; }
        }
    }

    ed_list[vuv_num-1] = t_length-1;

    for i in index .. t_length {
        if f0[i] != 0.0 && f0[i-1] == 0.0 { (ed_list[ed_count], st_list[st_count]) = (i-1, i); ed_count += 1; st_count += 1; }
        if f0[i] == 0.0 && f0[i-1] != 0.0 { (ed_list[ed_count], st_list[st_count]) = (i-1, i); ed_count += 1; st_count += 1; }
    }

    // Wedge List
    let (mut wedge_list, mut temp_wav): (Array1<usize>, Array1<f32>) = (Array1::<usize>::zeros(vuv_num), Array1::<f32>::zeros((fs * 2.0 / LOWEST_F0) as usize));
    let (mut peak_index, mut center, mut t0): (usize, usize, usize);
    let mut peak: f32;
    
    for i in 0 .. vuv_num {

        center = (st_list[i] + ed_list[i] + 1) / 2;
        
        (t0, peak_index) = (((fs / if f0[center] == 0.0 { DEFAULT_F0 } else { f0[center] }) + 0.5) as usize, ((center as f32 * frame_period * fs) + 0.5) as usize);

        temp_wav.slice_mut(s![0 .. t0 * 2 + 1]).assign(&(0 .. t0 * 2 + 1).map(|j: usize| x[0.max((x_length-1).min(peak_index-t0+j-1))]).collect::<Array1<f32>>());

        (peak, peak_index) = (0.0, 0);

        for j in 0 .. t0 * 2 + 1 { if libm::fabsf(temp_wav[j]) > peak { (peak, peak_index) = (temp_wav[j], j); } }

        wedge_list[i] = 0.max((x_length-1).min((0.5 + (center as f32 * frame_period * fs) - t0 as f32 + peak_index as f32 + 1.0) as usize - 1));

    }

    // Linear Interpolation
    let (fixed_f0, mut total_phase, signal_time, mut f0_interpolated_raw, mut pulse_locations): (Array1<f32>, Array1<f32>, Array1<f32>, Array1<f32>, Array1<f32>) = ((0 .. t_length).map(|i| if f0[i] == 0.0 { DEFAULT_F0 } else { f0[i] }).collect::<Array1<f32>>(), Array1::<f32>::zeros(x_length), Array1::<f32>::range(0.0, x_length as f32, 1.0) / fs, Array1::<f32>::zeros(x_length), Array1::<f32>::zeros(x_length));
    let (h, mut k): (Array1<f32>, Array1<usize>) = (t.slice(s![1..]).to_owned() - t.slice(s![..-1]), Array1::<usize>::ones(x_length));
    let (mut j, mut cnt): (usize, usize) = (0, 1);
    
    for i in 0 .. x_length { if signal_time[i] >= t[0] { j = i; break; } }
    
    for mut i in j .. x_length { if signal_time[i] < t[cnt] { k[i] = cnt; } else { k[i] = cnt; cnt += 1; i -= 1; } if cnt == t_length { j = i; break; } }

    cnt -= 1;

    k.slice_mut(s![j + 1 .. x_length]).fill(cnt);

    for i in 0 .. x_length { f0_interpolated_raw[i] = fixed_f0[k[i]-1] + (fixed_f0[k[i]] - fixed_f0[k[i]-1]) * ((signal_time[i] - t[k[i]-1]) / h[k[i]-1]); }

    let ratio: f32 = TAU / fs;

    total_phase[0] = f0_interpolated_raw[0] * ratio;
  
    for i in 1 .. x_length { total_phase[i] = total_phase[i-1] + f0_interpolated_raw[i] * ratio; }

    // Pulse Locations
    let pCount: usize = {
        
        let (mut base_phase, mut tmp_pulse_locations): (Array1<f32>, Array1<f32>) = (Array1::<f32>::zeros(x_length), Array1::<f32>::zeros(x_length));
        let (mut st_index, mut ed_index, mut number_of_location): (usize, usize, usize);
        let mut pulso_temp: usize = 0;
        let mut temp: f32;
        
        for i in 0 .. vuv_num {

            (st_index, ed_index) = (0.max((fs * st_list[i] as f32 * frame_period) as usize), (x_length-1).min((fs * (ed_list[i] + 1) as f32 * frame_period + 0.5) as usize - 1));

            temp = total_phase[wedge_list[i]];

            base_phase.slice_mut(s![st_index .. ed_index]).assign(&(total_phase.slice(s![st_index .. ed_index]).mapv(|x: f32| x - temp + PI_HALF) % TAU));

            number_of_location = 0;

            for (j, k) in (st_index .. ed_index - 1).zip(Array1::<f32>::range(st_index as f32, (ed_index - 1) as f32, 1.0)) {
                if base_phase[j+1] < base_phase[j] { tmp_pulse_locations[number_of_location] = k / fs; number_of_location += 1; }
            }

            pulse_locations.slice_mut(s![0 .. number_of_location]).assign(&tmp_pulse_locations.slice(s![0 .. number_of_location]));

            pulso_temp += number_of_location;

        }

        pulso_temp

    };

    let fftl: usize = 3 * (fs / FLOOR_F0) as usize + 1;

    // One Pulse Residual Signal
    let (mut residual_spectrum, mut residual_spectrum_index, mut residual_spectrum_length): (Array2<f32>, Array1<usize>, Array1<usize>) = (Array2::<f32>::zeros((pCount, fftl)), Array1::<usize>::zeros(t_length), Array1::<usize>::zeros(pCount));
    let (mut index, mut f0si, mut f0ei, mut w_length): (usize, usize, usize, usize);
    let (mut ff0, mut f0fi, mut t0): (f32, f32, f32);

    for j in 0 .. pCount - 1 {
        
        index = 1 + (0.5 + pulse_locations[j] * fs) as usize;
        
        f0fi = pulse_locations[j] / frame_period;
        f0si = f0fi as usize;
        f0ei = f0si + 1;

        ff0 = if f0[f0si] == 0.0 || f0[f0ei] == 0.0 {
            DEFAULT_F0 
        } else if f0[f0si] == 0.0 {
            f0[f0ei]
        } else if f0[f0ei] == 0.0 {
            f0[f0si]
        } else {
            f0[f0si] + (f0[f0ei] - f0[f0si]) * (f0fi - f0si as f32)
        };

        t0 = fs / ff0;
        w_length = (0.5 + t0 * 2.0) as usize;

        if w_length + index - (0.5 + t0) as usize >= x_length { residual_spectrum_length[j] = 0; continue; }

        residual_spectrum_length[j] = (fftl-1).min(w_length);

        residual_spectrum.slice_mut(s![j, 0 .. residual_spectrum_length[j]]).assign(&(0 .. residual_spectrum_length[j]).map(|i: usize| x[(x_length-1).min(if (i + index) < (0.5 + t0) as usize { 0 } else { 0.max((i + index) - (0.5 + t0) as usize) })]).collect::<Array1<f32>>());
    
        residual_spectrum.slice_mut(s![j, residual_spectrum_length[j] .. fftl]).fill(0.0);

    }

    // Frame Residual Index
    let (mut temp, mut tmp_value, mut tmp_index): (f32, f32, usize);

    for (j, k) in (0 .. t_length).zip(Array1::<f32>::range(0.0, t_length as f32, 1.0)) {

        (tmp_value, tmp_index) = (100000.0, pCount - 1);

        for i in 0 .. pCount - 1 {
            temp = libm::fabsf(pulse_locations[i] - k * frame_period);
            if temp < tmp_value { (tmp_value, tmp_index) = (temp, i); }
        }
        
        residual_spectrum_index[j] = tmp_index;

    }

    (pCount, residual_spectrum, residual_spectrum_index, residual_spectrum_length)

}
