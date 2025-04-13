#![allow(unused_assignments)]
use super::constants::{FRAMEPERIOD, CEIL_F0, CHANNELS_IN_OCTAVE, FFT3, FLOOR_F0, FLOOR_F0_W, LN_2, PI, SHRINK, SKIP_RATIO};
use super::matlabfftfunctions::fftfltQ;

use ndarray::{ArrayViewMut1, ArrayView1, Array1, Array2, Axis, s};

pub fn lowpass_f0_FIR_A(f0: &mut ArrayViewMut1<f32>, ratio: f32, fNum: usize) {
    
    let v: f32 = 1.0 / (1.0 + ratio);
    let mut f0s: Array1<f32> = f0.slice(s![0 .. fNum]).to_owned();

    for i in 1 .. fNum {
        if (FLOOR_F0 <= f0[i-1]) && (FLOOR_F0 <= f0[i]) { f0s[i] = (f0[i] + ratio * f0[i-1]) * v; }
    }

    f0.slice_mut(s![0 .. fNum]).assign(&f0s.slice(s![0 .. fNum]));

    for i in (0 ..= fNum-2).rev() {
        if (FLOOR_F0 <= f0[i]) && (FLOOR_F0 <= f0[i+1]) { f0[i] = (f0s[i] + ratio * f0s[i+1]) * v; }
    }

}

pub fn lowpass_f0_IIR_A(f0: &mut ArrayViewMut1<f32>, ratio: f32, fNum: usize) {
    
    let v: f32 = 1.0 / (1.0 + ratio);

    for i in 1 .. fNum {
        if (FLOOR_F0 <= f0[i-1]) && (FLOOR_F0 <= f0[i]) { f0[i] = (f0[i] + ratio * f0[i-1]) * v; }
    }

    for i in (0 ..= fNum-2).rev() {
        if (FLOOR_F0 <= f0[i]) && (FLOOR_F0 <= f0[i+1]) { f0[i] = (f0[i] + ratio * f0[i+1]) * v; }
    }

}

pub fn lowpass_f0_FIR(f0: &mut Array1<f32>, ratio: f32, fNum: usize) {
    
    let v: f32 = 1.0 / (1.0 + ratio);
    let mut f0s: Array1<f32> = f0.slice(s![0 .. fNum]).to_owned();

    for i in 1 .. fNum {
        if (FLOOR_F0 <= f0[i-1]) && (FLOOR_F0 <= f0[i]) { f0s[i] = (f0[i] + ratio * f0[i-1]) * v; }
    }

    f0.slice_mut(s![0 .. fNum]).assign(&f0s.slice(s![0 .. fNum]));

    for i in (0 ..= fNum-2).rev() {
        if (FLOOR_F0 <= f0[i]) && (FLOOR_F0 <= f0[i+1]) { f0[i] = (f0s[i] + ratio * f0s[i+1]) * v; }
    }

}

pub fn lowpass_f0_IIR(f0: &mut Array1<f32>, ratio: f32, fNum: usize) {
    
    let v: f32 = 1.0 / (1.0 + ratio);

    for i in 1 .. fNum {
        if (FLOOR_F0 <= f0[i-1]) && (FLOOR_F0 <= f0[i]) { f0[i] = (f0[i] + ratio * f0[i-1]) * v; }
    }

    for i in (0 ..= fNum-2).rev() {
        if (FLOOR_F0 <= f0[i]) && (FLOOR_F0 <= f0[i+1]) { f0[i] = (f0[i] + ratio * f0[i+1]) * v; }
    }

}


fn decimate_for_estimation(x: &Array1<f32>, y: &mut Array1<f32>, z: &mut Array1<f32>, decimate_num: usize, num_in: usize, num_work: usize) {

    let (index_bias, mut j): (usize, usize) = (SKIP_RATIO * 2 + SKIP_RATIO / 2, SKIP_RATIO * 4 + 1);
    let mut w: Array1<f32> = Array1::<f32>::zeros(FFT3);
    let i: f32 = j as f32;

    // Nuttall window
    if j & 1 == 0 {
        w.slice_mut(s![0 .. j]).assign(&Array1::<f32>::range(0.0, i, 1.0).mapv(|x: f32| 0.355768 - 0.487396 * libm::cosf(2.0 * ((PI * (x + 1.0)) / (i + 1.0))) + 0.144232 * libm::cosf(4.0 * ((PI * (x + 1.0)) / (i + 1.0))) - 0.012604 * libm::cosf(6.0 * ((PI * (x + 1.0)) / (i + 1.0)))));
    } else { 
        w.slice_mut(s![0 .. j]).assign(&Array1::<f32>::range(0.0, i, 1.0).mapv(|x: f32| 0.355768 - 0.487396 * libm::cosf(2.0 * ((PI * (x + 0.5)) / i))  + 0.144232 * libm::cosf(4.0 * ((PI * (x + 0.5)) / i)) - 0.012604 * libm::cosf(6.0 * ((PI * (x + 0.5)) / i))));
    }
    
    fftfltQ(x, &mut w, z, num_in, j, num_work);

    j = index_bias; for i in 0 .. decimate_num { y[i] = z[j]; j += SKIP_RATIO; }

    y.slice_mut(s![decimate_num .. num_work]).fill(0.0);
    
}

fn lowpass_for_estimation(y: &mut Array1<f32>, z: &mut Array1<f32>, fc: f32, fs: f32, num_in: usize, num_out: usize) {

    let half_average_length: usize = ((fs / fc) * 0.5 + 0.5) as usize;
    let (index_bias, j): (usize, usize) = (half_average_length * 2, half_average_length * 4);
    let mut w: Array1<f32> = Array1::<f32>::zeros(FFT3);
    let i: f32 = j as f32;
    
    // Nuttall window
    if j & 1 == 0 {
        w.slice_mut(s![0 .. j]).assign(&Array1::<f32>::range(0.0, i, 1.0).mapv(|x: f32| 0.355768 - 0.487396 * libm::cosf(2.0 * ((PI * (x + 1.0)) / (i + 1.0))) + 0.144232 * libm::cosf(4.0 * ((PI * (x + 1.0)) / (i + 1.0))) - 0.012604 * libm::cosf(6.0 * ((PI * (x + 1.0)) / (i + 1.0)))));
    } else { 
        w.slice_mut(s![0 .. j]).assign(&Array1::<f32>::range(0.0, i, 1.0).mapv(|x: f32| 0.355768 - 0.487396 * libm::cosf(2.0 * ((PI * (x + 0.5)) / i))  + 0.144232 * libm::cosf(4.0 * ((PI * (x + 0.5)) / i)) - 0.012604 * libm::cosf(6.0 * ((PI * (x + 0.5)) / i))));
    }
    
    fftfltQ(y, &mut w, z, num_in, j, num_out);

    for (i, k) in (0 .. num_in).zip(index_bias .. ) { z[i] = z[k]; }

    z.slice_mut(s![num_in .. num_out]).fill(0.0);
    
}

fn zero_crossing_engine(y: &mut Array1<f32>, fine_edges: &mut Array1<f32>, i_locations: &mut ArrayViewMut1<f32>, intervals: &mut ArrayViewMut1<f32>, i_num: &mut usize, fs: f32, s_num: usize) {
    
    let mut cnt: usize = 0;
    
    for (i, j) in (0 .. s_num - 1).zip(Array1::<f32>::range(0.0, (s_num - 1) as f32, 1.0).iter()) {
        if (0.0 < y[i]) && (y[i+1] <= 0.0) { fine_edges[cnt] = j + y[i] / (y[i] - y[i+1]); cnt += 1; }
    }

    if 3 <= cnt {
        for i in 0 .. cnt - 1 { (i_locations[i], intervals[i]) = (1000.0 * ((fine_edges[i] + fine_edges[i+1]) * 0.5) / fs, fs / (fine_edges[i+1] - fine_edges[i])); }  
        *i_num = cnt - 1;
    }

    else { *i_num = 0; }

}

fn interp1trim_clip(x: ArrayView1<f32>, y: ArrayView1<f32>, x0: &mut ArrayViewMut1<f32>, y0: &mut ArrayViewMut1<f32>, n: usize, m: usize) {

    let (mut a, mut b, mut c, mut h): (Array1<f32>, Array1<f32>, Array1<f32>, Array1<f32>) = (Array1::<f32>::zeros(n), Array1::<f32>::zeros(n), Array1::<f32>::zeros(n), Array1::<f32>::zeros(n));
    let (mut index, mut pnt0): (usize, usize) = (0, 0);
    let mut u: f32;

    if n < 4 {

        let (mut dy, mut vy): (f32, f32);

        for i in 0 .. n - 1 {
            dy = y[i+1] - y[i]; vy = 0.5 * dy;
            (a[i], b[i], c[i]) = (vy + vy - 2.0 * dy, 3.0 * dy - 2.0 * vy - vy, vy);
        }

        for i in 0 .. n - 1 { h[i] = x[i+1] - x[i]; }

        while pnt0 < m { if x[index] < x0[pnt0] { break; } y0[pnt0] = y[index]; pnt0 += 1; }

        'internal: while pnt0 < m {

            while x[index] < x0[pnt0] { index += 1; if n-1 < index { break 'internal; } }
    
            u = (x0[pnt0] - x[index-1]) / h[index-1];
            y0[pnt0] = y[index-1] + u * (c[index-1] + u * (b[index-1] + u * a[index-1]));
            pnt0 += 1;

        }

        index = n - 1;
        
        while pnt0 < m { y0[pnt0] = y[index]; pnt0 += 1; }

        return;

    }

    let (mut dy, mut vy, mut vy_a, mut vy_m): (Array1<f32>, Array1<f32>, Array1<f32>, Array1<f32>) = (Array1::<f32>::zeros(n), Array1::<f32>::zeros(n), Array1::<f32>::zeros(n), Array1::<f32>::zeros(n));
    let (mut vy1, mut vy2): (f32, f32);

    for i in 0 .. n - 1 { (h[i], dy[i]) = (x[i+1] - x[i], y[i+1] - y[i]);  vy[i] = dy[i] / h[i]; }

    (h[n-1], dy[n-1], vy[n-1]) = (h[n-2], dy[n-2] / 7.0, dy[n-1] / h[n-1]);
    (vy_a[0], vy_m[0]) = (vy[0] + (vy[0] / 7.0), vy[0] * (vy[0] / 7.0));

    for i in 1 .. n { (vy_a[i], vy_m[i]) = (vy[i] + vy[i-1], vy[i] * vy[i-1]); }

    if 0.0 == vy_m[0] { vy1 = 0.0; }
    
    else if 0.0 < vy_m[0] { vy1 = 4.0 * vy_m[0] / vy_a[0]; }

    else { vy1 = 0.0; }


    if 0.0 == vy_m[1] { vy2 = 0.0; }
    
    else if 0.0 < vy_m[1] { vy2 = 4.0 * vy_m[1] / vy_a[1]; }
    
    else if 0.0 < vy[1] * vy_a[1] { vy2 = -2.0 * vy_a[1] * vy_a[1] * vy[0] / (vy[1] * vy[1]); }
    
    else { vy2 = -2.0 * vy_a[1] * vy_a[1] * vy[1] / (vy[0] * vy[0]); }
    
    (a[0], b[0], c[0]) = (vy1 * h[0] + vy2 * h[0] - 2.0 * dy[0], 2.0 * dy[0] - 2.0 * vy1 * h[0] - vy2 * h[0], vy1 * h[0]);

    for i in 1 .. n - 1 {
        
        if 0.0 == vy_m[i] { vy1 = 0.0; }
        
        else if 0.0 < vy_m[i] { vy1 = 4.0 * vy_m[i] / vy_a[i]; }
        
        else if 0.0 < vy[i] * vy_a[i] { vy1 = -2.0 * vy_a[i] * vy_a[i] * vy[i-1] / (vy[i] * vy[i]); }
        
        else { vy1 = -2.0 * vy_a[i] * vy_a[i] * vy[i] / (vy[i-1] * vy[i-1]); }
        

        if 0.0 == vy_m[i+1] { vy2 = 0.0; }
        
        else if 0.0 < vy_m[i+1] { vy2 = 4.0 * vy_m[i+1] / vy_a[i+1]; }
        
        else if 0.0 < vy[i+1] * vy_a[i+1] { vy2 = -2.0 * vy_a[i+1] * vy_a[i+1] * vy[i] / (vy[i+1] * vy[i+1]); }
        
        else { vy2 = -2.0 * vy_a[i+1] * vy_a[i+1] * vy[i+1] / (vy[i] * vy[i]); }

        (a[i], b[i], c[i]) = (vy1 * h[i] + vy2 * h[i] - 2.0 * dy[i], 2.0 * dy[i] - 2.0 * vy1 * h[i] - vy2 * h[i], vy1 * h[i]);

    }

    while pnt0 < m { if x[index] < x0[pnt0] { break; } y0[pnt0] = y[index]; pnt0 += 1; }

    'internal: while pnt0 < m {

        while x[index] < x0[pnt0] { index += 1; if n-1 < index { break 'internal; } }

        u = (x0[pnt0] - x[index-1]) / h[index-1];
        y0[pnt0] = y[index-1] + u * (c[index-1] + u * (b[index-1] + u * a[index-1]));
        pnt0 += 1;

    }

    index = n - 1;
    
    while pnt0 < m { y0[pnt0] = y[index]; pnt0 += 1; }

}


fn make_interpolated_f0(boundary_f0: f32, t_in: &mut ArrayViewMut1<f32>, t_out: &mut Array1<f32>, f0_in: &mut ArrayViewMut1<f32>, f0_out: &mut ArrayViewMut1<f32>, num_in: usize, num_out: usize) {
    
    for i in 0 .. num_in {
        if f0_in[i] < boundary_f0 * 0.5 || boundary_f0 < f0_in[i] || f0_in[i] < FLOOR_F0 || CEIL_F0 < f0_in[i] {
            f0_in[i] = 0.0;
        }
    } 
    
    lowpass_f0_FIR_A(f0_in, 0.2, num_in);
    lowpass_f0_IIR_A(f0_in, 0.1, num_in);
    lowpass_f0_FIR_A(f0_in, 0.3, num_in);

    let (mut v_index, mut u_index): (Array1<usize>, Array1<usize>) = (Array1::<usize>::zeros(num_in), Array1::<usize>::zeros(num_in));
    let (mut v_cnt, mut u_cnt): (usize, usize) = (0, 0);

    let (mut st_index, mut ed_index): (usize, usize);
    let (mut left, mut mid, mut right): (usize, usize, usize);

    if FLOOR_F0 <= f0_in[0] { v_index[v_cnt] = 0; v_cnt += 1; }

    for i in 1 .. num_in {

        if (f0_in[i-1] < FLOOR_F0) && (FLOOR_F0 <= f0_in[i]) { v_index[v_cnt] = i; v_cnt += 1; }

        else if (FLOOR_F0 <= f0_in[i-1]) && (f0_in[i] < FLOOR_F0) { u_index[u_cnt] = i-1; u_cnt += 1;  }

    }

    u_index[u_cnt] = num_in-1;
    u_cnt += 1; 

    f0_out.slice_mut(s![0 .. num_out]).fill(0.0);

    (st_index, ed_index) = (0, 0);

    for i in 0 .. v_cnt {
        
        if u_index[i] - v_index[i] < 3 { continue; }

        t_in[v_index[i]] -= 500.0 / f0_in[v_index[i]];
        t_in[u_index[i]] += 500.0 / f0_in[u_index[i]];

        (left, right) = (st_index, num_out-1);

        while left < right {

            mid = (left + right) / 2;

            if t_out[mid] < t_in[v_index[i]] { left = mid + 1; }
            
            else { right = mid; }

        }

        st_index = left;

        if 0 < st_index { st_index -= 1; }
        if num_out-2 <= st_index { continue; }

        (left, right) = (st_index, num_out-1);

        while left < right {

            mid = (left + right) / 2;

            if t_out[mid] < t_in[u_index[i]] { left = mid + 1; }
            
            else { right = mid; }

        }

        ed_index = left;

        interp1trim_clip(t_in.slice(s![v_index[i] .. ]), f0_in.slice(s![v_index[i] .. ]), &mut t_out.slice_mut(s![st_index .. ]), &mut f0_out.slice_mut(s![st_index .. ]), u_index[i] - v_index[i] + 1, ed_index - st_index + 1);

    }

    lowpass_f0_IIR_A(f0_out, 0.1, num_out);
    lowpass_f0_IIR_A(f0_out, 0.1, num_out);
    
}

fn make_f0_map(x: &Array1<f32>, t: &mut Array1<f32>, f0_deviations: &mut Array1<f32>, interpolated_f0: &mut Array1<f32>, boundary_f0: f32, fs: f32, f_num: usize, s_num: usize, offset_time: f32) {

    let (mut m_intervals_locations, mut m_intervals, mut tmp_engine): (Array2<f32>, Array2<f32>, Array1<f32>) = (Array2::<f32>::zeros((4, s_num)), Array2::<f32>::zeros((4, s_num)), Array1::<f32>::zeros(s_num));
    let (mut tmp_wave, mut m_num): (Array1<f32>, [usize; 4]) = (x.slice(s![0 .. s_num]).to_owned(), [0, 0, 0, 0]);
    
    zero_crossing_engine(&mut tmp_wave, &mut tmp_engine, &mut m_intervals_locations.index_axis_mut(Axis(0), 0), &mut m_intervals.index_axis_mut(Axis(0), 0), &mut m_num[0], fs, s_num);
    
    tmp_wave = -x.slice(s![0 .. s_num]).to_owned();

    zero_crossing_engine(&mut tmp_wave, &mut tmp_engine, &mut m_intervals_locations.index_axis_mut(Axis(0), 1), &mut m_intervals.index_axis_mut(Axis(0), 1), &mut m_num[1], fs, s_num);
    
    tmp_wave[0] = (x[1] - x[0]) * 0.5;

    for i in 1 .. s_num - 1 { tmp_wave[i] = (x[i+1] - x[i-1]) * 0.5; }
    
    tmp_wave[s_num-1] = (x[s_num-1] - x[s_num-2]) * 0.5;

    zero_crossing_engine(&mut tmp_wave, &mut tmp_engine, &mut m_intervals_locations.index_axis_mut(Axis(0), 2), &mut m_intervals.index_axis_mut(Axis(0), 2), &mut m_num[2], fs, s_num);

    tmp_wave = -tmp_wave;

    zero_crossing_engine(&mut tmp_wave, &mut tmp_engine, &mut m_intervals_locations.index_axis_mut(Axis(0), 3), &mut m_intervals.index_axis_mut(Axis(0), 3), &mut m_num[3], fs, s_num);
    
    interpolated_f0.slice_mut(s![0 .. f_num]).fill(0.0);
    f0_deviations.slice_mut(s![0 .. f_num]).fill(1000000.0);

    // Usable Channel
    if m_num.iter().product::<usize>() == 0 { return; }

    let mut itpled_f0_set: Array2<f32> = Array2::<f32>::zeros((4, f_num));
    
    for i in 0 .. 4 {
        for j in 0 .. m_num[i] { m_intervals_locations[[i, j]] += offset_time; }
        make_interpolated_f0(boundary_f0, &mut m_intervals_locations.index_axis_mut(Axis(0), i),t, &mut m_intervals.index_axis_mut(Axis(0), i), &mut itpled_f0_set.index_axis_mut(Axis(0), i), m_num[i], f_num);
    }
    
    let mut itpled_f0: f32;

    for i in 0 .. f_num {

        if itpled_f0_set[[0, i]] < boundary_f0 * 0.5 || boundary_f0 < itpled_f0_set[[0, i]] || itpled_f0_set[[0, i]] < FLOOR_F0 || CEIL_F0 < itpled_f0_set[[0, i]] { continue; }

        if itpled_f0_set[[1, i]] < boundary_f0 * 0.5 || boundary_f0 < itpled_f0_set[[1, i]] || itpled_f0_set[[1, i]] < FLOOR_F0 || CEIL_F0 < itpled_f0_set[[1, i]] { continue; }

        if itpled_f0_set[[2, i]] < boundary_f0 * 0.5 || boundary_f0 < itpled_f0_set[[2, i]] || itpled_f0_set[[2, i]] < FLOOR_F0 || CEIL_F0 < itpled_f0_set[[2, i]] { continue; }

        if itpled_f0_set[[3, i]] < boundary_f0 * 0.5 || boundary_f0 < itpled_f0_set[[3, i]] || itpled_f0_set[[3, i]] < FLOOR_F0 || CEIL_F0 < itpled_f0_set[[3, i]] { continue; }

        itpled_f0 = (itpled_f0_set[[0, i]] + itpled_f0_set[[1, i]] + itpled_f0_set[[2, i]] + itpled_f0_set[[3, i]]) * 0.25;
        interpolated_f0[i] = itpled_f0;
        f0_deviations[i] = libm::sqrtf(((itpled_f0_set[[0, i]] - itpled_f0) * (itpled_f0_set[[0, i]] - itpled_f0) + (itpled_f0_set[[1, i]] - itpled_f0) * (itpled_f0_set[[1, i]] - itpled_f0) + (itpled_f0_set[[2, i]] - itpled_f0) * (itpled_f0_set[[2, i]] - itpled_f0) + (itpled_f0_set[[3, i]] - itpled_f0) * (itpled_f0_set[[3, i]] - itpled_f0)) / 3.0);

    }
    
    lowpass_f0_IIR(interpolated_f0, 0.1, f_num);
    lowpass_f0_IIR(interpolated_f0, 0.1, f_num);
    
}

fn combsort(x: &mut Array1<f32>, n: usize) {

    let mut gap: usize = n;
    let mut swapped: bool;

    while gap > 1 {

        gap = (gap as f32 / SHRINK) as usize;
        if gap < 1 { gap = 1; }
        swapped = false;

        for i in 0 .. (n - gap) {     
            if x[i] > x[i + gap] {
                x.swap(i, i + gap);
                swapped = true;
            }
        }

        if !swapped && gap <= 1 { break; }

    }

}


// F0 estimation by WHAM (Wave History Adaption Method)
// x            : Input signal (double *)
// t            : Time axis. f0[0] is the F0 estimated at timeAxis[0]. [msec]
// f0           : Estimated result [Hz].
// fs           : Sampling frequency [Hz]
// sNum         : number of signal [sample]
// tNum         : number of time axis [sample]
pub fn wham77(x: &Array1<f32>, t: &mut Array1<f32>, f0: &mut Array1<f32>, fs: f32, x_length: usize, t_length: usize) {

    let (downed_fs, n_bands): (f32, usize) = (fs / SKIP_RATIO as f32, (libm::logf(1.4 * CEIL_F0 / FLOOR_F0_W) / LN_2 * CHANNELS_IN_OCTAVE) as usize);
    let candidates: usize = n_bands * SKIP_RATIO * 2;

    let boundary_f0_list: Array1<f32> = Array1::<f32>::range(0.0, n_bands as f32, 1.0).mapv(|x: f32| FLOOR_F0_W * libm::powf(2.0, (x + 1.0) / CHANNELS_IN_OCTAVE));
    
    let (yn, zn): (usize, usize) = (x_length / SKIP_RATIO + 1, x_length + 8 * 4 * ((fs / boundary_f0_list[0]) * 0.5 + 0.5) as usize);
    
    let (mut f0_map, mut stability_map): (Array2<f32>, Array2<f32>) = (Array2::<f32>::zeros((candidates, t_length)), Array2::<f32>::zeros((candidates, t_length)));
    let (mut f0_deviations, mut interpolated_f0): (Array1<f32>, Array1<f32>) = (Array1::<f32>::zeros(t_length), Array1::<f32>::zeros(t_length));
    let (mut y, mut z): (Array1<f32>, Array1<f32>) = (Array1::<f32>::zeros(zn), Array1::<f32>::zeros(zn));
    let (mut shifted_x, mut tmp_shifted_x, yn_f): (Array1<f32>, Array1<f32>, f32) = (x.to_owned(), x.to_owned(), yn as f32);
    let mut mean_y: f32;
    
    for (cambio_num, cambio_num_f) in (0 .. SKIP_RATIO).zip(Array1::<f32>::range(0.0, SKIP_RATIO as f32, 1.0)) {
        
        if 1 < SKIP_RATIO { decimate_for_estimation(&shifted_x, &mut y, &mut z, yn, x_length, zn); }

        else { y = shifted_x.to_owned(); }

        mean_y = y.slice(s![0 .. yn]).sum() / yn_f;

        y.slice_mut(s![0 .. yn]).mapv_inplace(|x: f32| x - mean_y);

        for i in 0 .. n_bands {

            lowpass_for_estimation(&mut y, &mut z, boundary_f0_list[i], downed_fs, yn, zn);

            make_f0_map(&z, t, &mut f0_deviations, &mut interpolated_f0, boundary_f0_list[i], downed_fs, t_length, yn, 1000.0 * cambio_num_f / fs);

            f0_map.slice_mut(s![cambio_num * n_bands + i, 0 .. t_length]).assign(&interpolated_f0.slice(s![0 .. t_length]).to_owned());
            
            stability_map.slice_mut(s![cambio_num * n_bands + i, 0 .. t_length]).assign(&(f0_deviations.slice(s![0 .. t_length]).to_owned() / interpolated_f0. slice(s![0 .. t_length])));


        }

        tmp_shifted_x = shifted_x.to_owned();
        shifted_x.slice_mut(s![0 .. x_length - 1]).assign(&tmp_shifted_x.slice(s![1 .. x_length]));
    
    }

    shifted_x.slice_mut(s![0 .. x_length]).assign(&x.slice(s![..;-1]));

    for (cambio_num, cambio_num_f) in (0 .. SKIP_RATIO).zip(Array1::<f32>::range(0.0, SKIP_RATIO as f32, 1.0)) {

        if 1 < SKIP_RATIO { decimate_for_estimation(&shifted_x, &mut y, &mut z, yn, x_length, zn); }

        else { y = shifted_x.to_owned(); }

        mean_y = y.slice(s![0 .. yn]).sum() / yn_f;
        
        y.slice_mut(s![0 .. yn]).mapv_inplace(|x: f32| x - mean_y);

        for i in 0 .. n_bands {

            lowpass_for_estimation(&mut y, &mut z, boundary_f0_list[i], downed_fs, yn, zn);
            
            make_f0_map(&z, t, &mut f0_deviations, &mut interpolated_f0, boundary_f0_list[i], downed_fs, t_length, yn, 1000.0 * cambio_num_f / fs);

            f0_map.slice_mut(s![(SKIP_RATIO * n_bands) + cambio_num * n_bands + i, 0 .. t_length]).assign(&interpolated_f0.slice(s![0 .. t_length]).slice(s![..;-1]));

            stability_map.slice_mut(s![(SKIP_RATIO * n_bands) + cambio_num * n_bands + i, 0 .. t_length]).assign(
                &(f0_deviations.slice(s![0 .. t_length]).slice(s![..;-1]).to_owned() / interpolated_f0.slice(s![0 .. t_length]).slice(s![..;-1]))
            );

        }

        tmp_shifted_x = shifted_x.to_owned();
        shifted_x.slice_mut(s![0 .. x_length - 1]).assign(&tmp_shifted_x.slice(s![1 .. x_length]));

    }

    // Sift f0 map
    
    f0_map.slice_mut(s![0 .. candidates, 0]).fill(0.0);
    f0_map.slice_mut(s![0 .. candidates, t_length-1]).fill(0.0);
    stability_map.slice_mut(s![0 .. candidates, 0]).fill(1000000.0);
    stability_map.slice_mut(s![0 .. candidates, t_length-1]).fill(1000000.0);

    let (mut f0_step0, mut f0_step1, mut f0_step2, mut f0_step3, mut f0_step4, mut f0_step5): (Array1<f32>, Array1<f32>, Array1<f32>, Array1<f32>, Array1<f32>, Array1<f32>) = (Array1::<f32>::zeros(t_length), Array1::<f32>::zeros(t_length), Array1::<f32>::zeros(t_length), Array1::<f32>::zeros(t_length), Array1::<f32>::zeros(t_length), Array1::<f32>::zeros(t_length));
   
    // Sort Map
    let (mut f0s, mut gap, mut cnt): (Array1<f32>, usize, usize);
    let (mut avr_f0, mut l_lim, mut h_lim, mut tmp): (f32, f32, f32, f32);
    let mut swapped: bool;
    
    for i in 0 .. t_length {

        (gap, swapped) = (candidates, false);

        loop {

            gap = (gap as f32 / SHRINK) as usize;
            if gap < 1 { gap = 1; }
            swapped = false;
            
            for j in 0 .. candidates - gap {

                if stability_map[[gap+j, i]] < stability_map[[j, i]] {
                    
                    tmp = stability_map[[j, i]];

                    (stability_map[[j, i]], stability_map[[gap+j, i]]) = (stability_map[[gap+j, i]], tmp);

                    tmp = f0_map[[j, i]];

                    (f0_map[[j, i]], f0_map[[gap+j, i]]) = (f0_map[[gap+j, i]], tmp);

                    swapped = true;

                }
                
            }

            if !swapped || 1 >= gap { break; }

        }

    }

    // Median map
    (f0s, cnt) = (Array1::<f32>::zeros(5 * t_length), 0);
    for i in 0 .. 5 {
        for j in 0 .. t_length {
            if FLOOR_F0 <= f0_map[[i, j]] {
                f0s[cnt] = f0_map[[i, j]];
                cnt += 1;
            } 
        }
    }
    
    avr_f0 = 0.0;

    if cnt > 0 {
        
        combsort(&mut f0s, cnt);
        avr_f0 = if cnt & 1 == 0 {
            cnt /= 2;
            (f0s[cnt] + f0s[cnt-1]) * 0.5
        }

        else{
            cnt /= 2;
            f0s[cnt]
        };

    }

    // Filter
    (l_lim, h_lim) = (0.56123102415, 1.78179743628);

    for i in 0 .. t_length {
        for j in 0 .. candidates {
            if (f0_map[[j, i]] < l_lim * avr_f0) || (h_lim * avr_f0 < f0_map[[j, i]]) {
                (f0_map[[j, i]], stability_map[[j, i]]) = (0.0, 1000000.0);
            }
        }
    }

    // Sort Map
    for i in 0 .. t_length {

        (gap, swapped) = (candidates, false);

        loop {

            gap = (gap as f32 / SHRINK) as usize;
            if gap < 1 { gap = 1; }
            swapped = false;

            for j in 0 .. candidates - gap {

                if stability_map[[gap+j, i]] < stability_map[[j, i]] {
                    
                    tmp = stability_map[[j, i]];

                    (stability_map[[j, i]], stability_map[[gap+j, i]]) = (stability_map[[gap+j, i]], tmp);

                    tmp = f0_map[[j, i]];

                    (f0_map[[j, i]], f0_map[[gap+j, i]]) = (f0_map[[gap+j, i]], tmp);

                    swapped = true;

                }
                
            }

            if !swapped || 1 >= gap { break; }

        }

    }

    // Median map
    for i in 0 .. t_length {

        f0_step1[i] = {

            (f0s, cnt) = (Array1::<f32>::zeros(5), 0);

            for j in 0 .. 5 {
                if FLOOR_F0 <= f0_map[[j, i]] {
                    f0s[cnt] = f0_map[[j, i]];
                    cnt += 1;
                }
            }

            avr_f0 = 0.0;

            if cnt > 0 {
                
                combsort(&mut f0s, cnt);

                avr_f0 = if cnt & 1 == 0 {
                    cnt /= 2;
                    (f0s[cnt] + f0s[cnt-1]) * 0.5
                }

                else{
                    cnt /= 2;
                    f0s[cnt]
                };

            }

            avr_f0

        };

    }

    // Filter
    lowpass_f0_IIR(&mut f0_step1, 0.7, t_length);
    lowpass_f0_IIR(&mut f0_step1, 0.7, t_length);
    lowpass_f0_IIR(&mut f0_step1, 0.7, t_length);

    (l_lim, h_lim) = (0.79370052598, 1.25992104989);

    for i in 0 .. t_length {
        for j in 0 .. candidates {
            if (f0_map[[j, i]] < l_lim * f0_step1[i]) || (h_lim * f0_step1[i] < f0_map[[j, i]]) {
                (f0_map[[j, i]], stability_map[[j, i]]) = (0.0, 1000000.0);
            }
        }
    }

    // Sort map
    for i in 0 .. t_length {

        (gap, swapped) = (candidates, false);

        loop {

            gap = (gap as f32 / SHRINK) as usize;
            if gap < 1 { gap = 1; }
            swapped = false;

            for j in 0 .. candidates - gap {
                
                if stability_map[[gap+j, i]] < stability_map[[j, i]] {
                    
                    tmp = stability_map[[j, i]];

                    (stability_map[[j, i]], stability_map[[gap+j, i]]) = (stability_map[[gap+j, i]], tmp);

                    tmp = f0_map[[j, i]];

                    (f0_map[[j, i]], f0_map[[gap+j, i]]) = (f0_map[[gap+j, i]], tmp);

                    swapped = true;

                }
                
            }

            if !swapped || 1 >= gap { break; }

        }

    }

    // Median map
    for i in 0 .. t_length {

        f0_step2[i] = {

            (f0s, cnt) = (Array1::<f32>::zeros(5), 0);

            for j in 0 .. 5 {
                if FLOOR_F0 <= f0_map[[j, i]] {
                    f0s[cnt] = f0_map[[j, i]];
                    cnt += 1;
                }
            }

            avr_f0 = 0.0;

            if cnt > 0 {
                
                combsort(&mut f0s, cnt);

                avr_f0 = if cnt & 1 == 0 {
                    cnt /= 2;
                    (f0s[cnt] + f0s[cnt-1]) * 0.5
                }

                else{
                    cnt /= 2;
                    f0s[cnt]
                } 

            }

            avr_f0

        };

    }

    // Filter
    lowpass_f0_IIR(&mut f0_step1, 0.3, t_length);
    lowpass_f0_IIR(&mut f0_step1, 0.3, t_length);
    lowpass_f0_IIR(&mut f0_step1, 0.3, t_length);

    (l_lim, h_lim) = (0.89089871814, 1.12246204831);

    for i in 0 .. t_length {
        for j in 0 .. candidates {
            if (f0_map[[j, i]] < l_lim * f0_step2[i]) || (h_lim * f0_step2[i] < f0_map[[j, i]]) {
                (f0_map[[j, i]], stability_map[[j, i]]) = (0.0, 1000000.0);
            }
        }
    }

    // Sort map
    for i in 0 .. t_length {

        (gap, swapped) = (candidates, false);

        loop {

            gap = (gap as f32 / SHRINK) as usize;
            if gap < 1 { gap = 1; }
            swapped = false;

            for j in 0 .. candidates - gap {

                if stability_map[[gap+j, i]] < stability_map[[j, i]] {
                    
                    tmp = stability_map[[j, i]];

                    (stability_map[[j, i]], stability_map[[gap+j, i]]) = (stability_map[[gap+j, i]], tmp);

                    tmp = f0_map[[j, i]];

                    (f0_map[[j, i]], f0_map[[gap+j, i]]) = (f0_map[[gap+j, i]], tmp);

                    swapped = true;

                }
                
            }

            if !swapped || 1 >= gap { break; }

        }

    }

    // Median map
    for i in 0 .. t_length {

        f0[i] = {

            (f0s, cnt) = (Array1::<f32>::zeros(4), 0);

            for j in 0 .. 4 {
                if FLOOR_F0 <= f0_map[[j, i]] {
                    f0s[cnt] = f0_map[[j, i]];
                    cnt += 1;
                }
            }

            avr_f0 = 0.0;

            if cnt > 0 {
                
                combsort(&mut f0s, cnt);

                avr_f0 = if cnt & 1 == 0 {
                    cnt /= 2;
                    (f0s[cnt] + f0s[cnt-1]) * 0.5
                }

                else{
                    cnt /= 2;
                    f0s[cnt]
                };

            }

            avr_f0

        };

    }

    // Filter
    lowpass_f0_IIR(f0, 0.1, t_length);
    lowpass_f0_IIR(f0, 0.1, t_length);

    // Prune best f0

    // 1st step
    f0_step0.slice_mut(s![1 .. t_length]).assign(&f0.slice(s![1 .. t_length]));

    // Filter
    let (l_lim, h_lim): (f32, f32) = (libm::powf(2.0, (-200.0 * FRAMEPERIOD * 0.2) / 1200.0), libm::powf(2.0,  (200.0 * FRAMEPERIOD * 0.2) / 1200.0));

    f0_step1[0] = f0_step0[0];

    for i in 1 .. t_length {
        if f0_step1[i-1] < FLOOR_F0 { f0_step1[i] = f0_step0[i]; }
        else if (l_lim * f0_step1[i-1] <= f0_step0[i]) && (f0_step0[i] <= h_lim * f0_step1[i-1]) { f0_step1[i] = f0_step0[i]; }
    }

    // 2nd step
    let (mut v_index, mut u_index): (Array1<usize>, Array1<usize>) = (Array1::<usize>::zeros(t_length), Array1::<usize>::zeros(t_length));
    let (mut v_cnt, mut u_cnt): (usize, usize) = (0, 0);
    
    if FLOOR_F0 <= f0_step1[0] { v_index[v_cnt] = 0; v_cnt += 1; }

    for i in 1 .. t_length {
        if (f0_step1[i-1] < FLOOR_F0) && (FLOOR_F0 <= f0_step1[i]) { v_index[v_cnt] = i; v_cnt += 1; }
        else if (FLOOR_F0 <= f0_step1[i-1]) && (f0_step1[i] < FLOOR_F0) { u_index[u_cnt] = i-1; u_cnt += 1; }
    }

    u_index[u_cnt] = t_length-1;
    u_cnt += 1;

    f0_step2.slice_mut(s![0 .. t_length]).assign(&f0_step1.slice(s![0 .. t_length]));

    let len_min: usize = ((1000.0 / (FLOOR_F0 * FRAMEPERIOD * 0.5)) + 0.5) as usize;

    for i in 0 .. v_cnt {
        if (u_index[i] - v_index[i] + 1) <= len_min {
            f0_step2.slice_mut(s![v_index[i] ..= u_index[i]]).fill(0.0);
        }
    }

    // 3rd step
    f0_step3[t_length-1] = f0_step0[t_length-1];

    for i in (0 ..= t_length-2).rev() {
        if f0_step3[i+1] < FLOOR_F0 { f0_step3[i] = f0_step0[i]; }
        else if (l_lim * f0_step3[i+1] <= f0_step0[i]) && (f0_step0[i] <= h_lim * f0_step3[i+1]) { f0_step3[i] = f0_step0[i]; }
    }

    // 4th step
    (v_cnt, u_cnt) = (0, 0);

    if FLOOR_F0 < f0_step3[0] { v_index[v_cnt] = 0; v_cnt += 1; }

    for i in 1 .. t_length {
        if (f0_step3[i-1] < FLOOR_F0) && (FLOOR_F0 <= f0_step3[i]) { v_index[v_cnt] = i; v_cnt += 1; }
        else if (FLOOR_F0 <= f0_step3[i-1]) && (f0_step3[i] < FLOOR_F0) { u_index[u_cnt] = i-1; u_cnt += 1; }
    }

    u_index[u_cnt] = t_length-1; u_cnt += 1;

    f0_step4.slice_mut(s![0 .. t_length]).assign(&f0_step3.slice(s![0 .. t_length]));

    for i in 0 .. v_cnt {
        if (u_index[i] - v_index[i] + 1) <= len_min { f0_step4.slice_mut(s![v_index[i] ..= u_index[i]]).fill(0.0); }
    }

    // 5th step
    for i in 0 .. t_length {
        if (FLOOR_F0 <= f0_step2[i]) || (FLOOR_F0 <= f0_step4[i]) { f0_step5[i] = f0_step0[i]; }
    }

    f0.slice_mut(s![0 .. t_length]).assign(&f0_step5.slice(s![0 .. t_length]));

    (f0[0], f0[t_length-1]) = (0.0, 0.0);

    lowpass_f0_IIR(f0, 0.1, t_length);
    lowpass_f0_IIR(f0, 0.2, t_length);
    
}

