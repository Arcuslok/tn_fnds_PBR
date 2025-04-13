#![allow(unused_assignments)]
use super::constants::{LN_2, NMIN};

use ndrustfft::{ndfft_r2c, ndifft_r2c, Complex, R2cFftHandler, Normalization};
use ndarray::{Array1, concatenate, Axis, s};

pub fn fftfltQ(x: &Array1<f32>, h: &mut Array1<f32>, y: &mut Array1<f32>, i_num: usize, h_num: usize, o_num: usize) {

    let mut n: usize = libm::powf(2.0, ((libm::logf(2.0 * h_num as f32 + 1.0) / LN_2) as usize + 1) as f32) as usize;
    
    if n < NMIN { n = NMIN; }

    let n2: f32 = n as f32;
    let l: usize = n - h_num + 1;
    let frames: usize = i_num / l;

    let mut x_fft: Array1<f32> = Array1::<f32>::zeros(n);
    let mut manipulator: R2cFftHandler<f32> = R2cFftHandler::<f32>::new(n).normalization(Normalization::None);
    let (mut x_spec, mut b_spec): (Array1<Complex<f32>>, Array1<Complex<f32>>) = (Array1::<Complex<f32>>::zeros(n / 2 + 1), Array1::<Complex<f32>>::zeros(n / 2 + 1));

    ndfft_r2c(&concatenate![Axis(0), h.slice(s![0 .. h_num]).to_owned(), Array1::<f32>::zeros(n - h_num)], &mut b_spec, &mut manipulator, 0);

    *y = Array1::<f32>::zeros(o_num);

    let mut offset: usize = 0;

    for _ in 0 .. frames {
        
        for (n, m) in (0..l).zip((offset..).take(l)) { x_fft[n] = x[m] / n2; }

        x_fft.slice_mut(s![l .. n]).fill(0.0);

        ndfft_r2c(&x_fft, &mut x_spec, &mut manipulator, 0);

        x_spec = x_spec * b_spec.to_owned();
        
        ndifft_r2c(&x_spec, &mut x_fft, &mut manipulator, 0);

        for (m, n) in (offset .. ).zip(0 .. n) { y[m] += x_fft[n]; }

        offset += l;

    }

    let mut n_end: usize = n.min(i_num - offset);

    if n_end < 1 { return; }

    for (n, m) in (0 .. n_end).zip(offset ..) { x_fft[n] = x[m] / n2; }

    x_fft.slice_mut(s![n_end .. n]).fill(0.0);

    ndfft_r2c(&x_fft, &mut x_spec, &mut manipulator, 0);

    x_spec = x_spec * b_spec;

    ndifft_r2c(&x_spec, &mut x_fft, &mut manipulator, 0);

    n_end = n.min(o_num - offset);

    for (m, n) in (offset .. ).zip(0 .. n_end) { y[m] += x_fft[n]; }

}
