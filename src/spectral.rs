use nalgebra::{base::*, ComplexField};
use rand::prelude::*;
use num::complex::{Complex64, self};
use std::{f64::consts::PI, ops::Mul};
use num::zero;
use wasm_bindgen::prelude::*;
use web_sys;

use rand_distr::{StandardNormal, DistIter};

use crate::utils;


#[derive(Clone, Copy)]
#[wasm_bindgen]
pub struct ULA {
    n_sensor: usize,
    d: f64,
}
#[wasm_bindgen]
impl ULA {
    pub fn new(n_sensor: usize, d: f64) -> Self {
        Self {
            n_sensor,
            d,
        }
    }
}
#[wasm_bindgen(getter_with_clone)]
pub struct Signal {
    ula: ULA,
    amp: DVector<f64>,
    doA: DVector<f64>,
    n_sources: usize,
    samples: usize,
    pert: f64,
    coh: bool,
    data: DMatrix<Complex64>,    
}


impl Signal {
    pub fn generate_gaussian_vector(N: usize) -> RowDVector<Complex64> {
        let r_iter:DistIter<StandardNormal, ThreadRng, f64> = thread_rng().sample_iter(StandardNormal);
        let real_noise: RowDVector<f64>= RowDVector::from_iterator(N, r_iter.take(N));
        let r_iter:DistIter<StandardNormal, ThreadRng, f64> = thread_rng().sample_iter(StandardNormal);
        let complex_noise: RowDVector<f64> = RowDVector::from_iterator(N, r_iter.take(N));
        let iter = real_noise.iter().zip(complex_noise.iter());
        let noise: RowDVector<Complex64> = RowDVector::from_iterator(N, iter.map(|(&x, &y)| 
        Complex64::new(x, y).scale(1_f64/2_f64.sqrt())));
        return noise;

    }
    pub fn generate_gaussian_matrix(N: usize, M: usize) -> DMatrix<Complex64> {
        return DMatrix::from_fn(N, M, |_, _| {

            let real_noise = thread_rng().sample(StandardNormal);
            let complex_noise = thread_rng().sample(StandardNormal);
            Complex64::new(real_noise, complex_noise).scale(1_f64/2_f64.sqrt())
        })  
    }

    pub fn with_data(mut self, data: &DMatrix<Complex64>) -> Self {
        self.data = data.clone();
        self
    }

    pub fn get_data(&self) -> &DMatrix<Complex64> {
        &self.data
    }
    pub fn with_pert(mut self, pert: f64) -> Self {
        self.pert = pert;
        self
    }

    pub fn with_coh(mut self, coh: bool) -> Self {
        self.coh = coh;
        self
    }

    fn mut_amp(&mut self) -> &mut DVector<f64> {
        &mut self.amp
    }

    fn mut_doA(&mut self) -> &mut DVector<f64> {
        &mut self.doA
    }



}
#[wasm_bindgen]
impl Signal {
    pub fn new(ula: ULA, N: usize, amp: &[f64], doA: &[f64]) -> Signal {
        utils::set_panic_hook();
        Self {
            ula: ula,
            amp: DVector::from_row_slice(amp),
            doA: DVector::from_row_slice(doA),
            n_sources: amp.len(),
            samples: N,
            pert: 0.0,
            coh: false,
            data: DMatrix::<Complex64>::zeros(0, 0),

        }
        
    }

    pub fn add_source(&mut self, amp: f64, doa: f64) {
        let m: usize = self.ula.n_sensor;
        let N: usize = self.samples;
        let nrows = self.doA.nrows();
        self.mut_doA().resize_vertically_mut(nrows+1, doa);
        self.mut_amp().resize_vertically_mut(nrows+1, amp);
        let tmp = Signal::generate_gaussian_vector(N);
        let exp_comp = |x: usize| (-PI*Complex64::i()*(doa*PI/180_f64).sin()*(x as f64)).exp();
        let mut A: DVector<Complex64> = DVector::from_iterator(m, (0..m).map(|m| exp_comp(m)));
        if self.pert > 0_f64 {
            A = A+(Signal::generate_gaussian_vector(m).transpose()).scale(self.pert.sqrt());

        }
        if self.coh {
            let r_r: f64 = thread_rng().sample(StandardNormal);
            let r_c: f64 = thread_rng().sample(StandardNormal);
            let nn: Complex64 = Complex64::new(r_r, r_c);
            if !self.data.is_empty() {
                self.data.gemm(nn.scale(amp*2_f64.sqrt()), &A, &tmp, Complex64::new(1.0, 0.0));
               
            } else {
                let n = Signal::generate_gaussian_matrix(m, N).scale(2.0);
                let mut sig = DMatrix::<Complex64>::zeros(m, N);
                sig.gemm(nn.scale(amp*2_f64.sqrt()), &A, &tmp, Complex64::new(1.0, 0.0));
                sig = sig+n;
                self.data = sig;

            }

        } else {
            let nn = Signal::generate_gaussian_vector(N);
            let alpha = Complex64::new(1.0, 0.0).scale(amp);
            if !self.data.is_empty() {
                self.data.gemm(alpha, &A, &nn, Complex64::new(1.0, 0.0));
               
            } else {
                let n = Signal::generate_gaussian_matrix(m, N).scale(2.0);
                let mut sig = DMatrix::<Complex64>::zeros(m,N);
                sig.gemm(alpha, &A, &nn, Complex64::new(1.0, 0.0));
                sig = sig+n;
                self.data = sig;

            }

        }

        
    }
    pub fn with_random_signal(mut self) -> Signal {
        let m: usize = self.ula.n_sensor;
        let N: usize = self.samples;
        let mut signal = DMatrix::zeros(m, N);
        let tmp = Signal::generate_gaussian_vector(N);
        for k in 0..self.n_sources {
            let exp_comp = |x: usize| (-PI*Complex64::i()*(self.doA[k]*PI/180_f64).sin()*(x as f64)).exp();
            let mut A: DVector<Complex64> = DVector::from_iterator(m, (0..m).map(|m| exp_comp(m)));
            if self.pert > 0_f64 {
                A = A+(Signal::generate_gaussian_vector(m).transpose()).scale(self.pert.sqrt());

            }
            if self.coh {
                let r_r: f64 = thread_rng().sample(StandardNormal);
                let r_c: f64 = thread_rng().sample(StandardNormal);
                let nn: Complex64 = Complex64::new(r_r, r_c);
                signal.gemm(nn.scale(self.amp[k]*2_f64.sqrt()), &A, &tmp, Complex64::new(1.0, 0.0));

            } else {
                let nn = Signal::generate_gaussian_vector(N);
                let alpha = Complex64::new(1.0, 0.0).scale(self.amp[k]);
                signal.gemm(alpha, &A, &nn, Complex64::new(1.0, 0.0));

            }
        }
        let n = Signal::generate_gaussian_matrix(m, N).scale(2.0);
        signal = signal+n;
        self.data = signal;
        self

    }

    pub fn beamform(&self, L: usize) -> Vec<f64> {
        if !self.data.is_empty() {
            return super::beamform(&self.data, L, self.ula.d)
        } else {
           //console::error("Signal has no data");
           Vec::new()
        }
    }
}



pub fn beamform(Y: &DMatrix<Complex64>, L: usize, d: f64) -> Vec<f64> {
    let m = Y.nrows();
    let N = Y.ncols();
    let R = (Y*Y.adjoint()).scale(1_f64/N as f64);

    let mut phi = DVector::<f64>::zeros(L);
    for i in 0..L {
        let f = |m: f64| {(-2.0*PI*Complex64::i()*d*(-0.5*PI + PI*(i as f64)/L as f64).sin()*m).exp()};
        let a = DVector::from_fn(m, |i, _| f(i as f64));
        phi[i] = (a.adjoint()*&R*&a)[(0, 0)].re;
        
        
    }

    return phi.as_slice().to_vec();

    //let R = Y.gemm_ad(alpha, &Y, &Y, num::zero());
}




