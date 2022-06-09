use nalgebra::{base::*, ComplexField, Complex};
use rand::prelude::*;
use num::{complex::{Complex64, self}, Zero, One};
use std::{f64::consts::PI, ops::Mul};
use num::zero;
use wasm_bindgen::prelude::*;
use web_sys;

use rand_distr::{StandardNormal, DistIter};

use crate::utils;

// Struct representing a ULA
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
    pub fn n_sensor(&self) -> usize {
        self.n_sensor
    }
}
// Struct representing some signal(s) impinging on an array of sensors
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
    // Generates a normally distributed complex vector of size N
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
    // Generates a normally distributed complex matrix of size NxM
    pub fn generate_gaussian_matrix(N: usize, M: usize) -> DMatrix<Complex64> {
        return DMatrix::from_fn(N, M, |_, _| {

            let real_noise = thread_rng().sample(StandardNormal);
            let complex_noise = thread_rng().sample(StandardNormal);
            Complex64::new(real_noise, complex_noise).scale(1_f64/2_f64.sqrt())
        })  
    }

    pub fn get_data(&self) -> &DMatrix<Complex64> {
        &self.data
    }


    fn mut_amp(&mut self) -> &mut DVector<f64> {
        &mut self.amp
    }

    fn mut_doA(&mut self) -> &mut DVector<f64> {
        &mut self.doA
    }
    pub fn with_data(mut self, data: &DMatrix<Complex64>) -> Self {
        self.data = data.clone();
        self
    }


}
// WASM bindings
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

    pub fn with_pert(mut self, pert: f64) -> Self {
        self.pert = pert;
        self
    }

    pub fn with_coh(mut self, coh: bool) -> Self {
        self.coh = coh;
        self
    }

    pub fn set_coh(&mut self, coh: bool) {
        self.coh = coh;
    }

    pub fn set_pert(&mut self, pert: f64) {
        self.pert = pert;
    }

    pub fn set_n_sensor(&mut self, n_sensor: usize) {
        self.ula.n_sensor = n_sensor;
    }
    // Adds a source with a random signal with amplitude amp and direction doa.
    pub fn add_source(&mut self, amp: f64, doa: f64) {
        self.n_sources +=1;
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
                let n = Signal::generate_gaussian_matrix(m, N).scale(5.0);
                let mut sig = DMatrix::<Complex64>::zeros(m,N);
                sig.gemm(alpha, &A, &nn, Complex64::new(1.0, 0.0));
                sig = sig+n;
                self.data = sig;

            }

        }

        
    }

    pub fn n_sensor(&self) -> usize {
        self.ula.n_sensor()
    }
    pub fn n_sources(&self) -> usize {
        self.n_sources
    }

    pub fn clear(&mut self) {
        self.data =  DMatrix::<Complex64>::zeros(0, 0);
    }
    // generates random signals from the given signal data;
    pub fn random_signal(&mut self) {
        let m: usize = self.ula.n_sensor;
        let N: usize = self.samples;
        let mut signal = DMatrix::zeros(m, N);
        let tmp = Signal::generate_gaussian_vector(N);
        for k in 0..self.n_sources {
            let exp_comp = |x: usize| (-PI*Complex64::i()*(self.doA[k]*PI/180_f64).sin()*(x as f64)).exp();
            let mut A: DVector<Complex64> = DVector::from_iterator(m, (0..m).map(|m| exp_comp(m)));
            // Perturbation disturbance
            if self.pert > 0_f64 {
                A = A+(Signal::generate_gaussian_vector(m).transpose()).scale(self.pert.sqrt());

            }
            // Coherence disturbance
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
       
    }
    // Builder method
    pub fn with_random_signal(mut self) -> Signal {
        self.random_signal();
        self

    }
    // Classical beamformer
    pub fn beamform(&self, L: usize) -> Vec<f64> {
        if !self.data.is_empty() {
            return super::beamform(&self.data, L, self.ula.d)
        } else {
           //console::error("Signal has no data");
           Vec::new()
        }
    }

    pub fn capon(&self, L: usize) -> Vec<f64> {
        if !self.data.is_empty() {
            return super::capon(&self.data, L, self.ula.d)
        } else {
            Vec::new()
        }
    }

    pub fn s_apes(&self, L: usize) -> Vec<f64> {
        if !self.data.is_empty() {
            return super::s_apes(&self.data, L);
        } else {
            Vec::new()
        }       
    }
}


// Classical beamformer
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
// Capon estimate. Implementation referenced from "Spectral Analysis of Signals" by 
// P.Stoica and R.Moses
pub fn capon(Y: &DMatrix<Complex64>, L: usize, d: f64) -> Vec<f64> {
    let m: usize = Y.nrows();
    let N: usize = Y.ncols();
    let R = (Y*Y.adjoint()).scale(1_f64/N as f64);
    let CH = R.cholesky();
    let mut phi = DVector::<f64>::zeros(L);
    let phi: Vec<f64> = (0..L).map(|i| {
        let f = |m: f64| {(-2.0*PI*Complex64::i()*d*(-0.5*PI + PI*(i as f64)/L as f64).sin()*m).exp()};
        let a = DVector::from_fn(m, |i, _| f(i as f64));
        if let Some(IR) = &CH {
            1_f64/((a.adjoint()*IR.solve(&a))[(0,0)].re)
        } else {
            0.0
        }
        

    }).collect();
    return phi
}

/*
 A. Jakobsson and P. Stoica, "On the Forward-Backward Spatial APES",
% Signal Processing, Vol. 86, pp. 710-715, 2006.
*/
pub fn s_apes(Y: &DMatrix<Complex64>, L: usize) -> Vec<f64> {
    let m: usize = Y.nrows();
    let N: usize = Y.ncols();
    //println!("m  and n is {}, {}", m, N);
    let _J = DMatrix::<Complex64>::from_fn(m-1, m-1, |i, j| {
        if m-2-j == i {
            Complex64::new(1_f64, 0_f64)
        } else {
            Complex64::zero()
        }
    });
    let Rfb = (Y.index((0..m-1, ..))*Y.index((0..m-1, ..)).adjoint()
        + Y.index((1..m, ..))*Y.index((1..m, ..)).adjoint()).scale(1.0/(N as f64*2.0));
    let Rfb = (&Rfb + &_J*Rfb.transpose()*&_J).scale(0.5);

    let phi: Vec<f64> = (0..L).map(|i| {
        let omega_s = PI*(-0.5*PI + PI*(i as f64)/L as f64).sin();
        let mut Gf = DMatrix::<Complex64>::zeros(m-1,m-1);
        let mut Gb = Gf.clone();
        for t in 0..N {
            let gk = Y.index((0..m-1, t)) + Y.index((1..m, t)) * (Complex64::i()*omega_s).exp();
            let Yud_1 = DVector::<Complex64>::from_iterator(m-1, Y.index((0..m-1, t)).iter().cloned().rev());
            let Yud_2 = DVector::<Complex64>::from_iterator(m-1, Y.index((1..m, t)).iter().cloned().rev());
            let gt = Yud_2.conjugate() + Yud_1.conjugate()*(Complex64::i()*omega_s).exp();

            Gf.gerc(Complex64::one(), &gk, &gk, Complex64::one());

            Gb.gerc(Complex64::one(), &gt, &gt, Complex64::one());
        }
        Gb.scale_mut(1.0/(N as f64*4.0));
        Gf.scale_mut(1.0/(N as f64*4.0));
        let a1 = DVector::from_iterator(m-1, (0..m-1).
        map(|x| (-1.0*Complex64::i()*omega_s*x as f64).exp()));

        let iQ = (&Rfb-Gb.scale(0.5)-Gf.scale(0.5)).lu().try_inverse().unwrap();
        let iQa = &iQ*&a1;
        let hfb = &iQa*(1.0/(a1.adjoint()*&iQa)[0]);
        (hfb.adjoint()*Gf*hfb)[0].abs()

    }).collect();
    return phi;
}




