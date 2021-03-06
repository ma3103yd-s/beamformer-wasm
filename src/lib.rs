
extern crate web_sys;
extern crate wasm_bindgen;

// A macro to provide `println!(..)`-style syntax for `console.log` logging.

use wasm_bindgen::prelude::*;
mod utils;
mod spectral;
use spectral::{ULA, Signal, beamform, capon, s_apes};
use std::io::BufWriter;
use std::fs::File;

// When the `wee_alloc` feature is enabled, use `wee_alloc` as the global
// allocator.
#[cfg(feature = "wee_alloc")]
#[global_allocator]
static ALLOC: wee_alloc::WeeAlloc = wee_alloc::WeeAlloc::INIT;

#[wasm_bindgen]
extern {
    fn alert(s: &str);
}

#[wasm_bindgen]
pub fn greet() {
    alert("Hello, beamformer!");
}

#[cfg(test)]
mod test_beamformer {
    use std::io::Write;

    use nalgebra::{ComplexField, Normed};

    use super::*;

    #[test]
    fn test_random_noise() {
        let var = 1.0;
        let m = 0.0;

        let noise = Signal::generate_gaussian_vector(100);
        let test_m = noise.sum()/100.0;
        let test_var = noise.norm_squared()/100.0 - test_m.norm_squared();
        println!("Noise variance is {}", test_var);
        println!("Noise mean is {}", test_m);
        println!("Variance is {}", noise.variance().norm());

    }

    #[test]
    fn test_signal() {
        // test if generating sources works.
        let f = File::create("y_values.csv").unwrap();
        let mut writer = BufWriter::new(f);
        let sigAmp = [10.0, 20.0];
        let doa = [25.0, 36.0];
        let ula = ULA::new(10, 0.5);
        let test_signal = Signal::new(ula, 64, &sigAmp, &doa).with_random_signal();
        //println!("Dimension of test_signal is {:?}", test_signal.get_data().shape());
        for row in test_signal.get_data().row_iter() {
            let mut index = 1;
            for &val in row.iter() {
                let real = val.re;
                let c = val.im;
                let str = format!("{} {},", real, c);
                writer.write(str.as_bytes());
                index+=1;
            }
            writer.write("\n".as_bytes());
        }
    }
    #[test]
    fn test_add_source() {
        // Test if adding a source works
        let f = File::create("test_values.txt").unwrap();
        let mut writer = BufWriter::new(f);
        let sigAmp = [];
        let doa = [];
        let ula = ULA::new(10, 0.5);
        let mut test_signal = Signal::new(ula, 64, &sigAmp, &doa);
        test_signal.add_source(10.0, 25.0);
        test_signal.add_source(20.0, 36.0);
        println!("Dimension of test_signal is {:?}", test_signal.get_data().shape());
        for row in test_signal.get_data().row_iter() {
            let mut index = 1;
            for &val in row.iter() {
                println!("Col {}", index);
                let real = val.re;
                let c = val.im;
                let str = format!("{} {},", real, c);
                writer.write(str.as_bytes());
                index+=1;
            }
            writer.write("\n".as_bytes());
        }
    }
    #[test]
    fn test_beamform() {
        // Writes data to file to compare with MATLAB implementation
        let f = File::create("beamform.csv").unwrap();
        let mut writer = BufWriter::new(f);
        let sigAmp = [10.0, 20.0];
        let doa = [25.0, 36.0];
        let ula = ULA::new(50, 0.5);
        let test_signal = Signal::new(ula, 64, &sigAmp, &doa).with_random_signal();
        let phi = test_signal.beamform(512);
        for &val in phi.iter() {
            writer.write(format!("{}, ", val).as_bytes());
        }


    }
    #[test]
    fn test_capon() {
        // Writes data to file to compare with MATLAB implementation
        let f = File::create("capon.csv").unwrap();
        let mut writer = BufWriter::new(f);
        let sigAmp = [10.0, 20.0];
        let doa = [25.0, 36.0];
        let ula = ULA::new(50, 0.5);
        let test_signal = Signal::new(ula, 64, &sigAmp, &doa).with_random_signal();
        let phi = test_signal.capon(512);
        for &val in phi.iter() {
            writer.write(format!("{}, ", val).as_bytes());
        }
    }
    #[test]
    fn test_sapes() {
        // Writes useful data to files to compare with MATLAB implementation
        let f = File::create("sapes.csv").unwrap();
        let mut writer = BufWriter::new(f);
        let sigAmp = [10.0, 20.0];
        let doa = [25.0, 36.0];
        let ula = ULA::new(10, 0.5);
        let test_signal = Signal::new(ula, 64, &sigAmp, &doa).with_random_signal();
        {
            let f2 = File::create("yvalues.csv").unwrap();
            let mut y_writer = BufWriter::new(f2);
            for row in test_signal.get_data().row_iter() {
                let mut index = 1;
                for &val in row.iter() {
                    let real = val.re;
                    let c = val.im;
                    let str = format!("{} {},", real, c);
                    y_writer.write(str.as_bytes());
                    index+=1;
                }
                y_writer.write("\n".as_bytes());
            }
        }

        let phi = test_signal.s_apes(512);


        for &val in phi.iter() {
            writer.write(format!("{}, ", val).as_bytes());
        }
    }

    }
    
