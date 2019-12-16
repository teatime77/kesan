extern crate num;

extern crate wasm_bindgen;
extern crate num_traits;
extern crate rand;

use web_sys::console;

use num::traits::Num;
use num::complex::{Complex32, Complex64};

pub trait FromF64<T> {
    fn from_f64(f: f64) -> T;
}

impl FromF64<Complex32> for Complex32 {
    fn from_f64(f: f64) -> Complex32 {
        Complex32::new(f as f32, 0.0)
    }
}

impl FromF64<Complex64> for Complex64 {
    fn from_f64(f: f64) -> Complex64 {
        Complex64::new(f, 0.0)
    }
}

impl FromF64<f32> for f32 {
    fn from_f64(f: f64) -> f32 {
        f as f32
    }
}

impl FromF64<f64> for f64 {
    fn from_f64(f: f64) -> f64 {
        f
    }
}

pub struct MatC<T: Clone + Num + FromF64<T> > {
    nrow: u32,
    ncol: u32,
    pub dt: Vec<T>,
}

// #[wasm_bindgen]
impl<T: Clone + Copy + Num + FromF64<T> + std::string::ToString > MatC<T> {
    pub fn new(nrow:u32, ncol:u32) -> MatC<T> {
        let nrow:u32 = nrow;
        let ncol:u32 = ncol;

        let dt: Vec<T> = vec![ T::from_f64(0.0) ; (nrow * ncol) as usize];

        MatC::<T> {
            nrow,
            ncol,
            dt,
        }
    }

    #[inline(always)]
    fn at(&self, row:u32, col:u32) -> T {
        self.dt[(row * self.ncol + col) as usize]
    }

    pub fn print(&self, s: &str){
        console::log_1(&format!("{} = [", s).into());

        for row in 0..self.nrow {

            let v2: Vec<_> = (0..self.ncol).map(|col| self.at(row, col).to_string()).collect();
            let joined = v2.join(", ");
            console::log_1(&format!("\t{}", joined).into());
        }
        console::log_1(&"]".into());
    }
}

pub fn testC(){
    let m = MatC::<Complex32>::new(5, 5);
    m.print("mat-c");

    let m1 = MatC::<Complex64>::new(5, 5);
    m1.print("mat-cd");

    let m2 = MatC::<f32>::new(5, 5);
    m2.print("mat-32");

    let m3 = MatC::<f64>::new(5, 5);
    m3.print("mat-64");

}