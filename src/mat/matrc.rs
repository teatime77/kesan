#![allow(non_snake_case)]
#![allow(dead_code)]

extern crate wasm_bindgen;
extern crate num_traits;
extern crate rand;
extern crate num;

use web_sys::console;
use num::traits::Num;
use num::complex::{Complex32, Complex64};

pub trait FloatComplex<T> {
    fn from_f64(f: f64) -> T;
    fn max_abs(&self) -> f64;
    fn new_abs(&self) -> T;
    fn get_max(&self, x:T) -> T;
}

impl FloatComplex<f32> for f32 {
    fn from_f64(f: f64) -> f32 {
        f as f32
    }

    fn max_abs(&self) -> f64 {
        self.abs() as f64        
    }

    fn new_abs(&self) -> f32 {
        self.abs()
    }

    fn get_max(&self, x:f32) -> f32 {
        self.max(x)                
    }
}

impl FloatComplex<f64> for f64 {
    fn from_f64(f: f64) -> f64 {
        f
    }

    fn max_abs(&self) -> f64 {
        self.abs()
    }

    fn new_abs(&self) -> f64 {
        self.abs()
    }

    fn get_max(&self, x:f64) -> f64 {
        self.max(x)        
    }
}

impl FloatComplex<Complex32> for Complex32 {
    fn from_f64(f: f64) -> Complex32 {
        Complex32::new(f as f32, 0.0)
    }

    fn max_abs(&self) -> f64 {
        self.re.abs().max(self.im.abs()) as f64
    }

    fn new_abs(&self) -> Complex32 {
        Complex32::new(self.re.abs(), self.im.abs())
    }

    fn get_max(&self, x:Complex32) -> Complex32 {
        if self.max_abs() < x.max_abs() {
            x
        }
        else{
            *self
        }
    }
}

impl FloatComplex<Complex64> for Complex64 {
    fn from_f64(f: f64) -> Complex64 {
        Complex64::new(f, 0.0)
    }

    fn max_abs(&self) -> f64 {
        self.re.abs().max(self.im.abs())
    }

    fn new_abs(&self) -> Complex64 {
        Complex64::new(self.re.abs(), self.im.abs())
    }

    fn get_max(&self, x:Complex64) -> Complex64 {
        if self.max_abs() < x.max_abs() {
            x
        }
        else{
            *self
        }
    }
}

//--------------------------------------------------------------------------------
//
//  MatRC<T>
//
//--------------------------------------------------------------------------------

// #[wasm_bindgen]
pub struct MatRC<T> {
    nrow: u32,
    ncol: u32,
    pub dt: Vec<T>,
}

// #[wasm_bindgen]
impl<T: Copy + Num + FloatComplex<T> + std::string::ToString + std::ops::AddAssign + std::ops::MulAssign> MatRC<T> {
    pub fn new(nrow:u32, ncol:u32) -> MatRC<T> {
        let nrow:u32 = nrow;
        let ncol:u32 = ncol;

        let dt: Vec<T> = vec![ T::from_f64(0.0) ; (nrow * ncol) as usize];

        MatRC::<T> {
            nrow,
            ncol,
            dt,
        }
    }

    pub fn I(dim:u32) -> MatRC<T> {
        let mut dt: Vec<T> = vec![ T::from_f64(0.0) ; (dim * dim) as usize];

        for i in 0..dim {
            dt[(i * dim + i) as usize] = T::from_f64(1.0);
        }

        MatRC::<T> {
            nrow: dim,
            ncol: dim,
            dt,
        }
    }

    pub fn nrow(&self) -> u32 {
        self.nrow
    }

    pub fn ncol(&self) -> u32 {
        self.ncol
    }

    pub fn dt(&self) -> *const T {
        self.dt.as_ptr()
    }

    pub fn set2(mut self, dt: *const T){
        unsafe {
            for i in 0..self.dt.len(){
                self.dt[i] = * dt.offset(i as isize);
            }
        }
    }

//    #[inline]
    #[inline(always)]
    fn at(&self, row:u32, col:u32) -> T {
        self.dt[(row * self.ncol + col) as usize]
    }

//    #[inline]
    #[inline(always)]
    fn at2(&mut self, row:u32, col:u32) -> &mut T {
        &mut self.dt[(row * self.ncol + col) as usize]
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


    pub fn abs(&self) -> MatRC<T> {
        let dt: Vec<T> = self.dt.iter().map(|x| x.new_abs() ).collect();

        MatRC::<T> {
            nrow: self.nrow,
            ncol: self.ncol,
            dt,
        }
    }

    pub fn max(&self) -> T {
        self.dt.iter().fold( T::from_f64(0.0), |x, y| x.get_max(*y))
    }

    pub fn sub(&self, A: &MatRC<T>) -> MatRC<T> {
        let dt: Vec<T> = (0..self.dt.len()).map(|i| self.dt[i] - A.dt[i] ).collect();

        MatRC::<T> {
            nrow: self.nrow,
            ncol: self.ncol,
            dt,
        }
    }


    pub fn dot(&self,  A: &MatRC<T>) -> MatRC<T> {
//        console.assert(this.cols == A.rows);

        let nrow = self.nrow;
        let ncol = A.ncol;

        let mut dt: Vec<T> = vec![T::from_f64(0.0) ; (nrow * ncol) as usize];

        for row in 0..nrow {
            for col in 0..ncol {
                let mut sum:T = T::from_f64(0.0);

                for k in 0..self.ncol {
                    sum += self.at(row, k) * A.at(k, col);
                }

                dt[(row * self.ncol + col) as usize] = sum;
            }
        }

        MatRC::<T> {
            nrow,
            ncol,
            dt,
        }
    }

    pub fn cat(&self, A: &MatRC<T>) -> MatRC<T> {
        let mut B = MatRC::<T>::new(self.nrow, self.ncol + A.ncol);

        for row in 0..self.nrow {
            for col in 0..self.ncol {
//                B.set(row, col, self.at(row, col));
                * B.at2(row, col) = self.at(row, col);
            }

            for col in 0..A.ncol {
//                B.set(row, self.ncol + col, A.at(row, col));
                * B.at2(row, self.ncol + col) = A.at(row, col);
            }
        }

        B
    }

    fn selectPivot(&self, idx: u32) -> u32 {
        let mut maxRow = idx;
        let mut maxVal = self.at(idx, idx).max_abs();

        for row in idx + 1..self.nrow {
            let val = self.at(row, idx).max_abs();

            if maxVal < val {
                maxRow = row;
                maxVal = val;
            }
        }

        maxRow
    }

    pub fn swapRows(&mut self, row1: u32, row2: u32){
        for col in 0..self.ncol {
            let tmp = self.at(row1, col);

            * self.at2(row1, col) = self.at(row2, col);
            * self.at2(row2, col) = tmp;
        }
    }

    pub fn inv(&self) -> MatRC<T> {
//    console.assert(this.rows == this.cols);

        // let mut A = self.cat(&MatRC<T>::I(self.nrow));
        let mut A = MatRC::<T>::new(self.nrow, self.ncol + self.ncol);

        for row in 0..self.nrow {
            for col in 0..self.ncol {
                * A.at2(row, col) = self.at(row, col);
            }

            * A.at2(row, self.ncol + row) = T::from_f64(1.0);
        }


        for r1 in 0..A.nrow {
            if r1 + 1 < A.nrow {

                // ピボットの行を選ぶ。
                let r2 = A.selectPivot(r1);

                if r2 != r1 {
                    // ピボットの行が現在行と違う場合

                    // 行を入れ替える。
                    A.swapRows(r1, r2);
                }
            }

            // ピボットの逆数
            let  div = T::from_f64(1.0) / A.at(r1, r1);

            for c in 0..A.ncol {
                if c == r1 {
                    // 対角成分の場合

                    * A.at2(r1,c) = T::from_f64(1.0);
                }
                else{
                    // 対角成分でない場合

                    // ピボットの逆数をかける。
                    *A.at2(r1, c) *= div;
                }
            }

            for  r2 in 0..self.nrow{
                if r2 != r1 {
                    // ピボットの行でない場合

                    let a = T::from_f64(0.0) - A.at(r2, r1);
                    *A.at2(r2, r1) = T::from_f64(0.0);

                    for  c in r1 + 1..A.ncol {
                        let a2 = A.at(r1, c);
                        *A.at2(r2, c) += a * a2 ;
                    }
                }
            }
        }

        let mut B = MatRC::<T>::new( self.nrow, self.ncol);

        for  r in 0..self.nrow {
            for c in 0..self.ncol {
                *B.at2(r, c) = A.at(r, self.ncol + c);
            }
        }

        B
    }
}
