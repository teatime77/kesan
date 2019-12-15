#![allow(non_snake_case)]
#![allow(dead_code)]

extern crate wasm_bindgen;
extern crate num_traits;
extern crate rand;

use web_sys::console;

//--------------------------------------------------------------------------------
//
//  MatD
//
//--------------------------------------------------------------------------------

// #[wasm_bindgen]
pub struct MatD {
    nrow: u32,
    ncol: u32,
    pub dt: Vec<f64>,
}

// #[wasm_bindgen]
impl MatD {
    pub fn new(nrow:u32, ncol:u32) -> MatD {
        let nrow:u32 = nrow;
        let ncol:u32 = ncol;

        let dt: Vec<f64> = vec![0.0; (nrow * ncol) as usize];

        MatD {
            nrow,
            ncol,
            dt,
        }
    }

    pub fn I(dim:u32) -> MatD {
        let mut dt: Vec<f64> = vec![0.0; (dim * dim) as usize];

        for i in 0..dim {
            dt[(i * dim + i) as usize] = 1.0;
        }

        MatD {
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

    pub fn dt(&self) -> *const f64 {
        self.dt.as_ptr()
    }

    pub fn set2(mut self, dt: *const f64){
        unsafe {
            for i in 0..self.dt.len(){
                self.dt[i] = * dt.offset(i as isize);
            }
        }
    }

//    #[inline]
    #[inline(always)]
    fn at(&self, row:u32, col:u32) -> f64 {
        self.dt[(row * self.ncol + col) as usize]
    }

//    #[inline]
    #[inline(always)]
    fn at2(&mut self, row:u32, col:u32) -> &mut f64 {
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


    pub fn abs(&self) -> MatD {
        let dt: Vec<f64> = self.dt.iter().map(|x| x.abs() ).collect();

        MatD {
            nrow: self.nrow,
            ncol: self.ncol,
            dt,
        }
    }

    pub fn max(&self) -> f64 {
        self.dt.iter().fold(0.0, |x, y| x.max(*y))
    }

    pub fn sub(&self, A: &MatD) -> MatD {
        let dt: Vec<f64> = (0..self.dt.len()).map(|i| self.dt[i] - A.dt[i] ).collect();

        MatD {
            nrow: self.nrow,
            ncol: self.ncol,
            dt,
        }
    }


    pub fn dot(&self,  A: &MatD) -> MatD{
//        console.assert(this.cols == A.rows);

        let nrow = self.nrow;
        let ncol = A.ncol;

        let mut dt: Vec<f64> = vec![0.0; (nrow * ncol) as usize];

        for row in 0..nrow {
            for col in 0..ncol {
                let mut sum = 0.0;

                for k in 0..self.ncol {
                    sum += self.at(row, k) * A.at(k, col);
                }

                dt[(row * self.ncol + col) as usize] = sum;
            }
        }

        MatD {
            nrow,
            ncol,
            dt,
        }
    }

    pub fn cat(&self, A: &MatD) -> MatD {
        let mut B = MatD::new(self.nrow, self.ncol + A.ncol);

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
        let mut maxVal = num_traits::abs(self.at(idx, idx) );

        for row in idx + 1..self.nrow {
            let val = num_traits::abs(self.at(row, idx));

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

    pub fn inv(&self) -> MatD {
//    console.assert(this.rows == this.cols);

        // let mut A = self.cat(&MatD::I(self.nrow));
        let mut A = MatD::new(self.nrow, self.ncol + self.ncol);

        for row in 0..self.nrow {
            for col in 0..self.ncol {
                * A.at2(row, col) = self.at(row, col);
            }

            * A.at2(row, self.ncol + row) = 1.0;
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
            let  div = 1.0 / A.at(r1, r1);

            for c in 0..A.ncol {
                if c == r1 {
                    // 対角成分の場合

                    * A.at2(r1,c) = 1.0;
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

                    let a = - A.at(r2, r1);
                    *A.at2(r2, r1) = 0.0;

                    for  c in r1 + 1..A.ncol {
                        *A.at2(r2, c) += a * A.at(r1, c) ;
                    }
                }
            }
        }

        let mut B = MatD::new( self.nrow, self.ncol);

        for  r in 0..self.nrow {
            for c in 0..self.ncol {
                *B.at2(r, c) = A.at(r, self.ncol + c);
            }
        }

        B
    }
}
