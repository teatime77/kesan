#![allow(non_snake_case)]

extern crate wasm_bindgen;
extern crate num_traits;
extern crate rand;

mod mat;

use wasm_bindgen::prelude::*;
use web_sys::console;
use rand::prelude::*;

use self::mat::matf::MatF;
use self::mat::matd::MatD;
use self::mat::matg::Mat;

#[wasm_bindgen]
pub fn greet(name: &str) {
    console::log_1(&format!("Hello, {}!", name).into());

    console::log_1(&"今日はHello using web-sys".into());

    let js: JsValue = 4.into();
    console::log_2(&"Logging arbitrary values looks like".into(), &js);


    let v1: Vec<i32> = vec![1, 2, 3];

    let v2: Vec<_> = v1.iter().map(|x| x.to_string()).collect();
    let joined = v2.join("-");
    console::log_1(&joined.into());
}

#[wasm_bindgen]
pub struct PointsRs {
    size: u32,
    positions: Vec<f32>,
    colors: Vec<f32>,
    avg: f32
}

#[wasm_bindgen]
impl PointsRs {
    pub fn fff(){

    }

    pub fn new(size:u32) -> PointsRs {
        let mut positions: Vec<f32> = vec![0.0; (size * 3) as usize];
        let mut colors   : Vec<f32> = vec![0.0; (size * 4) as usize];

        let mut rng = thread_rng();
    
        let mut i1:usize = 0;
        let mut i2:usize = 0;
        let mut sum1:f32 = 0.0;
        let w:f32 = 16.0;

        for _idx in 0..size {
            loop {
                let x:f32 = -w + 2.0 * w * rng.gen::<f32>();
                let y:f32 = -w + 2.0 * w * rng.gen::<f32>();
                let z:f32 = -w + 2.0 * w * rng.gen::<f32>();
    
                let r:f32 = (x*x + y*y + z*z).sqrt();

                let theta:f32;
                let phi:f32;
                if r == 0.0 {
    
                    continue;
                }
                else{
    
                    // z = r cos θ
                    let z_r = z / r;
                    if z_r < -1.0 || 1.0 < z_r {
                        continue;
                    }
                    theta = z_r.acos();
    
                    // x = r sin θ cos φ
                    let r_sin:f32 = r * theta.sin();
                    if r_sin == 0.0 {

                        continue;
                    }
                    else{

                        let x2 = x / r_sin;
                        if x2 < -1.0 || 1.0 < x2 {
                            continue;
                        }
    
                        phi   = x2.acos();
                    }
                }
    
                let sth:f32 = theta.sin();
                let f:f32 = (- r / 3.0).exp() * sth * sth * (2.0 * phi).sin();

                let R:f32 = f * f;
                sum1 += R;
    
                if R < 1.0 * rng.gen::<f32>() {
                    continue;
                }
    
                positions[i1    ] = x / w;
                positions[i1 + 1] = y / w;
                positions[i1 + 2] = z / w;
                i1  += 3;
    
                let mut cr = 0.0;
                let cg = 0.0;
                let mut cb = 0.0;
                if 0.0 <= x * y {

                    cr = 1.0;
                }
                else{

                    cb = 1.0;
                }
    
                colors[i2] = cr;
                colors[i2 + 1] = cg;
                colors[i2 + 2] = cb;
                colors[i2 + 3] = 1.0;
                i2  += 4;
    
                break;
            }
        }

        let avg:f32 = sum1 / (size as f32);

        PointsRs {
            size,
            positions,
            colors,
            avg
        }
    }

    pub fn size(&self) -> u32 {
        self.size
    }

    pub fn positions(&self) -> *const f32 {
        self.positions.as_ptr()
    }

    pub fn colors(&self) -> *const f32 {
        self.colors.as_ptr()
    }

    pub fn avg(&self) -> f32 {
        self.avg
    }
}



//-------------------------------------------------- Test-Mat

#[wasm_bindgen]
pub struct TestMat {
    diff1: f32,
    diff2: f64,
}

#[wasm_bindgen]
impl TestMat {
    pub fn new(nrow:u32, ncol:u32) -> TestMat {
        let mut rng = thread_rng();

        let mut A1: Mat<f32> = Mat::<f32>::new(nrow, ncol);
        let mut A2: Mat<f64> = Mat::<f64>::new(nrow, ncol);

        for i in 0..(nrow*ncol) {
            A1.dt[i as usize] = rng.gen::<f32>();
            A2.dt[i as usize] = rng.gen::<f64>();
        }

        let B1: Mat<f32> = A1.inv();
        let B2: Mat<f64> = A2.inv();

        let C1: Mat<f32> = A1.dot(&B1);
        let C2: Mat<f64> = A2.dot(&B2);

        let I1 = Mat::<f32>::I(nrow);
        let I2 = Mat::<f64>::I(nrow);

        let D1 = I1.sub(&C1);
        let D2 = I2.sub(&C2);

        let E1 = D1.abs();
        let E2 = D2.abs();

        E1.print("E1 ");
        E2.print("E2 ");

        let diff1 = E1.max();
        let diff2 = E2.max();

        TestMat {
            diff1,
            diff2,
        }
    }

    pub fn test(nrow:u32, ncol:u32) {
        let mut rng = thread_rng();

        let mut A1: MatF = MatF::new(nrow, ncol);
        let mut A2: MatD = MatD::new(nrow, ncol);

        for i in 0..(nrow*ncol) {
            A1.dt[i as usize] = rng.gen::<f32>();
            A2.dt[i as usize] = rng.gen::<f64>();
        }

        let B1: MatF = A1.inv();
        let B2: MatD = A2.inv();

        let C1: MatF = A1.dot(&B1);
        let C2: MatD = A2.dot(&B2);

        let I1 = MatF::I(nrow);
        let I2 = MatD::I(nrow);

        let D1 = I1.sub(&C1);
        let D2 = I2.sub(&C2);

        let E1 = D1.abs();
        let E2 = D2.abs();

        E1.print("E1-F ");
        E2.print("E2-D ");

        let diff1 = E1.max();
        let diff2 = E2.max();

        let js1: JsValue = diff1.into();
        let js2: JsValue = diff2.into();

        console::log_2(&"diff1 = ".into(), &js1);
        console::log_2(&"diff2 = ".into(), &js2);
    }


    pub fn diff1(&self) -> f32 {
        self.diff1
    }

    pub fn diff2(&self) -> f64 {
        self.diff2
    }
}


