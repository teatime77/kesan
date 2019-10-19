extern crate wasm_bindgen;

use wasm_bindgen::prelude::*;
use web_sys::console;

#[wasm_bindgen]
pub fn greet(name: &str) {
    console::log_1(&format!("Hello, {}!", name).into());

    console::log_1(&"今日はHello using web-sys".into());

    let js: JsValue = 4.into();
    console::log_2(&"Logging arbitrary values looks like".into(), &js);
}
