import { CreateGPGPU, Color, Points, range } from "../../gpgpu.ts/ts/gpgpu.js";

function msg(text: string){
    console.log(text);
}

export function initKesan(greet, positions, colors){
    msg("計算の初期処理");
    greet("main-tsから今日は。");

    var canvas = document.getElementById("webgl-canvas") as HTMLCanvasElement;
    let mygpgpu = CreateGPGPU(canvas);

/*
    const sz = 40;
    const sz2 = sz * sz * sz;
    const positions = new Float32Array(sz2 * 3);
    const colors    = new Float32Array(sz2 * 4);

    let i1 = 0;
    let i2 = 0;
    let sum1 = 0, sum2 = 0;
    for(let idx of range(sz2)){
        while(true){

            let x = -1.5 + 3 * Math.random();
            let y = -1.5 + 3 * Math.random();
            let z = -1.5 + 3 * Math.random();

            let r = Math.sqrt(x*x + y*y + z*z);
            let theta, phi;
            if(r == 0){

                theta = 0;
                phi   = 0;
            }
            else{

                // z = r cos θ
                theta = Math.acos(z / r);

                // x = r sin θ cos φ
                phi   = Math.acos(x / (r * Math.sin(theta)));
            }

            let f = Math.exp(- r) * Math.cos(theta);
            let R = f * f;
            sum1 += R;
            sum2 += Math.random();

            if(R < 1.0 * Math.random()){
                continue;
            }

            positions[i1    ] = x;
            positions[i1 + 1] = y;
            positions[i1 + 2] = z;
            i1  += 3;

            let cr = R;
            let cg = R; // * theta / Math.PI;
            let cb = R; //* phi / Math.PI;

            colors[i2] = cr;
            colors[i2 + 1] = cg;
            colors[i2 + 2] = cb;
            colors[i2 + 3] = 1;
            i2  += 4;

            break;
        }
    }
    console.log(`平均 R = ${sum1/sz2} random:${sum2/sz2}`);
    */
    /*
    for(let ix of range(sz)){
        ix += Math.random() - 0.5;
        let x = -1 + 2.0 * ix / sz;
        let cr = 1.0 * ix / sz;
        for(let iy of range(sz)){
            iy += Math.random() - 0.5;

            let y = -1 + 2.0 * iy / sz;
            let cg = 1.0 * iy / sz;
            for(let iz of range(sz)){
                iz += Math.random() - 0.5;

                let z = -1 + 2.0 * iz / sz;
                let cb = 1.0 * iz / sz;

            }        
        }
    }
    */

    mygpgpu.startDraw3D([ 
        // new Points(new Float32Array([1.5, -1.3, 0, -1.5, -1.3, 0]), new Float32Array([1,0,0,1, 0,0,1,1]), 5)
        new Points(positions, colors, 1)
    ]);

    setInterval(()=>{

    mygpgpu.drawParam.yRot += Math.PI / 40;
    }, 100);
    mygpgpu.drawParam.z = -2.5;
}