import { CreateGPGPU, Color, Points } from "../../gpgpu.ts/ts/gpgpu.js";

function msg(text: string){
    console.log(text);
}

export function initKesan(){
    msg("計算の初期処理");

    var canvas = document.getElementById("webgl-canvas") as HTMLCanvasElement;
    let mygpgpu = CreateGPGPU(canvas);

    mygpgpu.startDraw3D([ 
        new Points(new Float32Array([1.5, -1.3, 0, -1.5, -1.3, 0]), Color.red, 5)
    ]);

}