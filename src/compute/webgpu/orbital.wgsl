struct Term {
    cx: f32, cy: f32, cz: f32, exponent: f32,
    coefficient: f32, ax: u32, ay: u32, az: u32,
};
struct Params {
    min_xyz: vec4<f32>, max_xyz: vec4<f32>,
    dims: vec4<u32>,
    base_xyz_count: vec4<u32>,
};
@group(0) @binding(0) var<storage, read> terms: array<Term>;
@group(0) @binding(1) var<uniform> params: Params;
@group(0) @binding(2) var<storage, read_write> result: array<f32>;

fn power4(x: f32, degree: u32) -> f32 {
    var value = 1.0;
    for (var i = 0u; i < degree; i = i + 1u) { value = value * x; }
    return value;
}

@compute @workgroup_size(64)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    if (i >= params.base_xyz_count.w) { return; }
    let nx = params.dims.x;
    let ny = params.dims.y;
    let sx = params.base_xyz_count.x + i;
    let ix = sx % nx;
    let sy = params.base_xyz_count.y + sx / nx;
    let iy = sy % ny;
    let iz = params.base_xyz_count.z + sy / ny;
    var fx = 0.0;
    var fy = 0.0;
    var fz = 0.0;
    if (nx > 1u) { fx = f32(ix) / f32(nx - 1u); }
    if (ny > 1u) { fy = f32(iy) / f32(ny - 1u); }
    if (params.dims.z > 1u) { fz = f32(iz) / f32(params.dims.z - 1u); }
    let xyz = params.min_xyz.xyz + vec3<f32>(fx, fy, fz) *
              (params.max_xyz.xyz - params.min_xyz.xyz);
    var sum = 0.0;
    for (var t = 0u; t < params.dims.w; t = t + 1u) {
        let term = terms[t];
        let dx = xyz.x - term.cx;
        let dy = xyz.y - term.cy;
        let dz = xyz.z - term.cz;
        sum = sum + term.coefficient * power4(dx, term.ax) *
              power4(dy, term.ay) * power4(dz, term.az) *
              exp(-term.exponent * (dx * dx + dy * dy + dz * dz));
    }
    result[i] = sum;
}
