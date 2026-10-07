#include <metal_stdlib>
using namespace metal;

struct CovGaussianTerm {
    float cx, cy, cz, exponent;
    float coefficient;
    uint ax, ay, az;
};
static_assert(sizeof(CovGaussianTerm) == 32);

struct CovGridRequest {
    float min_x, min_y, min_z;
    float max_x, max_y, max_z;
    uint nx, ny, nz;
    ulong first_point;
    uint point_count;
};
static_assert(sizeof(CovGridRequest) == 56);

float power4(float value, uint degree) {
    float result = 1.0f;
    for (uint i = 0; i < degree; ++i) result *= value;
    return result;
}

kernel void cov_orbital(
    device const CovGaussianTerm* terms [[buffer(0)]],
    constant uint& term_count [[buffer(1)]],
    constant CovGridRequest& grid [[buffer(2)]],
    device float* output [[buffer(3)]],
    uint i [[thread_position_in_grid]]) {
    if (i >= grid.point_count) return;
    const ulong index = grid.first_point + ulong(i);
    const uint ix = uint(index % ulong(grid.nx));
    const uint iy = uint((index / ulong(grid.nx)) % ulong(grid.ny));
    const uint iz = uint(index / (ulong(grid.nx) * ulong(grid.ny)));
    const float x = grid.min_x + (grid.max_x - grid.min_x) *
        (grid.nx == 1 ? 0.0f : float(ix) / float(grid.nx - 1));
    const float y = grid.min_y + (grid.max_y - grid.min_y) *
        (grid.ny == 1 ? 0.0f : float(iy) / float(grid.ny - 1));
    const float z = grid.min_z + (grid.max_z - grid.min_z) *
        (grid.nz == 1 ? 0.0f : float(iz) / float(grid.nz - 1));
    float sum = 0.0f;
    for (uint j = 0; j < term_count; ++j) {
        const CovGaussianTerm t = terms[j];
        if (t.coefficient == 0.0f) continue;
        const float dx = x - t.cx, dy = y - t.cy, dz = z - t.cz;
        const float radial = exp(-t.exponent * (dx * dx + dy * dy + dz * dz));
        if (radial == 0.0f) continue;
        sum += t.coefficient * power4(dx, t.ax) * power4(dy, t.ay) *
            power4(dz, t.az) * radial;
    }
    output[i] = sum;
}
