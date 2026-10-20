#include "cov/volume_renderer.hpp"
#include "cov/viewer_layout.hpp"

#include "cov/gl_api.hpp"


#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#ifndef GL_TEXTURE_3D
#define GL_TEXTURE_3D 0x806F
#endif
#ifndef GL_TEXTURE_WRAP_R
#define GL_TEXTURE_WRAP_R 0x8072
#endif
#ifndef GL_CLAMP_TO_EDGE
#define GL_CLAMP_TO_EDGE 0x812F
#endif
#ifndef GL_R32F
#define GL_R32F 0x822E
#endif
#ifndef GL_RED
#define GL_RED 0x1903
#endif
#ifndef GL_TEXTURE0
#define GL_TEXTURE0 0x84C0
#endif
#ifndef GL_VERTEX_SHADER
#define GL_VERTEX_SHADER 0x8B31
#endif
#ifndef GL_FRAGMENT_SHADER
#define GL_FRAGMENT_SHADER 0x8B30
#endif
#ifndef GL_COMPILE_STATUS
#define GL_COMPILE_STATUS 0x8B81
#endif
#ifndef GL_LINK_STATUS
#define GL_LINK_STATUS 0x8B82
#endif
#ifndef GL_INFO_LOG_LENGTH
#define GL_INFO_LOG_LENGTH 0x8B84
#endif

namespace cov {
namespace {

constexpr float kPi = 3.14159265358979323846f;
constexpr float kNearPlane = 0.03f;
constexpr float kFarPlane = 8.0f;

struct Vec3 {
    float x, y, z;
};

Vec3 operator+(const Vec3 a, const Vec3 b) { return {a.x+b.x,a.y+b.y,a.z+b.z}; }
Vec3 operator-(const Vec3 a, const Vec3 b) { return {a.x-b.x,a.y-b.y,a.z-b.z}; }
Vec3 operator*(const Vec3 a, const float s) { return {a.x*s,a.y*s,a.z*s}; }
Vec3 operator/(const Vec3 a, const float s) { return {a.x/s,a.y/s,a.z/s}; }

float dot(const Vec3 a, const Vec3 b) { return a.x*b.x + a.y*b.y + a.z*b.z; }

Vec3 cross(const Vec3 a, const Vec3 b) {
    return {a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x};
}

float length(const Vec3 v) { return std::sqrt(std::max(0.0f, dot(v, v))); }

Vec3 normalise(const Vec3 v) {
    const float n = std::sqrt(std::max(1.0e-20f, dot(v,v)));
    return v * (1.0f/n);
}

struct CameraBasis {
    Vec3 position;
    Vec3 forward;
    Vec3 right;
    Vec3 up;
};

CameraBasis camera_basis(const OrbitCamera& camera) {
    const Vec3 centre{0.5f, 0.5f, 0.5f};
    const float cp = std::cos(camera.pitch);
    const Vec3 offset{
        camera.distance * cp * std::cos(camera.yaw),
        camera.distance * std::sin(camera.pitch),
        camera.distance * cp * std::sin(camera.yaw)
    };
    CameraBasis b{};
    b.position = centre + offset;
    b.forward = normalise(centre - b.position);

    Vec3 world_up{0.0f, 1.0f, 0.0f};
    if (std::abs(dot(b.forward, world_up)) > 0.98f) world_up = {0.0f, 0.0f, 1.0f};
    b.right = normalise(cross(b.forward, world_up));
    b.up = normalise(cross(b.right, b.forward));
    return b;
}

GLuint compile_shader(const GLenum type, const char* source) {
    const GLuint shader = gl::CreateShader(type);
    gl::ShaderSource(shader, 1, &source, nullptr);
    gl::CompileShader(shader);

    GLint ok = 0;
    gl::GetShaderiv(shader, GL_COMPILE_STATUS, &ok);
    if (!ok) {
        GLint len = 0;
        gl::GetShaderiv(shader, GL_INFO_LOG_LENGTH, &len);
        std::vector<char> log(static_cast<std::size_t>(std::max(1, len)), '\0');
        gl::GetShaderInfoLog(shader, len, nullptr, log.data());
        const std::string message(log.data());
        gl::DeleteShader(shader);
        throw std::runtime_error("OpenGL shader compilation failed: " + message);
    }
    return shader;
}

GLuint make_program() {
    static constexpr const char* kVertex = R"GLSL(
#version 120
varying vec2 vUV;
void main() {
    gl_Position = gl_Vertex;
    vUV = gl_MultiTexCoord0.xy;
}
)GLSL";

    static constexpr const char* kFragment = R"GLSL(
#version 120
uniform sampler3D uVolume;
uniform float uIso;
uniform float uOpacity;
uniform int uMaterial;
uniform int uSurfaceMode;
uniform int uPhasePalette;
uniform vec3 uCameraPos;
uniform vec3 uCameraForward;
uniform vec3 uCameraRight;
uniform vec3 uCameraUp;
uniform float uAspect;
uniform float uTanHalfFov;
uniform vec3 uTexel;
varying vec2 vUV;

const float kNear = 0.03;
const float kFar = 8.0;

bool rayBox(vec3 ro, vec3 rd, out float tNear, out float tFar) {
    vec3 inv = 1.0 / rd;
    vec3 t0 = (vec3(0.0) - ro) * inv;
    vec3 t1 = (vec3(1.0) - ro) * inv;
    vec3 lo = min(t0, t1);
    vec3 hi = max(t0, t1);
    tNear = max(max(lo.x, lo.y), lo.z);
    tFar = min(min(hi.x, hi.y), hi.z);
    return tFar >= max(tNear, 0.0);
}

float field(vec3 p) {
    return texture3D(uVolume, clamp(p, vec3(0.0), vec3(1.0))).r;
}

vec3 gradient(vec3 p) {
    float gx = field(p + vec3(uTexel.x,0.0,0.0)) - field(p - vec3(uTexel.x,0.0,0.0));
    float gy = field(p + vec3(0.0,uTexel.y,0.0)) - field(p - vec3(0.0,uTexel.y,0.0));
    float gz = field(p + vec3(0.0,0.0,uTexel.z)) - field(p - vec3(0.0,0.0,uTexel.z));
    return normalize(vec3(gx,gy,gz));
}

float depthFromViewDistance(float viewDepth) {
    float z = max(viewDepth, kNear + 0.0001);
    float ndc = (kFar + kNear) / (kFar - kNear)
              - (2.0 * kFar * kNear) / ((kFar - kNear) * z);
    return clamp(0.5 * ndc + 0.5, 0.0, 1.0);
}

float edgeDistance(float x) {
    float f = fract(x);
    return min(f, 1.0 - f);
}

float triangularWireMask(vec3 p, vec3 n) {
    vec3 an = abs(n);
    vec2 uv;
    if (an.x >= an.y && an.x >= an.z) {
        uv = p.yz;
    } else if (an.y >= an.z) {
        uv = p.xz;
    } else {
        uv = p.xy;
    }

    uv *= 22.0;
    float d0 = edgeDistance(uv.x);
    float d1 = edgeDistance(uv.y);
    float d2 = edgeDistance(uv.x + uv.y);
    float d = min(d0, min(d1, d2));
    return 1.0 - smoothstep(0.035, 0.085, d);
}

void main() {
    vec2 ndc = vUV * 2.0 - 1.0;
    vec3 ro = uCameraPos;
    vec3 rd = normalize(
        uCameraForward +
        ndc.x * uAspect * uTanHalfFov * uCameraRight +
        ndc.y * uTanHalfFov * uCameraUp);

    float t0, t1;
    if (!rayBox(ro, rd, t0, t1)) discard;

    float startT = max(t0, 0.0);
    float travel = max(0.0001, t1 - startT);
    float dt = travel / 512.0;
    float t = startT;
    vec3 p = ro + rd * t;
    float prev = field(p);
    float prevF = abs(prev) - uIso;

    for (int i = 0; i < 512; ++i) {
        t += dt;
        if (t > t1) break;
        p = ro + rd * t;
        float cur = field(p);
        float curF = abs(cur) - uIso;

        if (prevF < 0.0 && curF >= 0.0) {
            vec3 hit = p;
            vec3 n = gradient(hit);
            if (dot(n, rd) > 0.0) n = -n;

            vec3 positive = uPhasePalette==0 ? vec3(0.92, 0.20, 0.17) : vec3(0.12,0.85,0.65);
            vec3 negative = uPhasePalette==0 ? vec3(0.16, 0.38, 0.95) : vec3(0.98,0.66,0.15);
            vec3 base = cur >= 0.0 ? positive : negative;

            // Camera-relative soft key + weak fill. It follows the view so the
            // orbital keeps readable form while rotating, without a harsh studio light.
            vec3 key = normalize(-0.34 * uCameraRight + 0.56 * uCameraUp - 0.75 * uCameraForward);
            vec3 fill = normalize(0.62 * uCameraRight + 0.16 * uCameraUp - 0.34 * uCameraForward);
            vec3 viewDir = normalize(-rd);
            vec3 halfVector = normalize(key + viewDir);
            float keyDiffuse = max(0.0, dot(n, key));
            float fillDiffuse = max(0.0, dot(n, fill));
            float specular = pow(max(0.0, dot(n, halfVector)), 28.0);
            float illumination = 0.42 + 0.40 * keyDiffuse + 0.10 * fillDiffuse;
            vec3 shaded = base * illumination + vec3(0.13 * specular);

            float wire = triangularWireMask(hit, n);
            if (uSurfaceMode == 1 && wire < 0.45) discard;
            if (uSurfaceMode == 1) {
                shaded = mix(base * 0.32, base * 1.04 + vec3(0.10), wire);
            } else if (uSurfaceMode == 2) {
                shaded = mix(shaded, base * 0.30 + vec3(0.16), 0.80 * wire);
            }

            float alpha = clamp(uOpacity, 0.02, 1.0);
            if (uMaterial == 1) {
                float facing = clamp(abs(dot(n, viewDir)), 0.0, 1.0);
                float fresnel = pow(1.0 - facing, 2.2);
                float glassSpec = pow(max(0.0, dot(n, halfVector)), 52.0);
                vec3 glassBody = mix(base * 0.48 + vec3(0.06),
                                     base * 0.88 + vec3(0.12), fresnel);
                shaded = mix(glassBody, vec3(1.0), 0.20 * glassSpec + 0.12 * fresnel);
                alpha *= 0.28 + 0.56 * fresnel;
                alpha = clamp(alpha, 0.02, 0.90);
            }

            if (uSurfaceMode == 1) {
                alpha *= 0.82 + 0.18 * wire;
            }

            float viewDepth = dot(hit - uCameraPos, uCameraForward);
            gl_FragDepth = depthFromViewDistance(viewDepth);
            gl_FragColor = vec4(clamp(shaded, vec3(0.0), vec3(1.0)), alpha);
            return;
        }
        prev = cur;
        prevF = curF;
    }

    discard;
}
)GLSL";

    const GLuint vs = compile_shader(GL_VERTEX_SHADER, kVertex);
    const GLuint fs = compile_shader(GL_FRAGMENT_SHADER, kFragment);
    const GLuint program = gl::CreateProgram();
    gl::AttachShader(program, vs);
    gl::AttachShader(program, fs);
    gl::LinkProgram(program);

    GLint ok = 0;
    gl::GetProgramiv(program, GL_LINK_STATUS, &ok);
    if (!ok) {
        GLint len = 0;
        gl::GetProgramiv(program, GL_INFO_LOG_LENGTH, &len);
        std::vector<char> log(static_cast<std::size_t>(std::max(1, len)), '\0');
        gl::GetProgramInfoLog(program, len, nullptr, log.data());
        const std::string message(log.data());
        gl::DeleteShader(vs);
        gl::DeleteShader(fs);
        gl::DeleteProgram(program);
        throw std::runtime_error("OpenGL program link failed: " + message);
    }

    gl::DeleteShader(vs);
    gl::DeleteShader(fs);
    return program;
}

Vec3 atom_colour(const int z) {
    switch (z) {
        case 1: return {0.95f,0.95f,0.95f};
        case 5: return {0.93f,0.62f,0.62f};
        case 6: return {0.30f,0.31f,0.34f};
        case 7: return {0.20f,0.35f,0.95f};
        case 8: return {0.92f,0.12f,0.12f};
        case 9: return {0.20f,0.85f,0.30f};
        case 14: return {0.62f,0.64f,0.67f};
        case 15: return {0.95f,0.55f,0.10f};
        case 16: return {0.95f,0.85f,0.12f};
        case 17: return {0.15f,0.80f,0.22f};
        case 35: return {0.65f,0.16f,0.16f};
        case 53: return {0.45f,0.18f,0.62f};
        default: return {0.65f,0.65f,0.70f};
    }
}

Vec3 to_texture_point(const double x, const double y, const double z, const GridBox& box) {
    const float dx = std::max(1.0e-8f, box.max_x - box.min_x);
    const float dy = std::max(1.0e-8f, box.max_y - box.min_y);
    const float dz = std::max(1.0e-8f, box.max_z - box.min_z);
    return {
        (static_cast<float>(x) - box.min_x) / dx,
        (static_cast<float>(y) - box.min_y) / dy,
        (static_cast<float>(z) - box.min_z) / dz
    };
}

Vec3 to_texture(const Atom& atom, const GridBox& box) {
    return to_texture_point(atom.x, atom.y, atom.z, box);
}

Vec3 shaded_colour(const Vec3 base,
                   const Vec3 normal,
                   const Vec3 point,
                   const CameraBasis& camera) {
    const Vec3 light = normalise(camera.right * -0.35f + camera.up * 0.58f + camera.forward * -0.74f);
    const Vec3 view = normalise(camera.position - point);
    const Vec3 half_vector = normalise(light + view);
    const float diffuse = std::max(0.0f, dot(normal, light));
    const float specular = std::pow(std::max(0.0f, dot(normal, half_vector)), 30.0f);
    const float illumination = 0.28f + 0.72f * diffuse;
    return {
        std::clamp(base.x * illumination + 0.34f * specular, 0.0f, 1.0f),
        std::clamp(base.y * illumination + 0.34f * specular, 0.0f, 1.0f),
        std::clamp(base.z * illumination + 0.34f * specular, 0.0f, 1.0f)
    };
}

void load_3d_camera(const CameraBasis& b,
                    const float aspect,
                    const float tan_half_fov) {
    const double top = static_cast<double>(kNearPlane) * static_cast<double>(tan_half_fov);
    const double right = top * static_cast<double>(std::max(0.01f, aspect));

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glFrustum(-right, right, -top, top, kNearPlane, kFarPlane);

    const GLfloat view[16] = {
        b.right.x, b.up.x, -b.forward.x, 0.0f,
        b.right.y, b.up.y, -b.forward.y, 0.0f,
        b.right.z, b.up.z, -b.forward.z, 0.0f,
        -dot(b.right, b.position), -dot(b.up, b.position), dot(b.forward, b.position), 1.0f
    };
    glMatrixMode(GL_MODELVIEW);
    glLoadMatrixf(view);
}

void draw_sphere(const Vec3 centre,
                 const float radius,
                 const Vec3 base_colour,
                 const float alpha,
                 const CameraBasis& camera) {
    constexpr int latitude_segments = 18;
    constexpr int longitude_segments = 28;

    for (int lat = 0; lat < latitude_segments; ++lat) {
        const float phi0 = -0.5f * kPi + kPi * static_cast<float>(lat) /
                                             static_cast<float>(latitude_segments);
        const float phi1 = -0.5f * kPi + kPi * static_cast<float>(lat + 1) /
                                             static_cast<float>(latitude_segments);
        const float y0 = std::sin(phi0);
        const float y1 = std::sin(phi1);
        const float r0 = std::cos(phi0);
        const float r1 = std::cos(phi1);

        glBegin(GL_TRIANGLE_STRIP);
        for (int lon = 0; lon <= longitude_segments; ++lon) {
            const float theta = 2.0f * kPi * static_cast<float>(lon) /
                                                static_cast<float>(longitude_segments);
            const float ct = std::cos(theta);
            const float st = std::sin(theta);

            const Vec3 n0{r0 * ct, y0, r0 * st};
            const Vec3 p0 = centre + n0 * radius;
            const Vec3 c0 = shaded_colour(base_colour, n0, p0, camera);
            glColor4f(c0.x, c0.y, c0.z, alpha);
            glVertex3f(p0.x, p0.y, p0.z);

            const Vec3 n1{r1 * ct, y1, r1 * st};
            const Vec3 p1 = centre + n1 * radius;
            const Vec3 c1 = shaded_colour(base_colour, n1, p1, camera);
            glColor4f(c1.x, c1.y, c1.z, alpha);
            glVertex3f(p1.x, p1.y, p1.z);
        }
        glEnd();
    }
}

void draw_cylinder(const Vec3 a,
                   const Vec3 b,
                   const float radius,
                   const Vec3 base_colour,
                   const float alpha,
                   const CameraBasis& camera) {
    const Vec3 axis_vector = b - a;
    const float axis_length = std::sqrt(std::max(1.0e-20f, dot(axis_vector, axis_vector)));
    if (axis_length <= 1.0e-6f) return;
    const Vec3 axis = axis_vector * (1.0f / axis_length);
    const Vec3 helper = std::abs(axis.y) < 0.90f ? Vec3{0.0f, 1.0f, 0.0f}
                                                  : Vec3{1.0f, 0.0f, 0.0f};
    const Vec3 u = normalise(cross(axis, helper));
    const Vec3 v = normalise(cross(axis, u));
    constexpr int slices = 18;

    glBegin(GL_TRIANGLE_STRIP);
    for (int i = 0; i <= slices; ++i) {
        const float theta = 2.0f * kPi * static_cast<float>(i) / static_cast<float>(slices);
        const Vec3 normal = u * std::cos(theta) + v * std::sin(theta);
        const Vec3 pa = a + normal * radius;
        const Vec3 pb = b + normal * radius;
        const Vec3 ca = shaded_colour(base_colour, normal, pa, camera);
        const Vec3 cb = shaded_colour(base_colour, normal, pb, camera);
        glColor4f(ca.x, ca.y, ca.z, alpha);
        glVertex3f(pa.x, pa.y, pa.z);
        glColor4f(cb.x, cb.y, cb.z, alpha);
        glVertex3f(pb.x, pb.y, pb.z);
    }
    glEnd();

    const Vec3 cap_a_colour = shaded_colour(base_colour, axis * -1.0f, a, camera);
    glColor4f(cap_a_colour.x, cap_a_colour.y, cap_a_colour.z, alpha);
    glBegin(GL_TRIANGLE_FAN);
    glVertex3f(a.x, a.y, a.z);
    for (int i = 0; i <= slices; ++i) {
        const float theta = -2.0f * kPi * static_cast<float>(i) / static_cast<float>(slices);
        const Vec3 rim = a + (u * std::cos(theta) + v * std::sin(theta)) * radius;
        glVertex3f(rim.x, rim.y, rim.z);
    }
    glEnd();

    const Vec3 cap_b_colour = shaded_colour(base_colour, axis, b, camera);
    glColor4f(cap_b_colour.x, cap_b_colour.y, cap_b_colour.z, alpha);
    glBegin(GL_TRIANGLE_FAN);
    glVertex3f(b.x, b.y, b.z);
    for (int i = 0; i <= slices; ++i) {
        const float theta = 2.0f * kPi * static_cast<float>(i) / static_cast<float>(slices);
        const Vec3 rim = b + (u * std::cos(theta) + v * std::sin(theta)) * radius;
        glVertex3f(rim.x, rim.y, rim.z);
    }
    glEnd();
}

void draw_dashed_cylinder(const Vec3 a,
                          const Vec3 b,
                          const float radius,
                          const Vec3 colour,
                          const float alpha,
                          const CameraBasis& camera,
                          const int dash_count = 8,
                          const float dash_fraction = 0.55f) {
    if (dash_count <= 0 || dash_fraction <= 0.0f) return;
    const Vec3 delta = b - a;
    for (int i = 0; i < dash_count; ++i) {
        const float t0 = static_cast<float>(i) / static_cast<float>(dash_count);
        const float t1 = std::min(1.0f, t0 + dash_fraction / static_cast<float>(dash_count));
        draw_cylinder(a + delta * t0, a + delta * t1,
                      radius, colour, alpha, camera);
    }
}

std::pair<std::size_t, std::size_t> ordered_atom_pair(const std::size_t a,
                                                       const std::size_t b) {
    return a < b ? std::pair{a, b} : std::pair{b, a};
}

Vec3 heavy_atom_centroid(const Wavefunction& wavefunction,
                         const std::vector<Vec3>& points) {
    Vec3 sum{0.0f, 0.0f, 0.0f};
    std::size_t count = 0;
    for (std::size_t i = 0; i < wavefunction.atoms.size(); ++i) {
        if (wavefunction.atoms[i].atomic_number == 1) continue;
        sum = sum + points[i];
        ++count;
    }
    if (count == 0) {
        for (const Vec3 p : points) sum = sum + p;
        count = points.size();
    }
    return count ? sum / static_cast<float>(count) : Vec3{0.5f, 0.5f, 0.5f};
}

Vec3 inward_delocalisation_offset(const Vec3 a,
                                  const Vec3 b,
                                  const Vec3 centre,
                                  const float distance) {
    const Vec3 axis = normalise(b - a);
    const Vec3 midpoint = (a + b) * 0.5f;
    Vec3 inward = centre - midpoint;
    inward = inward - axis * dot(inward, axis);
    if (length(inward) <= 1.0e-6f) return {0.0f, 0.0f, 0.0f};
    return normalise(inward) * distance;
}

} // namespace

VolumeRenderer::VolumeRenderer() {
    if (!gl::load()) {
        throw std::runtime_error(
            "Required OpenGL 2.1 shader/3D-texture entry points are unavailable");
    }
    program_ = make_program();

    glGenTextures(1, &texture_);
    glBindTexture(GL_TEXTURE_3D, texture_);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
    glBindTexture(GL_TEXTURE_3D, 0);
}

VolumeRenderer::~VolumeRenderer() {
    if (program_) gl::DeleteProgram(program_);
    if (texture_) glDeleteTextures(1, &texture_);
}

void VolumeRenderer::invalidate_geometry_cache() noexcept {
    geometry_cache_wavefunction_ = nullptr;
    geometry_graph_source_ = nullptr;
    geometry_bonds_.clear();
    geometry_bond_indices_.clear();
    geometry_interactions_ = {};
}

void VolumeRenderer::resize_volume(const int nx, const int ny, const int nz) {
    if (nx <= 0 || ny <= 0 || nz <= 0) {
        throw std::runtime_error("Invalid 3D texture dimensions");
    }
    GLint maximum = 0;
    glGetIntegerv(0x8073 /* GL_MAX_3D_TEXTURE_SIZE */, &maximum);
    if (nx > maximum || ny > maximum || nz > maximum)
        throw std::runtime_error("The requested grid exceeds the display texture limit");
    glBindTexture(GL_TEXTURE_3D, texture_);
    while (glGetError() != GL_NO_ERROR) {}
    gl::TexImage3D(GL_TEXTURE_3D, 0, gl::volume_internal_format(),
                   nx, ny, nz, 0, gl::volume_external_format(), GL_FLOAT, nullptr);
    const auto error = glGetError();
    glBindTexture(GL_TEXTURE_3D, 0);
    if (error != GL_NO_ERROR) throw std::runtime_error("Unable to allocate the orbital display texture");
    nx_ = nx;
    ny_ = ny;
    nz_ = nz;
}

void VolumeRenderer::render_volume(const int framebuffer_width,
                                   const int framebuffer_height,
                                   const float isovalue,
                                   const OrbitCamera& camera,
                                   const float opacity,
                                   const OrbitalMaterial material,
                                   const OrbitalSurfaceMode surface_mode,
                                   const int phase_palette) {
    if (nx_ <= 0 || ny_ <= 0 || nz_ <= 0) return;

    const CameraBasis b = camera_basis(camera);
    const auto projection = scene_projection(framebuffer_width, framebuffer_height, camera.fov_degrees);
    const float aspect = projection.aspect;
    const float tan_half_fov = projection.tan_half_vertical_fov;

    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);
    glDepthMask(GL_TRUE);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    gl::UseProgram(program_);

    gl::Uniform1i(gl::GetUniformLocation(program_, "uVolume"), 0);
    gl::Uniform1i(gl::GetUniformLocation(program_, "uPhasePalette"), phase_palette);
    gl::Uniform1f(gl::GetUniformLocation(program_, "uIso"), isovalue);
    gl::Uniform1f(gl::GetUniformLocation(program_, "uOpacity"),
                  std::clamp(opacity, 0.02f, 1.0f));
    gl::Uniform1i(gl::GetUniformLocation(program_, "uMaterial"),
                  static_cast<int>(material));
    gl::Uniform1i(gl::GetUniformLocation(program_, "uSurfaceMode"),
                  static_cast<int>(surface_mode));
    gl::Uniform3f(gl::GetUniformLocation(program_, "uCameraPos"),
                  b.position.x, b.position.y, b.position.z);
    gl::Uniform3f(gl::GetUniformLocation(program_, "uCameraForward"),
                  b.forward.x, b.forward.y, b.forward.z);
    gl::Uniform3f(gl::GetUniformLocation(program_, "uCameraRight"),
                  b.right.x, b.right.y, b.right.z);
    gl::Uniform3f(gl::GetUniformLocation(program_, "uCameraUp"),
                  b.up.x, b.up.y, b.up.z);
    gl::Uniform1f(gl::GetUniformLocation(program_, "uAspect"), aspect);
    gl::Uniform1f(gl::GetUniformLocation(program_, "uTanHalfFov"), tan_half_fov);
    gl::Uniform3f(gl::GetUniformLocation(program_, "uTexel"),
                  1.0f / static_cast<float>(nx_),
                  1.0f / static_cast<float>(ny_),
                  1.0f / static_cast<float>(nz_));

    gl::ActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_3D, texture_);

    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();

    glBegin(GL_QUADS);
    glTexCoord2f(0.0f, 0.0f); glVertex2f(-1.0f, -1.0f);
    glTexCoord2f(1.0f, 0.0f); glVertex2f( 1.0f, -1.0f);
    glTexCoord2f(1.0f, 1.0f); glVertex2f( 1.0f,  1.0f);
    glTexCoord2f(0.0f, 1.0f); glVertex2f(-1.0f,  1.0f);
    glEnd();

    glPopMatrix();
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);

    glBindTexture(GL_TEXTURE_3D, 0);
    gl::UseProgram(0);
    glDisable(GL_BLEND);
}

void VolumeRenderer::render_geometry(const Wavefunction& wavefunction,
                                     const GridBox& box,
                                     const int framebuffer_width,
                                     const int framebuffer_height,
                                     const OrbitCamera& camera,
                                     const MoleculeRenderSettings& settings,
                                     const MoleculeOverlay* overlay,
                                     const InteractionGraph* routed_graph) {
    geometry_targets_.clear();
    geometry_width_=std::max(1,framebuffer_width);
    geometry_height_=std::max(1,framebuffer_height);
    if (wavefunction.atoms.empty()) return;

    const CameraBasis b = camera_basis(camera);
    const auto projection = scene_projection(framebuffer_width, framebuffer_height, camera.fov_degrees);
    const float aspect = projection.aspect;
    const float tan_half_fov = projection.tan_half_vertical_fov;

    std::vector<Vec3> points;
    points.reserve(wavefunction.atoms.size());
    for (const Atom& atom : wavefunction.atoms) points.push_back(to_texture(atom, box));
    if (geometry_cache_wavefunction_ != &wavefunction ||
        geometry_graph_source_ != routed_graph) {
        auto bonds = analyse_bonds(wavefunction);
        auto routed_interactions = routed_graph?*routed_graph:build_interaction_graph(wavefunction);
        geometry_bonds_ = std::move(bonds);
        geometry_interactions_ = std::move(routed_interactions);
        geometry_bond_indices_.clear();
        for (std::size_t index=0;index<geometry_bonds_.size();++index) {
            const auto& bond=geometry_bonds_[index];
            geometry_bond_indices_[ordered_atom_pair(
                bond.atom_a,bond.atom_b)]=index;
        }
        geometry_cache_wavefunction_ = &wavefunction;
        geometry_graph_source_ = routed_graph;
    }
    const auto& bonds = geometry_bonds_;
    const auto& interactions = geometry_interactions_;
    const Vec3 molecular_centre = heavy_atom_centroid(wavefunction, points);

    gl::UseProgram(0);
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);
    glDepthMask(GL_TRUE);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glShadeModel(GL_SMOOTH);

    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    load_3d_camera(b, aspect, tan_half_fov);

    const float opacity = std::clamp(settings.molecule_opacity, 0.05f, 1.0f);
    const float bond_radius = (settings.style == MoleculeStyle::MediumBallAndStick ? 0.00540f : 0.00415f) *
                              std::clamp(settings.bond_scale, 0.35f, 3.0f);
    const Vec3 bond_colour{0.73f, 0.76f, 0.81f};
    const Vec3 delocalised_colour{0.46f, 0.69f, 0.98f};
    const Vec3 coordination_colour{0.35f, 0.77f, 0.92f};
    const Vec3 multicentre_colour{0.98f, 0.66f, 0.24f};
    const Vec3 polyhedral_cage_colour{0.42f, 0.88f, 0.62f};
    const Vec3 hydrogen_bond_colour{0.43f, 0.84f, 1.00f};
    const Vec3 noncovalent_colour{0.61f, 0.69f, 0.81f};
    const Vec3 ionic_colour{0.77f, 0.55f, 0.98f};

    auto project=[&](Vec3 p) {
        const auto delta=p-b.position;
        const float z=dot(delta,b.forward);
        return Vec3{0.5f+dot(delta,b.right)/(2*std::max(z,kNearPlane)*tan_half_fov*aspect),
                    0.5f-dot(delta,b.up)/(2*std::max(z,kNearPlane)*tan_half_fov),z};
    };
    auto target=[&](GeometryTargetKind kind,std::size_t index,Vec3 a,Vec3 c,float radius,bool segment) {
        const auto pa=project(a),pc=project(c);
        if(pa.z<=kNearPlane || pc.z<=kNearPlane)return;
        geometry_targets_.push_back({kind,index,pa.x,pa.y,pc.x,pc.y,radius,(pa.z+pc.z)/2,segment});
    };
    auto center=[&](const std::vector<std::size_t>& ids) {
        Vec3 result{};std::size_t count=0;
        for(auto id:ids)if(id<points.size()){result=result+points[id];++count;}
        return count?result/static_cast<float>(count):Vec3{0.5f,0.5f,0.5f};
    };
    std::map<std::pair<std::size_t,std::size_t>,const MoleculeOverlayBond*> enhanced;
    if(overlay) {
        for(const auto& bond:overlay->bonds) {
            if(bond.style==OverlayBondStyle::Coordination && !settings.show_coordination_contacts)continue;
            if(bond.style==OverlayBondStyle::Multicentre && !settings.show_multicentre_support)continue;
            enhanced[ordered_atom_pair(bond.atom_a,bond.atom_b)]=&bond;
        }
    }

    for (const auto& interaction : interactions.edges) {
        if (interaction.atom_a >= wavefunction.atoms.size() ||
            interaction.atom_b >= wavefunction.atoms.size()) {
            continue;
        }
        const auto pair=ordered_atom_pair(interaction.atom_a,interaction.atom_b);
        // A hyperedge is an additional analysis layer, not an all-pairs mask
        // that erases ordinary connectivity between its members.
        if(const auto replacement=enhanced.find(pair);replacement!=enhanced.end()) {
            const auto style=replacement->second->style;
            if((interaction.kind==InteractionKind::CovalentConnectivity &&
                (style==OverlayBondStyle::Covalent || style==OverlayBondStyle::Delocalised)) ||
               (interaction.kind==InteractionKind::CoordinationContact && style==OverlayBondStyle::Coordination))continue;
        }
        const auto visual_style = interaction_visual_style(interaction.kind, settings);
        if (visual_style == InteractionVisualStyle::Hidden) continue;

        const Atom& atom_a = wavefunction.atoms[interaction.atom_a];
        const Atom& atom_b = wavefunction.atoms[interaction.atom_b];
        if (!settings.show_hydrogens &&
            (atom_a.atomic_number == 1 || atom_b.atomic_number == 1)) {
            continue;
        }

        const Vec3 a = points[interaction.atom_a];
        const Vec3 c = points[interaction.atom_b];
        switch (visual_style) {
            case InteractionVisualStyle::OrdinaryBond: {
                draw_cylinder(a, c, bond_radius, bond_colour, opacity * 0.94f, b);
                const auto found = geometry_bond_indices_.find(ordered_atom_pair(
                    interaction.atom_a, interaction.atom_b));
                if (found != geometry_bond_indices_.end() &&
                    found->second<geometry_bonds_.size() &&
                    geometry_bonds_[found->second].delocalised) {
                    const Vec3 inward = inward_delocalisation_offset(
                        a, c, molecular_centre, bond_radius * 2.35f);
                    draw_dashed_cylinder(a + inward, c + inward,
                                         bond_radius * 0.58f,
                                         delocalised_colour, opacity * 0.96f, b);
                }
                break;
            }
            case InteractionVisualStyle::CoordinationDash:
                // Long, close-spaced segments read as a definite structural
                // connection while remaining distinct from a covalent stick.
                draw_dashed_cylinder(a, c, bond_radius * 0.78f,
                                     coordination_colour, opacity * 0.92f, b,
                                     6, 0.76f);
                break;
            case InteractionVisualStyle::MulticentreDash:
                // A thinner/finer support line avoids reducing a true
                // multicentre hyperedge to an ordinary pairwise bond.
                draw_dashed_cylinder(a, c, bond_radius * 0.46f,
                                     multicentre_colour, opacity * 0.82f, b,
                                     11, 0.43f);
                break;
            case InteractionVisualStyle::PolyhedralCageDash:
                // Cage support is a global skeletal/topological inference,
                // not an individual 3c electron channel.
                draw_dashed_cylinder(a,c,bond_radius*0.38f,
                                     polyhedral_cage_colour,opacity*0.72f,b,
                                     8,0.58f);
                break;
            case InteractionVisualStyle::HydrogenBondDots:
                draw_dashed_cylinder(a, c, bond_radius * 0.38f,
                                     hydrogen_bond_colour, opacity * 0.64f, b,
                                     13, 0.25f);
                break;
            case InteractionVisualStyle::NoncovalentDots:
                draw_dashed_cylinder(a, c, bond_radius * 0.32f,
                                     noncovalent_colour, opacity * 0.50f, b,
                                     12, 0.22f);
                break;
            case InteractionVisualStyle::IonicDash:
                draw_dashed_cylinder(a, c, bond_radius * 0.42f,
                                     ionic_colour, opacity * 0.60f, b,
                                     9, 0.34f);
                break;
            case InteractionVisualStyle::Hidden:
                break;
        }
    }

    if(overlay) {
        for(const auto& bond:overlay->bonds) {
            if(bond.style==OverlayBondStyle::Coordination && !settings.show_coordination_contacts)continue;
            if(bond.style==OverlayBondStyle::Multicentre && !settings.show_multicentre_support)continue;
            if(bond.atom_a>=points.size() || bond.atom_b>=points.size())continue;
            if(!settings.show_hydrogens && (wavefunction.atoms[bond.atom_a].atomic_number==1 ||
                                           wavefunction.atoms[bond.atom_b].atomic_number==1))continue;
            const auto a=points[bond.atom_a],c=points[bond.atom_b];
            const auto existing=geometry_bond_indices_.find(ordered_atom_pair(bond.atom_a,bond.atom_b));
            const bool delocalised=bond.style==OverlayBondStyle::Delocalised ||
                ((bond.style==OverlayBondStyle::Covalent || bond.style==OverlayBondStyle::Unresolved) &&
                 existing!=geometry_bond_indices_.end() &&
                 geometry_bonds_[existing->second].delocalised);
            const Vec3 colour=bond.selected?Vec3{1,0.85f,0.2f}:bond.style==OverlayBondStyle::Coordination?
                coordination_colour:delocalised?delocalised_colour:bond_colour;
            const float r=bond_radius*(bond.selected?1.25f:1.0f);
            Vec3 transverse=cross(c-a,b.forward);
            if(length(transverse)<1e-6f)transverse=b.right;
            transverse=normalise(transverse);
            if(bond.style==OverlayBondStyle::Multicentre ||
               (bond.style==OverlayBondStyle::Unresolved && !delocalised)) {
                draw_dashed_cylinder(a,c,r,colour,opacity,b,8,
                    bond.style==OverlayBondStyle::Unresolved?0.35f:0.6f);
            } else {
                // Preserve the validated delocalised-network representation.
                // A particular NBO Lewis structure remains inspectable as evidence.
                const int count=delocalised?1:std::clamp(bond.multiplicity,1,6);
                for(int line=0;line<count;++line){
                    const auto offset=transverse*((line-0.5f*(count-1))*r*3.5f);
                    draw_cylinder(a+offset,c+offset,r,colour,opacity,b);
                }
                if(delocalised)
                    draw_dashed_cylinder(a+transverse*(r*3.5f),c+transverse*(r*3.5f),r*0.6f,
                                         delocalised_colour,opacity,b,10,0.45f);
            }
            target(GeometryTargetKind::Bond,bond.evidence_index,a,c,10.0f,true);
        }
        for(const auto& group:overlay->multicentre) {
            if(!settings.show_multicentre_support)continue;
            if(group.atoms.size()<3)continue;
            const auto hub=center(group.atoms);
            const auto colour=group.selected?Vec3{1,0.85f,0.2f}:multicentre_colour;
            for(auto atom:group.atoms)if(atom<points.size() &&
                (settings.show_hydrogens || wavefunction.atoms[atom].atomic_number!=1))
                draw_dashed_cylinder(points[atom],hub,bond_radius*0.7f,colour,opacity,b,9,0.55f);
            draw_sphere(hub,bond_radius*1.8f,colour,opacity,b);
            target(GeometryTargetKind::Multicentre,group.evidence_index,hub,hub,14.0f,false);
        }
        for(const auto& relation:overlay->relations) {
            if(relation.donor_atoms.empty() || relation.acceptor_atoms.empty())continue;
            auto a=center(relation.donor_atoms),c=center(relation.acceptor_atoms);
            // A visible offset separates an interaction arrow from the skeleton.
            const auto offset=b.up*(bond_radius*7.0f);
            a=a+offset;c=c+offset;
            if(length(c-a)<bond_radius*3)c=c+b.right*(bond_radius*14);
            const auto direction=normalise(c-a);
            auto side=normalise(cross(direction,b.forward));
            if(length(side)<0.1f)side=b.right;
            const float head=std::min(length(c-a)*0.22f,bond_radius*9);
            const auto colour=relation.selected?Vec3{1,0.86f,0.18f}:Vec3{0.2f,0.92f,0.65f};
            draw_dashed_cylinder(a,c,bond_radius*0.55f,colour,opacity,b,12,0.75f);
            draw_cylinder(c,c-direction*head+side*head*0.4f,bond_radius*0.55f,colour,opacity,b);
            draw_cylinder(c,c-direction*head-side*head*0.4f,bond_radius*0.55f,colour,opacity,b);
            target(GeometryTargetKind::Relation,relation.evidence_index,a,c,9.0f,true);
        }
    }

    {
        for (std::size_t i = 0; i < wavefunction.atoms.size(); ++i) {
            const Atom& atom = wavefunction.atoms[i];
            if (!settings.show_hydrogens && atom.atomic_number == 1) continue;

            const float covalent = static_cast<float>(covalent_radius_angstrom(atom.atomic_number));
            const float sphere_radius = std::clamp(
                (0.0275f + 0.0200f * covalent) * settings.atom_scale,
                0.030f, 0.108f);
            Vec3 colour=atom_colour(atom.atomic_number);
            if(overlay && overlay->colour_mode!=AtomScalarMode::Element &&
               i<overlay->atom_values.size() && overlay->atom_values[i]) {
                const auto v=std::clamp(*overlay->atom_values[i]/std::max(1e-12,overlay->scalar_range),-1.0,1.0);
                const Vec3 neutral{0.88f,0.88f,0.80f};
                const Vec3 end=v<0?Vec3{0.52f,0.23f,0.78f}:Vec3{0.96f,0.70f,0.13f};
                colour=neutral*(1-static_cast<float>(std::abs(v)))+end*static_cast<float>(std::abs(v));
            }
            const bool selected=overlay && std::find(overlay->selected_atoms.begin(),overlay->selected_atoms.end(),i)!=overlay->selected_atoms.end();
            if(selected)draw_sphere(points[i],sphere_radius*1.17f,{1,0.83f,0.15f},opacity*0.7f,b);
            if(settings.style==MoleculeStyle::MediumBallAndStick || selected ||
               (overlay && overlay->colour_mode!=AtomScalarMode::Element))
                draw_sphere(points[i],sphere_radius,colour,opacity,b);
            const auto projected=project(points[i]);
            const float hit_radius=std::max(7.0f,sphere_radius*framebuffer_height/(2*std::max(projected.z,kNearPlane)*tan_half_fov));
            target(GeometryTargetKind::Atom,i,points[i],points[i],hit_radius,false);
        }
    }

    glPopMatrix();
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);

    glDisable(GL_BLEND);
    glDisable(GL_DEPTH_TEST);
    glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
}

std::optional<GeometryTarget> VolumeRenderer::pick_geometry(float x,float y) const {
    std::optional<GeometryTarget> result;
    float best=std::numeric_limits<float>::infinity();
    for(const auto& target:geometry_targets_) {
        const float ax=target.x*geometry_width_,ay=target.y*geometry_height_;
        const float dx=(target.end_x-target.x)*geometry_width_,dy=(target.end_y-target.y)*geometry_height_;
        const float px=x*geometry_width_,py=y*geometry_height_;
        const float t=target.segment?std::clamp(((px-ax)*dx+(py-ay)*dy)/std::max(dx*dx+dy*dy,1e-12f),0.0f,1.0f):0.0f;
        const float distance=std::hypot(px-ax-t*dx,py-ay-t*dy);
        if(distance>target.radius)continue;
        const float priority=target.kind==GeometryTargetKind::Atom?0.0f:target.kind==GeometryTargetKind::Relation?2.0f:4.0f;
        const float score=priority+distance/std::max(target.radius,1.0f)+target.depth*0.001f;
        if(score<best){best=score;result=target;}
    }
    return result;
}

} // namespace cov
