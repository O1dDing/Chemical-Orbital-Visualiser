#include "cov/gl_api.hpp"

#define GLFW_INCLUDE_NONE
#include <GLFW/glfw3.h>
#include <cstdlib>
#include <cstring>

namespace cov::gl {

GLenum volume_external_format() {
    const auto* version = reinterpret_cast<const char*>(glGetString(GL_VERSION));
    const auto* extensions = reinterpret_cast<const char*>(glGetString(GL_EXTENSIONS));
    // Legacy macOS OpenGL provides floating luminance textures, but may lack
    // the later texture_rg extension. Both formats sample through .r.
    const bool red = (version && std::atoi(version) >= 3) ||
        (extensions && std::strstr(extensions, "GL_ARB_texture_rg"));
    return red ? 0x1903 /* GL_RED */ : GL_LUMINANCE;
}
GLint volume_internal_format() {
    return volume_external_format() == GL_LUMINANCE ? 0x8818 /* GL_LUMINANCE32F_ARB */ : 0x822E /* GL_R32F */;
}

GLuint (APIENTRY* CreateShader)(GLenum) = nullptr;
void (APIENTRY* ShaderSource)(GLuint, GLsizei, const char* const*, const GLint*) = nullptr;
void (APIENTRY* CompileShader)(GLuint) = nullptr;
void (APIENTRY* GetShaderiv)(GLuint, GLenum, GLint*) = nullptr;
void (APIENTRY* GetShaderInfoLog)(GLuint, GLsizei, GLsizei*, char*) = nullptr;
GLuint (APIENTRY* CreateProgram)() = nullptr;
void (APIENTRY* AttachShader)(GLuint, GLuint) = nullptr;
void (APIENTRY* LinkProgram)(GLuint) = nullptr;
void (APIENTRY* GetProgramiv)(GLuint, GLenum, GLint*) = nullptr;
void (APIENTRY* GetProgramInfoLog)(GLuint, GLsizei, GLsizei*, char*) = nullptr;
void (APIENTRY* UseProgram)(GLuint) = nullptr;
void (APIENTRY* DeleteShader)(GLuint) = nullptr;
void (APIENTRY* DeleteProgram)(GLuint) = nullptr;
GLint (APIENTRY* GetUniformLocation)(GLuint, const char*) = nullptr;
void (APIENTRY* Uniform1i)(GLint, GLint) = nullptr;
void (APIENTRY* Uniform1f)(GLint, GLfloat) = nullptr;
void (APIENTRY* Uniform2f)(GLint, GLfloat, GLfloat) = nullptr;
void (APIENTRY* Uniform3f)(GLint, GLfloat, GLfloat, GLfloat) = nullptr;
void (APIENTRY* ActiveTexture)(GLenum) = nullptr;
void (APIENTRY* TexImage3D)(GLenum, GLint, GLint, GLsizei, GLsizei, GLsizei,
                            GLint, GLenum, GLenum, const void*) = nullptr;
void (APIENTRY* TexSubImage3D)(GLenum, GLint, GLint, GLint, GLint,
                             GLsizei, GLsizei, GLsizei, GLenum, GLenum, const void*) = nullptr;

template <typename T>
bool load_one(T& target, const char* name) {
    target = reinterpret_cast<T>(glfwGetProcAddress(name));
    return target != nullptr;
}

bool load() {
    bool ok = true;
    ok &= load_one(CreateShader, "glCreateShader");
    ok &= load_one(ShaderSource, "glShaderSource");
    ok &= load_one(CompileShader, "glCompileShader");
    ok &= load_one(GetShaderiv, "glGetShaderiv");
    ok &= load_one(GetShaderInfoLog, "glGetShaderInfoLog");
    ok &= load_one(CreateProgram, "glCreateProgram");
    ok &= load_one(AttachShader, "glAttachShader");
    ok &= load_one(LinkProgram, "glLinkProgram");
    ok &= load_one(GetProgramiv, "glGetProgramiv");
    ok &= load_one(GetProgramInfoLog, "glGetProgramInfoLog");
    ok &= load_one(UseProgram, "glUseProgram");
    ok &= load_one(DeleteShader, "glDeleteShader");
    ok &= load_one(DeleteProgram, "glDeleteProgram");
    ok &= load_one(GetUniformLocation, "glGetUniformLocation");
    ok &= load_one(Uniform1i, "glUniform1i");
    ok &= load_one(Uniform1f, "glUniform1f");
    ok &= load_one(Uniform2f, "glUniform2f");
    ok &= load_one(Uniform3f, "glUniform3f");
    ok &= load_one(ActiveTexture, "glActiveTexture");
    ok &= load_one(TexImage3D, "glTexImage3D");
    ok &= load_one(TexSubImage3D, "glTexSubImage3D");
    return ok;
}

} // namespace cov::gl
