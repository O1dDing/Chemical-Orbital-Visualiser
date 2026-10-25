/*
 * SPDX-License-Identifier: CPAL-1.0
 *
 * The contents of this file are subject to the Common Public Attribution License
 * Version 1.0 (the “License”); you may not use this file except in compliance with the
 * License. You may obtain a copy of the License at
 * https://opensource.org/license/CPAL-1.0. The License is based on the Mozilla Public
 * License Version 1.1 but Sections 14 and 15 have been added to cover use of software
 * over a computer network and provide for limited attribution for the Original
 * Developer. In addition, Exhibit A has been modified to be consistent with Exhibit B.
 *
 * Software distributed under the License is distributed on an “AS IS” basis, WITHOUT
 * WARRANTY OF ANY KIND, either express or implied. See the License for the specific
 * language governing rights and limitations under the License.
 *
 * The Original Code is Chemical Orbital Visualiser (COV).
 *
 * The Original Developer is O1dDing, who is also the Initial Developer.
 *
 * The Initial Developer of the Original Code is O1dDing. All portions of the code
 * written by O1dDing are Copyright (c) 2026 O1dDing. All Rights Reserved.
 *
 * Contributor O1dDing.
 */

#pragma once

#ifdef _WIN32
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <Windows.h>
#endif
#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif

#ifndef APIENTRY
#define APIENTRY
#endif

namespace cov::gl {

bool load();
GLenum volume_external_format();
GLint volume_internal_format();

extern GLuint (APIENTRY* CreateShader)(GLenum);
extern void (APIENTRY* ShaderSource)(GLuint, GLsizei, const char* const*, const GLint*);
extern void (APIENTRY* CompileShader)(GLuint);
extern void (APIENTRY* GetShaderiv)(GLuint, GLenum, GLint*);
extern void (APIENTRY* GetShaderInfoLog)(GLuint, GLsizei, GLsizei*, char*);
extern GLuint (APIENTRY* CreateProgram)();
extern void (APIENTRY* AttachShader)(GLuint, GLuint);
extern void (APIENTRY* LinkProgram)(GLuint);
extern void (APIENTRY* GetProgramiv)(GLuint, GLenum, GLint*);
extern void (APIENTRY* GetProgramInfoLog)(GLuint, GLsizei, GLsizei*, char*);
extern void (APIENTRY* UseProgram)(GLuint);
extern void (APIENTRY* DeleteShader)(GLuint);
extern void (APIENTRY* DeleteProgram)(GLuint);
extern GLint (APIENTRY* GetUniformLocation)(GLuint, const char*);
extern void (APIENTRY* Uniform1i)(GLint, GLint);
extern void (APIENTRY* Uniform1f)(GLint, GLfloat);
extern void (APIENTRY* Uniform2f)(GLint, GLfloat, GLfloat);
extern void (APIENTRY* Uniform3f)(GLint, GLfloat, GLfloat, GLfloat);
extern void (APIENTRY* ActiveTexture)(GLenum);
extern void (APIENTRY* TexImage3D)(GLenum, GLint, GLint, GLsizei, GLsizei, GLsizei,
                                   GLint, GLenum, GLenum, const void*);
extern void (APIENTRY* TexSubImage3D)(GLenum, GLint, GLint, GLint, GLint,
                                    GLsizei, GLsizei, GLsizei, GLenum, GLenum, const void*);

} // namespace cov::gl
