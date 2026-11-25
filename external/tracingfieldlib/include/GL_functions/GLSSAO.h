/*
 * //Non-Commercial License Agreement
 * //Copyright (c) 2025, Nico Pietroni
 * //
 * //1. Permission is hereby granted to use, copy, modify, and distribute
 * //   this software for academic research, teaching, and personal
 * //   non-commercial purposes, provided that this copyright notice is
 * //   included in all copies or substantial portions of the software.
 * //
 * //2. Commercial use of this software, including but not limited to use in
 * //   commercial products, services, tools, or consulting, is strictly
 * //   prohibited without prior written permission from the copyright holder.
 * //
 * //3. This software is provided "as is", without warranty of any kind,
 * //   express or implied, including but not limited to the warranties of
 * //   merchantability, fitness for a particular purpose, and
 * //   noninfringement. In no event shall the authors or copyright holders
 * //   be liable for any claim, damages, or other liability, arising from,
 * //   out of, or in connection with the software or the use or other
 * //   dealings in the software.
 * //
 * //4. For commercial licensing, please contact:
 * //   Nico Pietroni
 * //   Email: nico.pietroni@gmail.com
 */

#ifndef GL_SSAO
#define GL_SSAO

#include <GL/glew.h>
#include <iostream>
#include <string>

namespace GLDraw {

// SCREEN SPACE AMBIENT OCCLUSION (legacy OpenGL 2 + GLSL 1.20, framebuffer
// objects of EXT_framebuffer_object): the scene is drawn into an offscreen
// buffer (color + depth textures), then copied on the screen darkened where
// the depth around a pixel is closer to the camera (creases, folds, surfaces
// close to each other). Usage:
//   SSAO.Begin(width, height);  // before drawing the scene (viewport, clear)
//   ... draw the scene ...
//   SSAO.End();                 // the scene on the screen, then the GUI
// zNear / zFar: the ones of the projection. When disabled (or not supported)
// Begin / End do nothing.
class GLScreenSpaceAO {
public:
  bool Enabled = true;
  // projection planes (gluPerspective)
  float zNear = 0.1f, zFar = 100.f;
  // radius of the samples in pixels
  float RadiusPix = 25.f;
  // depth range of the occlusion (scene units): farther occluders ignored
  float Range = 0.3f;
  // darkening (0 none, 1 black where fully occluded)
  float Strength = 0.8f;

  ~GLScreenSpaceAO() { Release(); }

  void Begin(int W, int H) {
    Active = false;
    if (!Enabled || (W <= 0) || (H <= 0))
      return;
    if (!Init(W, H))
      return;
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, FBO);
    glViewport(0, 0, W, H);
    Active = true;
  }

  void End() {
    if (!Active)
      return;
    Active = false;
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
    glPushAttrib(GL_ALL_ATTRIB_BITS);
    glViewport(0, 0, Width, Height);
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_LIGHTING);
    glDisable(GL_BLEND);
    glDisable(GL_CULL_FACE);
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();

    glUseProgram(Program);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, ColorTex);
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, DepthTex);
    glUniform1i(glGetUniformLocation(Program, "colorTex"), 0);
    glUniform1i(glGetUniformLocation(Program, "depthTex"), 1);
    glUniform2f(glGetUniformLocation(Program, "texel"), 1.f / Width,
                1.f / Height);
    glUniform1f(glGetUniformLocation(Program, "zNear"), zNear);
    glUniform1f(glGetUniformLocation(Program, "zFar"), zFar);
    glUniform1f(glGetUniformLocation(Program, "radius"), RadiusPix);
    glUniform1f(glGetUniformLocation(Program, "range"), Range);
    glUniform1f(glGetUniformLocation(Program, "strength"), Strength);

    glBegin(GL_QUADS);
    glTexCoord2f(0, 0);
    glVertex2f(-1, -1);
    glTexCoord2f(1, 0);
    glVertex2f(1, -1);
    glTexCoord2f(1, 1);
    glVertex2f(1, 1);
    glTexCoord2f(0, 1);
    glVertex2f(-1, 1);
    glEnd();

    glUseProgram(0);
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, 0);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, 0);
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
    glPopAttrib();
  }

private:
  GLuint FBO = 0, ColorTex = 0, DepthTex = 0, Program = 0;
  int Width = 0, Height = 0;
  bool Active = false;
  bool Failed = false;

  void Release() {
    if (FBO)
      glDeleteFramebuffersEXT(1, &FBO);
    if (ColorTex)
      glDeleteTextures(1, &ColorTex);
    if (DepthTex)
      glDeleteTextures(1, &DepthTex);
    FBO = ColorTex = DepthTex = 0;
  }

  // buffers (again when the size changes) and shader (once)
  bool Init(int W, int H) {
    if (Failed)
      return false;
    if (!GLEW_EXT_framebuffer_object || !GLEW_VERSION_2_0) {
      std::cout << "SSAO: framebuffer objects or GLSL not supported"
                << std::endl;
      Failed = true;
      return false;
    }
    if ((Program == 0) && !BuildProgram()) {
      Failed = true;
      return false;
    }
    if (FBO && (W == Width) && (H == Height))
      return true;
    Release();
    Width = W;
    Height = H;
    glGenTextures(1, &ColorTex);
    glBindTexture(GL_TEXTURE_2D, ColorTex);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, W, H, 0, GL_RGBA,
                 GL_UNSIGNED_BYTE, nullptr);
    glGenTextures(1, &DepthTex);
    glBindTexture(GL_TEXTURE_2D, DepthTex);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_DEPTH_TEXTURE_MODE, GL_LUMINANCE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_MODE, GL_NONE);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT24, W, H, 0,
                 GL_DEPTH_COMPONENT, GL_UNSIGNED_INT, nullptr);
    glBindTexture(GL_TEXTURE_2D, 0);
    glGenFramebuffersEXT(1, &FBO);
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, FBO);
    glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT,
                              GL_TEXTURE_2D, ColorTex, 0);
    glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT,
                              GL_TEXTURE_2D, DepthTex, 0);
    GLenum Status = glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
    if (Status != GL_FRAMEBUFFER_COMPLETE_EXT) {
      std::cout << "SSAO: framebuffer not complete (" << Status << ")"
                << std::endl;
      Release();
      Failed = true;
      return false;
    }
    return true;
  }

  static GLuint Compile(GLenum Type, const char *Src) {
    GLuint S = glCreateShader(Type);
    glShaderSource(S, 1, &Src, nullptr);
    glCompileShader(S);
    GLint Ok = 0;
    glGetShaderiv(S, GL_COMPILE_STATUS, &Ok);
    if (!Ok) {
      char Log[2048];
      glGetShaderInfoLog(S, sizeof(Log), nullptr, Log);
      std::cout << "SSAO shader: " << Log << std::endl;
      glDeleteShader(S);
      return 0;
    }
    return S;
  }

  bool BuildProgram() {
    const char *VS = "#version 120\n"
                     "void main() {\n"
                     "  gl_TexCoord[0] = gl_MultiTexCoord0;\n"
                     "  gl_Position = gl_Vertex;\n"
                     "}\n";
    // occlusion from the depth only: samples on a spiral around the pixel,
    // a sample closer to the camera (within range) occludes
    const char *FS =
        "#version 120\n"
        "uniform sampler2D colorTex;\n"
        "uniform sampler2D depthTex;\n"
        "uniform vec2 texel;\n"
        "uniform float zNear, zFar, radius, range, strength;\n"
        "float lin(float d) {\n"
        "  float z = d * 2.0 - 1.0;\n"
        "  return 2.0 * zNear * zFar / (zFar + zNear - z * (zFar - zNear));\n"
        "}\n"
        "void main() {\n"
        "  vec2 uv = gl_TexCoord[0].st;\n"
        "  vec4 col = texture2D(colorTex, uv);\n"
        "  float d = texture2D(depthTex, uv).r;\n"

        "  if (d >= 1.0) { gl_FragColor = col; return; }\n"
        "  float z = lin(d);\n"
        "  float occ = 0.0;\n"
        "  const int N = 24;\n"
        "  for (int i = 0; i < N; i++) {\n"
        "    float a = float(i) * 2.399963;\n"
        "    float r = radius * sqrt((float(i) + 0.5) / float(N));\n"
        "    vec2 off = vec2(cos(a), sin(a)) * r * texel;\n"
        "    float zs = lin(texture2D(depthTex, uv + off).r);\n"
        "    float diff = z - zs;\n"
        "    occ += step(0.002 * z, diff) * (1.0 - smoothstep(0.0, range, "
        "diff));\n"
        "  }\n"
        "  occ /= float(N);\n"
        "  gl_FragColor = vec4(col.rgb * (1.0 - strength * occ), col.a);\n"
        "}\n";
    GLuint V = Compile(GL_VERTEX_SHADER, VS);
    GLuint F = Compile(GL_FRAGMENT_SHADER, FS);
    if (!V || !F)
      return false;
    Program = glCreateProgram();
    glAttachShader(Program, V);
    glAttachShader(Program, F);
    glLinkProgram(Program);
    GLint Ok = 0;
    glGetProgramiv(Program, GL_LINK_STATUS, &Ok);
    glDeleteShader(V);
    glDeleteShader(F);
    if (!Ok) {
      char Log[2048];
      glGetProgramInfoLog(Program, sizeof(Log), nullptr, Log);
      std::cout << "SSAO program: " << Log << std::endl;
      glDeleteProgram(Program);
      Program = 0;
      return false;
    }
    return true;
  }
};

} // namespace GLDraw

#endif
