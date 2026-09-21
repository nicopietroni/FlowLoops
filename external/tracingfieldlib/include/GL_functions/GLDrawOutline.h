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

#ifndef GL_DRAW_OUTLINE
#define GL_DRAW_OUTLINE

#include "GLDrawBase.h"
#include "../point3.h"
#include <iostream>

namespace GLDraw {

template <class ScalarType>
void DrawOutline(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                 const std::vector<std::vector<int>> &Faces,
                 ScalarType thickness,ScalarType threshold=0.002) 
{
  glPushAttrib(GL_ALL_ATTRIB_BITS);
  // Disable color writes to prevent the mesh from being drawn
  glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);

  // Render the mesh to update the depth buffer

  glDisable(GL_LIGHTING);
  glDisable(GL_LIGHT0);
  glBegin(GL_TRIANGLES);
  for (size_t i = 0; i < Faces.size(); i++) {
    for (size_t j = 0; j < Faces[i].size(); j++) {
      int IndexV = Faces[i][j];
      Geo::Point3<ScalarType> Pos = VertPos[IndexV];
      glVertex(Pos);
    }
  }
  glEnd();

  // Re-enable color writes
  glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);

  // Get the matrices and viewport
  GLdouble modelview[16];
  GLdouble projection[16];
  GLint viewport[4];

  glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
  glGetDoublev(GL_PROJECTION_MATRIX, projection);
  glGetIntegerv(GL_VIEWPORT, viewport);

  int width = viewport[2];
  int height = viewport[3];

  // Allocate memory for the depth buffer
  float *depthBuffer = new float[width * height];
  glReadPixels(0, 0, width, height, GL_DEPTH_COMPONENT, GL_FLOAT, depthBuffer);

  // Edge detection parameters
  std::vector<std::pair<int, int>> edgePixels;
  // float threshold = 0.005f; // Adjust this value based on your mesh scale

  // Perform edge detection on the depth buffer
  //float MaxG=0;
  for (int y = 1; y < height - 1; y++) {
    for (int x = 1; x < width - 1; x++) {
      int index = y * width + x;
      float depth = depthBuffer[index];

      float depthRight = depthBuffer[y * width + (x + 1)];
      float depthLeft = depthBuffer[y * width + (x - 1)];
      float depthUp = depthBuffer[(y + 1) * width + x];
      float depthDown = depthBuffer[(y - 1) * width + x];

      float dx = fabs(depthRight - depthLeft);
      float dy = fabs(depthUp - depthDown);

      float gradient = sqrt(dx * dx + dy * dy);
      //MaxG=std::max(MaxG,gradient);
     
      if (gradient > threshold) {
        // Edge detected at (x, y)
        edgePixels.push_back(std::make_pair(x, y));
      }
    }
  }
  //std::cout<<" MAx G:"<<MaxG<<std::endl;
  
  // Convert edge pixels to object coordinates
  std::vector<GLdouble> edgeVertices;

  for (const auto &pixel : edgePixels) {
    int x = pixel.first;
    int y = pixel.second;

    // Read the depth value at (x, y)
    //float depth;
    //glReadPixels(x, y, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &depth);
    int index = y * width + x;
    float depth = depthBuffer[index];
      
    GLdouble winX = (GLdouble)x;
    GLdouble winY = (GLdouble)y;
    GLdouble winZ = (GLdouble)depth;

    // Adjust winY if necessary (OpenGL's coordinate origin is at the
    // bottom-left corner) winY = viewport[3] - winY; // Uncomment if necessary

    // Variables to store the object coordinates
    GLdouble objX, objY, objZ;

    // Unproject the window coordinates to object coordinates
    gluUnProject(winX, winY, winZ, modelview, projection, viewport, &objX,
                 &objY, &objZ);

    // Store the object coordinates
    edgeVertices.push_back(objX);
    edgeVertices.push_back(objY);
    edgeVertices.push_back(objZ);
  }

  delete[] depthBuffer; // Free the depth buffer memory

 // Clear the depth buffer if necessary
  glClear(GL_DEPTH_BUFFER_BIT);

  // Now draw the edge vertices in object space
  glColor3f(0.0f, 0.0f, 0.0f); // Set edge color to black
  glPointSize(10.0f);
  glBegin(GL_POINTS);
  for (size_t i = 0; i < edgeVertices.size(); i += 3) {
    glVertex3d(edgeVertices[i], edgeVertices[i + 1], edgeVertices[i + 2]);
  }

  glEnd();
  glPopAttrib();
}

}; // namespace GLDraw

#endif
