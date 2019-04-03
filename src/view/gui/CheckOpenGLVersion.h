/**
\file CheckOpenGLVersion.h

Checks for relevant OpenGL version on the running platform

Reference: https://github.com/Kitware/VTK/blob/master/Rendering/OpenGL2/vtkTestOpenGLVersion.cxx
*/

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#if WIN32
#include <Windows.h>
#include <GL/gl.h>

typedef HGLRC(WINAPI * PFNWGLCREATECONTEXTATTRIBSARBPROC) (HDC hDC, HGLRC hShareContext, const int* attribList);
#define WGL_CONTEXT_MAJOR_VERSION_ARB 0x2091
#define WGL_CONTEXT_MINOR_VERSION_ARB 0x2092
#define WGL_CONTEXT_FLAGS_ARB 0x2094
#define GL_MAJOR_VERSION 0x821B
#define GL_MINOR_VERSION 0x821C
#endif

class CheckOpenGLVersion {

public:

  CheckOpenGLVersion();
#if WIN32
  CheckOpenGLVersion(HINSTANCE hInstance);
#endif
  bool hasVersion_3_2();

  std::string version,
    renderer,
    vendor;

private:

#if WIN32
  HINSTANCE hInstance;

  PIXELFORMATDESCRIPTOR pixelFormatDescriptor;
#endif

  int glMajorVersion = 1;

  int glMinorVersion = 0;

};