# Raytracer

## Getting Started

1. Install QT Creator  
2. Install CMake  
3. Configure QT Creator:  
  - Qt `Creator->Preferences`, select `Build & Run`, go to CMake tab, and point to the installed CMake binary  
  - go to Compilers, select `gcc 4.7.2, clang 3.1, mingw gcc 4.7.2`,  
  - go to Kits, select `MinGW64` as the compiler and `GDB64` as the debugger, and `none` as the Qt version
4. Open this project, and select the CMakeLists.txt  
  - In the wizard write:  
    - Clang: `-Wdev -G "Unix Makefiles" -DCMAKE_C_COMPILER=/usr/bin/clang -DCMAKE_CXX_COMPILER=/usr/bin/clang++ -DUSE_CLANG=on`  
    - GCC Macports: `-Wdev -G "Unix Makefiles" -DCMAKE_C_COMPILER=/opt/local/bin/gcc-mp-4.8 -DCMAKE_CXX_COMPILER=/opt/local/bin/g++-mp-4.8`  
  - Run Cmake  
5. Choose `Build->Build all`  

## Running the Program

To load an render a scene fro 01_textured.json, write:
`$ ../bin/mk/04_pathtrace 01_textured.json`

To generate and view a test scene, write:
`$ ../bin/mk/04_pathtrace testscene0`

## Features
- Blurry Reflections
- Texture Tiling
- Bilinear filtering
- Mip map filtering
- Area lights
- Sphere lights
- Environment Illumination
- Indirect Lighting
- One final cool scene

## Authors

Created by Joon Cho and Josh Utter back for CS 77: Computer Graphics (Spring 2015)

### Contributions

We implemented these features on top of a Raytracer framework provided by Professor Emily Whiting of Dartmouth
