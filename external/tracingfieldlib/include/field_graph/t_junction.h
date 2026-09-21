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

#ifndef T_JUNCTION
#define T_JUNCTION

namespace Geo {

//     VT
//      |
//      |
// V0---V---V1
struct TJunction 
{
  int V,V0,V1,VT;

  TJunction(){
    V=-1;
    V0=-1;
    V1=-1;
    VT=-1;
  }
  
  TJunction(int _V,int _V0,int _V1,int _VT)
  {
    V=_V;
    V0=_V0;
    V1=_V1;
    VT=_VT;
  }
};

} // namespace Geo
#endif
