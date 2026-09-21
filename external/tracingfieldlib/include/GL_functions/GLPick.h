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

#ifndef GL_PICK
#define GL_PICK

#include <GL/glew.h>
#include <vector>
#include <algorithm>
#include <GL/glew.h>
#include <point3.h>
#include <Eigen/Dense>

//assert(Connectivity.size()==FaceNormal.size());
//glEnable(GL_LIGHTING);
//glEnable(GL_LIGHT0);
//for (size_t i=0;i<Connectivity.size();i++)
//{
//    glBegin(GL_POLYGON);
//    if (!SmoothShade)
//        glNormal(FaceNormal[i]);
//    for (size_t j=0;j<Connectivity[i].size();j++)
//    {
//        int IndexV=Connectivity[i][j];
//        assert(IndexV>=0);
//        assert(IndexV<VertPos.size());
//        if (SmoothShade)
//            glNormal(VertNormal[IndexV]);
//        Geo::Point3<ScalarType> Pos=VertPos[IndexV];
//        glVertex(Pos);
//    }
//    glEnd();
//}

namespace GLDraw
{

template <class ScalarType>
class GLPick
{
    //    static int PickFace(const int & x, const int &y,
    //                        std::vector<Geo::Point3<ScalarType > > VertPos,
    //                        std::vector<std::vector<int> > Connectivity,
    //                        std::vector<int> &PickedFaces,
    //                        int width=4,int height=4)
    //        {
    //            PickedFaces.clear();
    //            long hits;
    //            int sz = int(Connectivity.size())*5;
    //            GLuint *selectBuf =new GLuint[sz];
    //            glSelectBuffer(sz, selectBuf);
    //            glRenderMode(GL_SELECT);
    //            glInitNames();

    //            /* Because LoadName() won't work with no names on the stack */
    //            glPushName(-1);
    //            double mp[16];

    //            GLint viewport[4];
    //            glGetIntegerv(GL_VIEWPORT,viewport);
    //            //glPushAttrib(GL_TRANSFORM_BIT);
    //            glMatrixMode(GL_PROJECTION);
    //            glGetDoublev(GL_PROJECTION_MATRIX ,mp);
    //            glPushMatrix();
    //            glLoadIdentity();
    //            //gluPickMatrix(x, viewport[3]-y, 4, 4, viewport);
    //            gluPickMatrix(x, y, width, height, viewport);
    //            glMultMatrixd(mp);

    //            glMatrixMode(GL_MODELVIEW);
    //            glPushMatrix();
    //            int cnt=0;

    //            glDisable(GL_LIGHTING);
    //            for (size_t i=0;i<Connectivity.size();i++)
    //            {
    //                glLoadName(cnt);
    //                glBegin(GL_POLYGON);
    //                for (size_t j=0;j<Connectivity[i].size();j++)
    //                {
    //                    int IndexV=Connectivity[i][j];
    //                    assert(IndexV>=0);
    //                    assert(IndexV<VertPos.size());
    ////                    if (SmoothShade)
    ////                        glNormal(VertNormal[IndexV]);
    //                    Geo::Point3<ScalarType> Pos=VertPos[IndexV];
    //                    glVertex(Pos);
    //                }
    //                glEnd();
    //                cnt++;
    //            }

    ////            for(size_t i=0;i<Connectivity.size();i++)
    ////            {

    ////                    glLoadName(cnt);
    ////                    draw_func(*ei);
    ////                    cnt++;
    ////            }

    //            glPopMatrix();
    //            glMatrixMode(GL_PROJECTION);
    //            glPopMatrix();
    //            glMatrixMode(GL_MODELVIEW);
    //            hits = glRenderMode(GL_RENDER);

    //            if (hits <= 0)
    //                return 0;

    //            std::vector< std::pair<double,unsigned int> > H;
    //            for(int ii=0;ii<hits;ii++){
    //                H.push_back( std::pair<double,unsigned int>(selectBuf[ii*4+1]/4294967295.0,selectBuf[ii*4+3]));
    //            }
    //            std::sort(H.begin(),H.end());

    //            PickedFaces.resize(H.size());
    //            for(int ii=0;ii<hits;ii++){
    //                typename TO_PICK_CONT_TYPE::iterator ei=m.begin();
    //                std::advance(ei ,H[ii].second);
    //                result[ii]=&*ei;
    //            }
    //            glPopAttrib();
    //            delete [] selectBuf;
    //            return int(result.size());
    //        }

    static int ClosestVertIteratively(const Geo::Point3<ScalarType> &TestPos,
                                      const std::vector<Geo::Point3<ScalarType > > &VertPos)
    {
        ScalarType minD=std::numeric_limits<ScalarType>::max();
        int minIdx=-1;
        for (size_t i=0;i<VertPos.size();i++)
        {
            ScalarType CurrD=(TestPos-VertPos[i]).Norm();
            if (CurrD>minD)continue;
            minIdx=i;
            minD=CurrD;
        }
        return minIdx;
    }

//    static int ClosestFaceIteratively(const Geo::Point3<ScalarType> &TestPos,
//                                      const std::vector<std::vector<int > > &Connectivity)
//    {
////        ScalarType minD=std::numeric_limits<ScalarType>::max();
////        int minIdx=-1;
////        for (size_t i=0;i<VertPos.size();i++)
////        {
////            ScalarType CurrD=(TestPos-VertPos[i]).Norm();
////            if (CurrD>minD)continue;
////            minIdx=i;
////            minD=CurrD;
////        }
////        return minIdx;
//    }

public:

    static bool Pick(const int & x,const int &y,
                     Geo::Point3<ScalarType> &pp,
                     ScalarType &DepthVal)
    {
        GLdouble res[3];
        GLdouble mm[16],pm[16]; GLint vp[4];
        glGetDoublev(GL_MODELVIEW_MATRIX,mm);
        glGetDoublev(GL_PROJECTION_MATRIX,pm);
        glGetIntegerv(GL_VIEWPORT,vp);

        GLfloat   pix;
        glReadPixels(x,y,1,1,GL_DEPTH_COMPONENT,GL_FLOAT,&pix);
        GLfloat depthrange[2]={0,0};
        glGetFloatv(GL_DEPTH_RANGE,depthrange);
        if(pix==depthrange[1]) return false;

        DepthVal=pix;
        gluUnProject(x,y,pix,mm,pm,vp,&res[0],&res[1],&res[2]);
        pp=Geo::Point3<ScalarType>((ScalarType)res[0],(ScalarType)res[1],(ScalarType)res[2]);
        return true;
    }

    static Geo::Point3<ScalarType> Unproject(const int & x,const int &y,
                                             const ScalarType & DepthVal)
    {
        GLdouble res[3];
        GLdouble mm[16],pm[16]; GLint vp[4];
        glGetDoublev(GL_MODELVIEW_MATRIX,mm);
        glGetDoublev(GL_PROJECTION_MATRIX,pm);
        glGetIntegerv(GL_VIEWPORT,vp);

        gluUnProject(x,y,DepthVal,mm,pm,vp,&res[0],&res[1],&res[2]);
        Geo::Point3<ScalarType> pp=Geo::Point3<ScalarType>((ScalarType)res[0],(ScalarType)res[1],(ScalarType)res[2]);
        return pp;
    }

    static int PickVert(const int & x, const int &y,
                        const std::vector<Geo::Point3<ScalarType > > &VertPos,
                        ScalarType &DepthVal)
    {
        Geo::Point3<ScalarType> pp;
        bool has_pick=Pick(x,y,pp,DepthVal);
        if (!has_pick)return -1;
        return (ClosestVertIteratively(pp,VertPos));
    }

};

} // end namespace vcg

#endif
