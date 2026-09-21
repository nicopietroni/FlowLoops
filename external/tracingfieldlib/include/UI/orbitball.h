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

#ifndef TRACKBALL_H
#define TRACKBALL_H

#include "eigen_interface.h"
#include "GL_functions/GLDrawBase.h"
#include <GL/glew.h>
#include <list>
#include <vector>
#include <map>
#include "point3.h"
#include "plane3.h"
#include "sphere3.h"
#include "math_utils.h"
#include "quaternion.h"
#include "similarity.h"
#include "GL_functions/GLView.h"

namespace  UI {


class OrbitBall
{
public:

    enum OrbitMode{OMInactive, OMScale, OMPan, OMSphere};

    /// The componibile states of the manipulator system.
    enum Button { BUTTON_NONE   = 0x0000, ///< No button or key pressed.
                  BUTTON_LEFT   = 0x0001, ///< Left mouse button pressed.
                  BUTTON_RIGHT  = 0x0004, ///< Right mouse button pressed.
                  WHEEL         = 0x0008, ///< Mouse wheel activated.
                  KEY_CTRL      = 0x0020, ///< Ctrl key pressed.
                  MODIFIER_MASK = 0x00FF, ///< (mask to get modifiers only)
                };


public:

    Geo::Similarity<float> last_track;

    Geo::Point3<float> last_point;

    Geo::Similarity<float> track;
    Geo::Point3<float> center;
    float radius;

    unsigned int last_time;

    friend class TrackMode;

    int current_button;

    std::map<int, OrbitMode> OrbitModes;
    OrbitMode CurrOmode;

    GLDraw::View<float> camera;

    void Translate(Geo::Point3<float> tr)
    {
        Geo::Quaternion<float> irot = track.rot;
        irot.Invert();
        Geo::Point3<float> rotP=irot.Rotate(tr)/track.sca;
        track.tra = last_track.tra + rotP;
    }

    void GetView()
    {
        camera.GetView();
    }

    void Apply ()
    {
        GLDraw::glTranslate(center);
        Eigen::Matrix4f eigMatr=track.Matrix();
        GLDraw::glMultMatrix(eigMatr);
        GLDraw::glTranslate (-center);
    }

#define LineWidthInactive 0.9f
#define LineWidthActive 1.8f

    void DrawPostApply()
    {

        if(CurrOmode==OMInactive)
            GLDraw::GLDrawOrbitSphere<float>(InverseMatrix(),center,radius,track.sca,LineWidthInactive);
        //GLDrawOrbitSphere(false);

        if(CurrOmode==OMSphere)
            GLDraw::GLDrawOrbitSphere<float>(InverseMatrix(),center,radius,track.sca,LineWidthActive);

        if(CurrOmode==OMPan)
            GLDraw::GLDrawSphereAxis<float>(InverseMatrix(),center,radius,track.sca,LineWidthActive);
        //GLDrawOrbitSphere(true);

        if(CurrOmode==OMScale)
            GLDraw::GLDrawOrbitSphere<float>(InverseMatrix(),center,radius,track.sca,LineWidthActive);

    }

    void ApplyInverse()
    {
        GLDraw::glTranslate(Geo::Point3<float>(center));
        GLDraw::glMultMatrix (track.InverseMatrix());
        GLDraw::glTranslate(-center);
    }

    Eigen::Matrix4f Matrix() {return track.Matrix();}
    Eigen::Matrix4f InverseMatrix() {return track.InverseMatrix();}

    void Reset()
    {
        track.SetIdentity();
    }

    //interface
    void MouseDown(int button) {
        //undo_track = track;
        current_button |= button;
        SetCurrentAction();
        //Hits.clear();
    }

    void MouseDown(int x, int y, int button) {
        //undo_track = track;
        current_button |= button;
        SetCurrentAction();
        last_point = Geo::Point3<float>((float)x, (float)y, 0);
        //Hits.clear();
    }

    Geo::Plane3<float> GetViewPlane (const GLDraw::View<float> &camera,
                                     const Geo::Point3<float> & center)
    {
        Geo::Point3<float> vp =camera.ViewPoint ();
        Geo::Point3<float> plnorm = vp - center;
        plnorm.Normalize();

        Geo::Plane3<float> pl(Dot(plnorm,center),Geo::Point3<float>(plnorm.X,plnorm.Y,plnorm.Z));
        return pl;
    }


    Geo::Point3<float> HitViewPlane (const Geo::Point3<float> & p)
    {
        Geo::Plane3<float> vp = GetViewPlane (camera,center);

        Geo::Point3<float> OrPlGeo=vp.Direction()*vp.Offset();
        Geo::Point3<float> OrPl(OrPlGeo.X,OrPlGeo.Y,OrPlGeo.Z);
        Geo::Line3<float> Ln = camera.ViewLineFromWindow (Geo::Point3<float> (p.X, p.Y, 0));


        Geo::Point3<float> DirPlgeo(vp.Direction().X,vp.Direction().Y,vp.Direction().Z);
        Geo::Plane3<float> PlVp(vp.Offset(),DirPlgeo);

        Geo::Point3<float> PonVP;
        PlVp.Intersection(Ln,PonVP);
        Geo::Point3<float> retP(PonVP.X,PonVP.Y,PonVP.Z);

        return retP;
    }

    bool HitHyper (Geo::Point3<float> center,
                   float radius,
                   Geo::Point3<float> viewpoint,
                   Geo::Plane3<float> viewplane,
                   Geo::Point3<float> hitOnViewplane,
                   Geo::Point3<float> & hit)
    {
        float hitplaney = (center- hitOnViewplane).Norm();
        float viewpointx =(center- viewpoint).Norm();

        float a = hitplaney / viewpointx;
        float b = -hitplaney;
        float c = radius * radius / 2.0f;
        float delta = b * b - 4 * a * c;
        //float x1, x2, xval, yval;
        float x1, xval, yval;

        if (delta > 0) {
            x1 = (-b - sqrt (delta)) / (2.0f * a);
            //x2 = (-b + sqrt (delta)) / (2.0f * a);

            xval = x1;                  // always take the minimum value solution
            yval = c / xval;            //  alternatively it also could be the other part of the equation yval=-(hitplaney/viewpointx)*xval+hitplaney;
        }
        else {
            return false;
        }
        // Computing the result in 3d space;
        Geo::Point3<float> dirRadial = hitOnViewplane - center;
        dirRadial.Normalize ();
        Geo::Point3<float> dirViewGeo=viewplane.Direction ();
        Geo::Point3<float> dirView(dirViewGeo.X,dirViewGeo.Y,dirViewGeo.Z);
        dirView.Normalize ();
        hit = center + dirRadial * yval + dirView * xval;
        return true;
    }

    Geo::Point3<float> HitSphere (const Geo::Point3<float> & p)
    {
        Geo::Line3<float> Ln= camera.ViewLineFromWindow ( Geo::Point3<float>(p.X, p.Y, 0));
        Geo::Plane3<float> vp = GetViewPlane (camera, center);
        Geo::Point3<float> hitPlane(0,0,0), //intersection view plane with point touched
                hitSphere(0,0,0),
                hitSphere1(0,0,0),
                hitSphere2(0,0,0),
                hitHyper(0,0,0);

        Geo::Sphere3<float> sphere(center,radius);

        Geo::Point3<float> hitSphere1_geo,hitSphere2_geo;
        bool resSp = sphere.Intersection(Ln, hitSphere1_geo, hitSphere2_geo);

        hitSphere1=Geo::Point3<float>(hitSphere1_geo.X,hitSphere1_geo.Y,hitSphere1_geo.Z);
        hitSphere2=Geo::Point3<float>(hitSphere2_geo.X,hitSphere2_geo.Y,hitSphere2_geo.Z);

        //sphere.
        Geo::Point3<float> viewpoint = camera.ViewPoint ();
        if (resSp == true) {
            if ((viewpoint - hitSphere1).Norm() < (viewpoint- hitSphere2).Norm())
                hitSphere = hitSphere1;
            else
                hitSphere = hitSphere2;
        }

        Geo::Point3<float> PlaneDir(vp.Direction().X,
                                    vp.Direction().Y,
                                    vp.Direction().Z);

        Geo::Plane3<float> Pl(vp.Offset(),PlaneDir);
        bool resHp;

        Geo::Point3<float> HitGeo;
        Pl.Intersection(Ln,HitGeo);

        hitPlane=Geo::Point3<float>(HitGeo.X,HitGeo.Y,HitGeo.Z);

        //  if(tb->camera.isOrtho)
        //  {
        //    exit(0);
        //    resHp= HitHyperOrtho (center, tb->radius, viewpoint, vp, hitPlane, hitHyper);
        //  }
        //  else
        resHp= HitHyper (center,radius, viewpoint, vp, hitPlane, hitHyper);

        if ((!resSp && !resHp)) {
            Geo::Point3<float> hitGeo=Ln.ClosestPoint (center);
            Geo::Point3<float> hit(hitGeo.X,hitGeo.Y,hitGeo.Z);

            //printf("closest point to line %f\n",Distance(hit,tb->center));
            return hit;
        }
        if ((resSp && !resHp))
            return hitSphere;

        if ((!resSp && resHp))
            return hitHyper;

        float angleDeg = (Geo::AngleDeg((viewpoint - center),(hitSphere - center)));


        if (angleDeg < 45)
            return hitSphere;
        else
            return hitHyper;

    }

    float getDeltaY(UI::OrbitBall * tb, Geo::Point3<float> new_point)
    {
        float ScreenHeight = float (tb->camera.viewport[3] - tb->camera.viewport[1]);
        return (new_point.Y - tb->last_point.Y) / ScreenHeight;
    }


    void MouseMove(int x, int y) {

        if(last_point.Z == -1) {//changed mode in the middle of moving
            last_point = Geo::Point3<float>((float)x, (float)y, 0);
            return;
        }
        if(CurrOmode == OMInactive) return;
        if(CurrOmode==OMScale) return;

        Geo::Point3<float> new_point(float(x), float(y), 0);
        if(CurrOmode==OMSphere)
        {
            Geo::Point3<float> hitOld = HitSphere (last_point);

            Geo::Point3<float> hitNew = HitSphere (new_point);

            Geo::Point3<float> axis = Cross((hitNew - center),(hitOld - center));

            axis.Normalize();

            float phi = std::max(Geo::AngleRad((hitNew - center),(hitOld - center)),
                                 ((hitNew-hitOld).Norm()/radius)) ;
            track.rot = Geo::Quaternion<float> (-phi, axis) * last_track.rot;
        }
        if(CurrOmode==OMPan)
        {
            Geo::Point3<float> hitOld = HitViewPlane (last_point);
            Geo::Point3<float> hitNew = HitViewPlane (new_point);
            Translate (hitNew - hitOld);
        }
    }

    void MouseUp(int /* x */, int /* y */, int button)
    {
        //undo_track = track;
        ButtonUp(UI::OrbitBall::Button(button));

        //current_button &= (~button);
        //SetCurrentAction();
    }


    // it assumes that a notch of 1.0 is a single step of the wheel
    void MouseWheel(float notch)
    {
        //undo_track = track;
        int buttons = current_button;
        current_button = WHEEL;// | (buttons&(KEY_SHIFT|KEY_CTRL|KEY_ALT));
        SetCurrentAction();

        if (CurrOmode == OMScale)
            track.sca *= std::pow (1.2f, -notch);
        //current_mode->Apply(this, notch);

        current_button = buttons;
        SetCurrentAction();
    }

    void ButtonDown(UI::OrbitBall::Button button)
    {
        current_button |= button;
        SetCurrentAction();
    }

    void ButtonUp(UI::OrbitBall::Button button)
    {
        current_button &= (~button);
        SetCurrentAction();
    }


    void SetCurrentAction ()
    {
        if (OrbitModes.count (current_button & MODIFIER_MASK)==0) {
            CurrOmode = OMInactive;
        } else {
            CurrOmode = OrbitModes[current_button & MODIFIER_MASK];
        }
        last_point = Geo::Point3<float> (0, 0, -1);
        last_track = track;
    }

    void setDefaultMapping()
    {
        OrbitModes[BUTTON_NONE] = OMInactive;

        OrbitModes[BUTTON_LEFT] = OMSphere;

        OrbitModes[BUTTON_LEFT | KEY_CTRL] = OMPan;

        OrbitModes[WHEEL] = OMScale;
    }


    OrbitBall()
    {
        //current_mode= new Geo::InactiveMode();
        CurrOmode=OMInactive;
        current_button=0;
        last_time=0;
        setDefaultMapping ();
        track.SetIdentity();
        radius=1.0f;
        center=Geo::Point3<float>(0,0,0);
    }

    ~OrbitBall(){}
};


}//namespace

#endif
