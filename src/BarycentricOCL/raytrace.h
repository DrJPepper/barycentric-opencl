#pragma once
#include "util.h"

class RTSettings;
class Intersectable;
class Box;
class SMFModel;
class ACGBVH;
class SubSMF;

enum Shape { SMF, SUBSMF, BOX, HIERARCH };

class RTSettings {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        RTSettings(Vector4d &rayBox);
        void updateRayBox(Vector4d &rayBox);
        int x, y, maxObjCount;
        double theta, d, h, Sj, Sk, uStart, vStart, uLength, vLength, scale;
        Vector3d cameraLoc, Zv, Vup, Xv, Yv, P00, rotate;
        MatrixXd vertices;
        MatrixXi faces;
        string inputFileName;
};

// Interface for objects in scene
class Intersectable {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        virtual bool intersectPixel(Vector3d ray, RTSettings *setts, int pixelNum,
                                    int objNum) = 0;
        double scale;
        Shape type;
        Vector3d rotate;
        Vector3d minVec, maxVec;

    protected:
        void initSetts(RTSettings *settings);
};

// Bounding box implementation
class Box : public Intersectable {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        Box(SubSMF *subsmf);
        Box(SMFModel *smf);
        bool intersectPixel(Vector3d ray, RTSettings *setts, int pixelNum,
                            int objNum);
        Vector3d maxVec, minVec;
};

// Triangle mesh class for intersecting
class SMFModel : public Intersectable {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        SMFModel(RTSettings *settings, std::ifstream *inSMF);
        bool intersectPixel(Vector3d ray, RTSettings *setts, int pixelNum,
                            int objNum);
        bool intersectPixel(Vector3d ray, RTSettings *setts, MatrixXi *facesIn,
                            int pixelNum, int objNum);
        MatrixXd vertices;
        MatrixXi facesStored;
        // Vector3d minVec, maxVec;

    private:
        MatrixXd normals;
};

// Workaround class to only store some faces in a node of the ACGBVH tree
// without copying everything stored in the SMFModel
class SubSMF : public Intersectable {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        SubSMF(SMFModel *smfIn, MatrixXi facesIn);
        bool intersectPixel(Vector3d ray, RTSettings *setts, int pixelNum,
                            int objNum);
        SMFModel *smf;
        MatrixXi faces;
};

// Class to set up tree structure and hold other objects involved
class ACGBVH : public Intersectable {
    public:
        ACGBVH(Intersectable *inObj, int inMaxObjCount);
        Intersectable *obj;
        bool intersectPixel(Vector3d ray, RTSettings *setts, int pixelNum,
                            int objNum);

    private:
        int maxObjCount;
        bool isLeaf;
        Box *boundingVol;
        vector<ACGBVH *> children;
        void createChildren();
};

static vector<ACGBVH *> objects;
static MatrixXd lights;
static MatrixXd tBuffer;
static VectorXi tObjsBuffer;

void intersectForAllPixels(vector<ACGBVH *> objects, RTSettings *setts,
                           MatrixXd *result);
std::tuple<int, int, MatrixXd&> rayTrace(RTSettings *settings);
