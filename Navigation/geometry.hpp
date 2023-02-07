#pragma once
#include <cmath>
#include <iostream>
using namespace std;
//#include <algorithm> ???

struct BaryCoords {
    float a, b, c;
    BaryCoords() = default;

    BaryCoords(float i, float j, float k) {
        a = i; b = j; c = k;
    }

    BaryCoords(const BaryCoords& B) {
        a = B.a; b = B.b; c = B.c;
    }

    BaryCoords& operator *=(float s) {
        a *= s; b *= s; c *= s;
        return *this;
    }

    BaryCoords& operator /=(float s) {
        s = 1.0F / s;
        a *= s; b *= s; c *= s;
        return *this;
    }

    float& operator [](int i) {
        return (&a)[i];
    }

    const float& operator [](int i) const {
        return (&a)[i];
    }
};

inline BaryCoords operator *(const BaryCoords& B, float s) {
    return BaryCoords(B.a * s, B.b * s, B.c * s);
}

inline BaryCoords operator /(const BaryCoords& B, float s) {
    s = 1.0F / s;
    return BaryCoords(B.a * s, B.b * s, B.c * s);
}

ostream& operator << (ostream& os, const BaryCoords& B) {
    return os << "{ " << B.a << ", " << B.b << ", " << B.c << " }";
}

namespace gmtry3D {
    // int array pointers with lf at the start have length stored in [0]
    // int array pointers with nl at the start do not have a stored length
    // int array pointers without lf or nl are assumed to be nl

    const int uBvctrs[3][2] = { {2, 1}, {0, 2}, {1, 0} };
    const int otherUVctr[3][3] = { {3, 2, 1}, {2, 3, 0}, {1, 0, 3} };

    struct Vector3 {
        float x, y, z;
        Vector3() = default;

        Vector3(float i, float j, float k) {
            x = i; y = j; z = k;
        }

        Vector3 clone() const {
            return Vector3(x, y, z);
        }

        Vector3& operator *=(float s) {
            x *= s; y *= s; z *= s;
            return *this;
        }

        Vector3& operator /=(float s) {
            s = 1.0F / s;
            x *= s; y *= s; z *= s;
            return *this;
        }

        Vector3& operator +=(const Vector3& b) {
            x += b.x; y += b.y; z += b.z;
            return *this;
        }

        Vector3& operator -=(const Vector3& b) {
            x -= b.x; y -= b.y; z -= b.z;
            return *this;
        }

        float& operator [](int i) {
            return (&x)[i]; // or *((&x)+i)
        }

        const float& operator [](int i) const {
            return (&x)[i];
        }
    };

    inline Vector3 operator +(const Vector3& a, const Vector3& b) {
        return Vector3(a.x + b.x, a.y + b.y, a.z + b.z);
    }

    inline Vector3 operator -(const Vector3& a, const Vector3& b) {
        return Vector3(a.x - b.x, a.y - b.y, a.z - b.z);
    }

    inline Vector3 operator -(const Vector3& a) {
        return Vector3(-(a.x), -(a.y), -(a.z));
    }

    inline Vector3 operator *(const Vector3& a, float s) {
        return Vector3(a.x * s, a.y * s, a.z * s);
    }

    inline Vector3 operator /(const Vector3& a, float s) {
        s = 1.0F / s;
        return Vector3(a.x * s, a.y * s, a.z * s);
    }

    inline float Magnitude(const Vector3& a) {
        return sqrt((a.x * a.x) + (a.y * a.y) + (a.z * a.z));
    }

    inline Vector3 Normalize(const Vector3& a) {
        return a / Magnitude(a);
    }

    ostream& operator << (ostream& os, const Vector3& a) {
        return os << "{ " << a.x << ", " << a.y << ", " << a.z << " }";
    }

    struct Matrix3 {
        float n[3][3];
        Matrix3() {
            n[0][0] = 1; n[0][1] = 0; n[0][2] = 0;
            n[1][0] = 0; n[1][1] = 1; n[1][2] = 0;
            n[2][0] = 0; n[2][1] = 0; n[2][2] = 1;
        }

        Matrix3(float n00, float n01, float n02,
            float n10, float n11, float n12,
            float n20, float n21, float n22) {
            n[0][0] = n00; n[0][1] = n01; n[0][2] = n02;
            n[1][0] = n10; n[1][1] = n11; n[1][2] = n12;
            n[2][0] = n20; n[2][1] = n21; n[2][2] = n22;
        }

        Matrix3(const Vector3& a, const Vector3& b, const Vector3& c) {
            n[0][0] = a.x; n[0][1] = a.y; n[0][2] = a.z;
            n[1][0] = b.x; n[1][1] = b.y; n[1][2] = b.z;
            n[2][0] = c.x; n[2][1] = c.y; n[2][2] = c.z;
        }

        const Matrix3 clone() const {
            return Matrix3(n[0][0], n[0][1], n[0][2],
                n[1][0], n[1][1], n[1][2],
                n[2][0], n[2][1], n[2][2]);
        }

        float& operator ()(int i, int j) {
            return n[j][i];
        }

        const float& operator ()(int i, int j) const {
            return n[i][j];
        }

        Vector3& operator ()(int i) {
            return *(reinterpret_cast <Vector3*> (n[i]));
        }

        const Vector3& operator ()(int i) const {
            return *(reinterpret_cast <const Vector3*> (n[i]));
        }

        Matrix3 T() {
            return Matrix3(n[0][0], n[1][0], n[2][0],
                n[0][1], n[1][1], n[2][1],
                n[0][2], n[1][2], n[2][2]);
        }
    };

    Matrix3 operator *(const Matrix3& A, const Matrix3& B) {
        return Matrix3(
            A.n[0][0] * B.n[0][0] + A.n[1][0] * B.n[0][1] + A.n[2][0] * B.n[0][2],
            A.n[0][1] * B.n[0][0] + A.n[1][1] * B.n[0][1] + A.n[2][1] * B.n[0][2],
            A.n[0][2] * B.n[0][0] + A.n[1][2] * B.n[0][1] + A.n[2][2] * B.n[0][2],
            A.n[0][0] * B.n[1][0] + A.n[1][0] * B.n[1][1] + A.n[2][0] * B.n[1][2],
            A.n[0][1] * B.n[1][0] + A.n[1][1] * B.n[1][1] + A.n[2][1] * B.n[1][2],
            A.n[0][2] * B.n[1][0] + A.n[1][2] * B.n[1][1] + A.n[2][2] * B.n[1][2],
            A.n[0][0] * B.n[2][0] + A.n[1][0] * B.n[2][1] + A.n[2][0] * B.n[2][2],
            A.n[0][1] * B.n[2][0] + A.n[1][1] * B.n[2][1] + A.n[2][1] * B.n[2][2],
            A.n[0][2] * B.n[2][0] + A.n[1][2] * B.n[2][1] + A.n[2][2] * B.n[2][2]);
    }

    Vector3 operator*(const Matrix3& M, const Vector3& v) {
        return Vector3(M.n[0][0] * v.x + M.n[1][0] * v.y + M.n[2][0] * v.z,
            M.n[0][1] * v.x + M.n[1][1] * v.y + M.n[2][1] * v.z,
            M.n[0][2] * v.x + M.n[1][2] * v.y + M.n[2][2] * v.z);
    }

    Vector3 Dot(const Matrix3& M, const Vector3& v) {
        return Vector3(M.n[0][0] * v.x + M.n[0][1] * v.y + M.n[0][2] * v.z,
            M.n[1][0] * v.x + M.n[1][1] * v.y + M.n[1][2] * v.z,
            M.n[2][0] * v.x + M.n[2][1] * v.y + M.n[2][2] * v.z);
    }

    inline float Dot(const Vector3& a, const Vector3& b) {
        return (a.x * b.x + a.y * b.y + a.z * b.z);
    }

    inline Vector3 Cross(const Vector3& a, const Vector3& b) {
        return Vector3(a.y * b.z - a.z * b.y,
            a.z * b.x - a.x * b.z,
            a.x * b.y - a.y * b.x);
    }

    inline Vector3 Project(const Vector3& a, const Vector3& b) {
        return (b * (Dot(a, b) / Dot(b, b)));
    }

    inline Vector3 Reject(const Vector3& a, const Vector3& b) {
        return a - (b * (Dot(a, b) / Dot(b, b)));
    }

    inline Vector3 UnitProject(const Vector3& v, const Vector3& u) {
        return (u * Dot(v, u));
    }

    inline Vector3 UnitReject(const Vector3& v, const Vector3& u) {
        return v - (u * Dot(v, u));
    }

    Matrix3 Normalize(const Matrix3& M) {
        return Matrix3(Normalize(M(0)), Normalize(M(1)), Normalize(M(2)));
    }

    Matrix3 MakeDirection(const Vector3& dir, int dirAxis, int orthgnlAxis) {
        Vector3 orthgnlAxsVctr(0, 0, 0);
        (&(orthgnlAxsVctr.x))[orthgnlAxis] = 1;
        const int* orthgnlBvctr = uBvctrs[dirAxis];
        Matrix3 M;
        M(dirAxis) = dir;
        M(orthgnlBvctr[otherUVctr[dirAxis][orthgnlAxis]]) = Normalize(Cross(dir, orthgnlAxsVctr));
        M(orthgnlAxis) = Cross(M(orthgnlBvctr[otherUVctr[dirAxis][orthgnlAxis]]), dir);
        return M;
    }

    Matrix3 MakeRotation(int axis, float radians) {
        float sine = sin(radians);
        float cosine = cos(radians);
        const int* bvctr = uBvctrs[axis];
        Matrix3 R;
        R.n[bvctr[0]][bvctr[0]] = cosine;
        R.n[bvctr[0]][bvctr[1]] = sine;
        R.n[bvctr[0]][axis] = 0;
        R.n[bvctr[1]][bvctr[0]] = -sine;
        R.n[bvctr[1]][bvctr[1]] = cosine;
        R.n[bvctr[1]][axis] = 0;
        R.n[axis][bvctr[0]] = 0;
        R.n[axis][bvctr[1]] = 0;
        R.n[axis][axis] = 1;
        return R;
    }

    Matrix3 MakeRotation(const Vector3& axis, float radians) {
        Matrix3 O = MakeDirection(axis, 2, axis.x > axis.y);
        Matrix3 R = MakeRotation(2, radians);
        return O * R * (O.T());
    }

    ostream& operator << (ostream& os, const Matrix3& M) {
        os << "{{ " << M.n[0][0] << ", " << M.n[1][0] << ", " << M.n[2][0] << " }\n";
        os << " { " << M.n[0][1] << ", " << M.n[1][1] << ", " << M.n[2][1] << " }\n";
        os << " { " << M.n[0][2] << ", " << M.n[1][2] << ", " << M.n[2][2] << " }}\n";
        return os;
    }

    struct Transform3 {
        Matrix3 R;
        Vector3 t;
        Transform3() {
            R = Matrix3();
            t = Vector3();
        }
        Transform3(float n00, float n01, float n02,
            float n10, float n11, float n12,
            float n20, float n21, float n22,
            float n30, float n31, float n32) {
            R = Matrix3(n00, n01, n02, n10, n11, n12, n20, n21, n22);
            t = Vector3(n30, n31, n32);
        }

        Transform3(const Vector3& a, const Vector3& b, const Vector3& c, const Vector3& d) {
            R = Matrix3(a, b, c);
            t = Vector3(d.x, d.y, d.z);
        }

        Transform3(const Matrix3& M, const Vector3& v) {
            R = M.clone();
            t = v.clone();
        }

        Transform3 T() {
            const Matrix3 RT = R.T();
            return Transform3(RT, RT * -t);
        }
    };

    inline Vector3 operator *(const Transform3& T, const Vector3 p) {
        return T.R * p + T.t;
    }

    inline Transform3 operator *(const Transform3& A, const Transform3& B) {
        return Transform3(A.R * B.R, A.R * B.t + A.t);
    }

    inline Transform3 operator *(const Transform3& T, const Matrix3& M) {
        return Transform3(T.R * M, T.t);
    }

    inline Transform3 operator *(const Matrix3& M, const Transform3& T) {
        return Transform3(M * T.R, M * T.t);
    }

    inline Transform3 operator +(const Transform3& T, const Vector3& v) {
        return Transform3(T.R, T.t + v);
    }

    inline Transform3 operator -(const Transform3& T, const Vector3& v) {
        return Transform3(T.R, T.t - v);
    }

    inline Vector3 GetTranslation(const Transform3& T) {
        return T.t;
    }

    inline Matrix3 GetRotation(const Transform3& T) {
        return T.R;
    }

    struct Line3 {
        Vector3 d, t;
        Line3() = default;

        Line3(float d_x, float d_y, float d_z, float t_x, float t_y, float t_z) {
            d = Vector3(d_x, d_y, d_z);
            t = Vector3(t_x, t_y, t_z);
        }

        Line3(const Vector3& direction, const Vector3& translation) {
            d = direction;
            t = translation;
        }
    };

    inline Line3 Join(Vector3 a, Vector3 b) { 
        return Line3(Normalize(b - a), a);
    }

    inline Vector3 Project(const Vector3& a, const Line3& l) {
        return l.t + l.d * Dot(a - l.t, l.d);
    }

    inline Line3 operator *(const Matrix3& R, const Line3& l) {
        return Line3(R * l.d, R * l.t);
    }

    inline Line3 operator *(const Transform3& T, const Line3& l) {
        return Line3(T.R * l.d, T * l.t);
    }

    struct Plane3 {
        float x, y, z, w;

        Plane3() = default;
        
        Plane3(const float newX, const float newY, const float newZ, const float newW) {
            x = newX; y = newY; z = newZ; w = newW;
        }

        Plane3(const Vector3& normal, const float newW) {
            x = normal.x; y = normal.y; z = normal.z; w = newW;
        }

        Plane3(const Vector3& normal, const Vector3& point) {
            x = normal.x; y = normal.y; z = normal.z; w = -Dot(normal, point);
        }

        Vector3& GetNormal() {
            return *reinterpret_cast<Vector3 *>(&x);
        }

        const Vector3& GetNormal() const {
            return *reinterpret_cast<const Vector3*>(&x);
        }

        float& operator[](int i) {
            return (&x)[i];
        }

        const float& operator[](int i) const {
            return (&x)[i];
        }
    };

    inline Plane3 operator *(const Plane3& P, float s) {
        return Plane3(P.x * s, P.y * s, P.z * s, P.w * s);
    }

    inline Plane3 operator /(const Plane3& P, float s) {
        s = 1.0F / s;
        return Plane3(P.x * s, P.y * s, P.z * s, P.w * s);
    }

    inline float Dot(const Plane3& P, const Vector3& v) {
        return P.x * v.x + P.y * v.y + P.z * v.z + P.w;
    }

    inline Vector3 GetPoint(const Plane3& P) {
        float s = -P.w;
        return Vector3(P.x * s, P.y * s, P.z * s);
    }

    inline Plane3 operator *(const Matrix3& R, const Plane3& P) {
        return Plane3(R.n[0][0] * P.x + R.n[1][0] * P.y + R.n[2][0] * P.z, 
                      R.n[0][1] * P.x + R.n[1][1] * P.y + R.n[2][1] * P.z,
                      R.n[0][2] * P.x + R.n[1][2] * P.y + R.n[2][2] * P.z,
                      P.w);
    }

    inline Plane3 operator *(const Transform3& T, const Plane3& P) {
        Vector3 normal(T.R.n[0][0] * P.x + T.R.n[1][0] * P.y + T.R.n[2][0] * P.z,
                       T.R.n[0][1] * P.x + T.R.n[1][1] * P.y + T.R.n[2][1] * P.z,
                       T.R.n[0][2] * P.x + T.R.n[1][2] * P.y + T.R.n[2][2] * P.z);
        return Plane3(normal, T.t - normal * P.w);
    }

    inline Plane3 operator +(const Plane3& P, const Vector3& t) {
        return Plane3(P.GetNormal(), t - P.GetNormal() * P.w);
    }

    inline Plane3 Join(const Vector3& a, const Vector3& b, const Vector3& c) {
        return Plane3(Normalize(Cross(b - a, c - a)), a);
    }

    inline Plane3 Join(const Vector3& a, const Line3& l) {
        return Plane3(Normalize(Cross(l.t - a, l.d)), a);
    }

    inline Vector3 Project(const Vector3& a, const Plane3& P) {
        return a - P.GetNormal() * Dot(P, a);
    }

    inline Line3 Project(const Line3& l, const Plane3& P) {
        return Line3(Normalize(UnitReject(l.d, P.GetNormal())), l.t - P.GetNormal() * Dot(P, l.d));
    }

    inline Vector3 Meet(const Line3& l, const Plane3 P) {
        return l.t - l.d * (Dot(P, l.t) / Dot(l.d, P.GetNormal()));
    }

    inline Line3 Meet(const Plane3& A, const Plane3& B) {
        const Vector3 direction = Normalize(Cross(A.GetNormal(), B.GetNormal()));
        const Vector3 pointA = GetPoint(A);
        const Vector3 disp_direction = Cross(A.GetNormal(), direction);
        return Line3(direction, pointA - disp_direction * (Dot(B, pointA) / Dot(B.GetNormal(), disp_direction)));
    }

    inline Plane3 Normalize(const Plane3& P) {
        return P / Magnitude(P.GetNormal());
    }

    ostream& operator << (ostream& os, const Plane3& P) {
        os << "{ " << P.x << ", " << P.y << ", " << P.z << ", " << P.w << " }";
        return os;
    }

    inline BaryCoords Barycentric(const Vector3& a, const Vector3& b, const Vector3& c, const Plane3 P, const Vector3& p) {
        Vector3 normal = P.GetNormal();
        normal /= Dot(Cross(b - a, c - a), normal);
        const Vector3 pa = p - a;
        const Vector3 pb = p - b;
        const Vector3 pc = p - c;
        return BaryCoords(Dot(Cross(pb, pc), normal), Dot(Cross(pc, pa), normal), Dot(Cross(pa, pb), normal));
    }

    struct AlignedBox3 {
        Vector3 minV, maxV;

        AlignedBox3() = default;

        AlignedBox3(const Vector3& newMin, const Vector3& newMax) {
            minV = newMin;
            maxV = newMax;
        }

        AlignedBox3(const Vector3* includedPoints, const int ptsLen) {
            minV = includedPoints[0];
            maxV = includedPoints[0];
            for (int i = 0; i < ptsLen; i++) {
                minV.x = min(minV.x, includedPoints[i].x);
                minV.y = min(minV.y, includedPoints[i].y);
                minV.z = min(minV.z, includedPoints[i].z);
                maxV.x = max(maxV.x, includedPoints[i].x);
                maxV.y = max(maxV.y, includedPoints[i].y);
                maxV.z = max(maxV.z, includedPoints[i].z);
            }
        }

        Vector3& operator [](int extreme) {
            return (&minV)[extreme];
        }

        const Vector3& operator[](int extreme) const {
            return (&minV)[extreme];
        }
    };

    AlignedBox3 Combine(const AlignedBox3& aBox, const Vector3& point) {
        Vector3 points[3] = { aBox.minV, aBox.maxV, point };
        return AlignedBox3(points, 3);
    }

    AlignedBox3 Combine(const AlignedBox3& aBox1, const AlignedBox3& aBox2) {
        Vector3 points[4] = { aBox1.minV, aBox1.maxV, aBox2.minV, aBox2.maxV };
        return AlignedBox3(points, 4);
    }

    AlignedBox3 Combine(const AlignedBox3* aBoxes, const int len) {
        Vector3* points = new Vector3[len * 2];
        int nextPoint = 0;
        for (int i = 0; i < len; i++) {
            points[nextPoint++] = aBoxes[i].minV;
            points[nextPoint++] = aBoxes[i].maxV;
        }
        AlignedBox3 aBox(points, len * 2);
        delete[] points;
        return aBox;
    }

    bool Contains(const AlignedBox3& abox, const Vector3& p) {
        return abox.minV.x < p.x && abox.maxV.x > p.x && 
               abox.minV.y < p.y && abox.maxV.y > p.y && 
               abox.minV.z < p.z && abox.maxV.z > p.z;
    }

    struct CornerBox3 {
        Vector3 corners[8];

        CornerBox3() = default;

        CornerBox3(const Vector3* newCorners) {
            for (int i = 0; i < 8; i++) {
                corners[i] = newCorners[i].clone();
            }
        }

        Vector3& operator [](int i) {
            return corners[i];
        }

        const Vector3& operator[](int i) const {
            return corners[i];
        }
    };

    AlignedBox3 MakeAlignedBox(const CornerBox3& cBox) {
        return AlignedBox3(cBox.corners, 8);
    }

    CornerBox3 MakeCornerBox(const AlignedBox3& aBox) {
        // fills points in order of positive rotation around the xvz plane/bivector, moving bottom to top
        CornerBox3 cBox;
        cBox.corners[0] = Vector3(aBox.maxV.x, aBox.minV.y, aBox.maxV.z);
        cBox.corners[1] = Vector3(aBox.minV.x, aBox.minV.y, aBox.maxV.z);
        cBox.corners[2] = Vector3(aBox.minV.x, aBox.minV.y, aBox.minV.z);
        cBox.corners[3] = Vector3(aBox.maxV.x, aBox.minV.y, aBox.minV.z);
        cBox.corners[4] = Vector3(aBox.maxV.x, aBox.maxV.y, aBox.maxV.z);
        cBox.corners[5] = Vector3(aBox.minV.x, aBox.maxV.y, aBox.maxV.z);
        cBox.corners[6] = Vector3(aBox.minV.x, aBox.maxV.y, aBox.minV.z);
        cBox.corners[7] = Vector3(aBox.maxV.x, aBox.maxV.y, aBox.minV.z);
        return cBox;
    }

    CornerBox3 operator *(const Transform3& T, const CornerBox3& cBox) {
        CornerBox3 tcBox = CornerBox3();
        tcBox.corners[0] = T * cBox.corners[0];
        tcBox.corners[1] = T * cBox.corners[1];
        tcBox.corners[2] = T * cBox.corners[2];
        tcBox.corners[3] = T * cBox.corners[3];
        tcBox.corners[4] = T * cBox.corners[4];
        tcBox.corners[5] = T * cBox.corners[5];
        tcBox.corners[6] = T * cBox.corners[6];
        tcBox.corners[7] = T * cBox.corners[7];
        return tcBox;
    }

    CornerBox3 operator *(const Transform3& T, const AlignedBox3& aBox) {
        return T * MakeCornerBox(aBox);
    }

    AlignedBox3 Combine(const CornerBox3* cBoxes, const int len) {
        Vector3* points = new Vector3[len * 8];
        int nextPoint = 0;
        for (int i = 0; i < len; i++) {
            points[nextPoint++] = cBoxes[i].corners[0];
            points[nextPoint++] = cBoxes[i].corners[1];
            points[nextPoint++] = cBoxes[i].corners[2];
            points[nextPoint++] = cBoxes[i].corners[3];
            points[nextPoint++] = cBoxes[i].corners[4];
            points[nextPoint++] = cBoxes[i].corners[5];
            points[nextPoint++] = cBoxes[i].corners[6];
            points[nextPoint++] = cBoxes[i].corners[7];
        }
        AlignedBox3 aBox(points, len * 8);
        delete[] points;
        return aBox;
    }

    float* Dot(const CornerBox3& cBox, const Vector3& v) { // sketchy, consider removal
        return new float[8]{ Dot(v, cBox[0]), Dot(v, cBox[1]),
                                       Dot(v, cBox[2]), Dot(v, cBox[3]),
                                       Dot(v, cBox[4]), Dot(v, cBox[5]),
                                       Dot(v, cBox[6]), Dot(v, cBox[7]) };
    }
}

namespace gmtry3DInts {
    struct IntVector3 {
        int x, y, z;
        IntVector3() = default;
        IntVector3(int i, int j, int k) {
            x = i; y = j; z = k;
        }
        IntVector3(const IntVector3& v) {
            x = v.x; y = v.y; z = v.z;
        }
        IntVector3& operator *=(int s) {
            x *= s; y *= s; z *= s;
            return *this;
        }
        IntVector3& operator /=(int s) {
            s = 1.0F / s;
            x *= s; y *= s; z *= s;
            return *this;
        }
        IntVector3& operator +=(const IntVector3& b) {
            x += b.x; y += b.y; z += b.z;
            return *this;
        }
        IntVector3& operator -=(const IntVector3& b) {
            x -= b.x; y -= b.y; z -= b.z;
            return *this;
        }
        int& operator [](int i) {
            return (&x)[i]; // or *((&x)+i)
        }
        const int& operator [](int i) const {
            return (&x)[i];
        }
    };

    inline IntVector3 operator +(const IntVector3& a, const IntVector3& b) {
        return IntVector3(a.x + b.x, a.y + b.y, a.z + b.z);
    }

    struct IntAlignedBox3 {
        IntVector3 minV, maxV;

        IntAlignedBox3() = default;

        IntAlignedBox3(const IntVector3& newMin, const IntVector3& newMax) {
            minV = newMin;
            maxV = newMax;
        }

        IntAlignedBox3(const IntVector3* includedPoints, const int ptsLen) {
            minV = includedPoints[0];
            maxV = includedPoints[0];
            for (int i = 0; i < ptsLen; i++) {
                minV.x = min(minV.x, includedPoints[i].x);
                minV.y = min(minV.y, includedPoints[i].y);
                minV.z = min(minV.z, includedPoints[i].z);
                maxV.x = max(maxV.x, includedPoints[i].x);
                maxV.y = max(maxV.y, includedPoints[i].y);
                maxV.z = max(maxV.z, includedPoints[i].z);
            }
        }

        IntVector3& operator [](int extreme) {
            return (&minV)[extreme];
        }

        const IntVector3& operator[](int extreme) const {
            return (&minV)[extreme];
        }
    };

    IntAlignedBox3 Combine(const IntAlignedBox3& aBox, const IntVector3& point) {
        IntVector3 points[3] = { aBox.minV, aBox.maxV, point };
        return IntAlignedBox3(points, 3);
    }

    IntAlignedBox3 Combine(const IntAlignedBox3& aBox1, const IntAlignedBox3& aBox2) {
        IntVector3 points[4] = { aBox1.minV, aBox1.maxV, aBox2.minV, aBox2.maxV };
        return IntAlignedBox3(points, 4);
    }

    IntAlignedBox3 Combine(const IntAlignedBox3* aBoxes, const int len) {
        IntVector3* points = new IntVector3[len * 2];
        int nextPoint = 0;
        for (int i = 0; i < len; i++) {
            points[nextPoint++] = aBoxes[i].minV;
            points[nextPoint++] = aBoxes[i].maxV;
        }
        IntAlignedBox3 aBox(points, len * 2);
        delete[] points;
        return aBox;
    }

    bool Contains(const IntAlignedBox3& abox, const IntVector3& p) {
        return abox.minV.x < p.x&& abox.maxV.x > p.x &&
            abox.minV.y < p.y&& abox.maxV.y > p.y &&
            abox.minV.z < p.z&& abox.maxV.z > p.z;
    }
}

namespace gmtry2D {

    struct Vector2 {
        float x, y;
        Vector2() = default;

        Vector2(float i, float j) {
            x = i; y = j;
        }

        Vector2(const Vector2& v) {
            x = v.x; y = v.y;
        }

        Vector2& operator *=(float s) {
            x *= s; y *= s;
            return *this;
        }

        Vector2& operator /=(float s) {
            s = 1.0F / s;
            x *= s; y *= s;
            return *this;
        }

        Vector2& operator +=(const Vector2& b) {
            x += b.x; y += b.y;
            return *this;
        }

        Vector2& operator -=(const Vector2& b) {
            x -= b.x; y -= b.y; 
            return *this;
        }

        float& operator [](int i) {
            return (&x)[i];
        }

        const float& operator [](int i) const {
            return (&x)[i];
        }
    };

    inline Vector2 MakeVector2(const gmtry3D::Vector3& v) {
        return Vector2(v.x, v.y);
    }

    inline gmtry3D::Vector3 MakeVector3(const Vector2& v, float z) {
        return gmtry3D::Vector3(v.x, v.y, z);
    }

    inline Vector2 operator +(const Vector2& a, const Vector2& b) {
        return Vector2(a.x + b.x, a.y + b.y);
    }

    inline Vector2 operator -(const Vector2& a, const Vector2& b) {
        return Vector2(a.x - b.x, a.y - b.y);
    }

    inline Vector2 operator -(const Vector2& a) {
        return Vector2(-(a.x), -(a.y));
    }

    inline Vector2 operator *(const Vector2& a, float s) {
        return Vector2(a.x * s, a.y * s);
    }

    inline Vector2 operator /(const Vector2& a, float s) {
        s = 1.0F / s;
        return Vector2(a.x * s, a.y * s);
    }

    inline float Magnitude(const Vector2& a) {
        return sqrt((a.x * a.x) + (a.y * a.y));
    }

    inline Vector2 Normalize(const Vector2& a) {
        return a / Magnitude(a);
    }

    ostream& operator << (ostream& os, const Vector2& a) {
        return os << "{ " << a.x << ", " << a.y << " }";
    }

    struct Matrix2 {
        float n[2][2];
        Matrix2() {
            n[0][0] = 1; n[0][1] = 0;
            n[1][0] = 0; n[1][1] = 1;
        }

        Matrix2(const Matrix2& M) {
            n[0][0] = M.n[0][0]; n[0][1] = M.n[0][1];
            n[1][0] = M.n[1][0]; n[1][1] = M.n[1][1];
        }

        Matrix2(float n00, float n01,
                float n10, float n11) {
            n[0][0] = n00; n[0][1] = n01;
            n[1][0] = n10; n[1][1] = n11;
        }

        Matrix2(const Vector2& a, const Vector2& b) {
            n[0][0] = a.x; n[0][1] = a.y;
            n[1][0] = b.x; n[1][1] = b.y;
        }

        float& operator ()(int i, int j) {
            return n[j][i];
        }

        const float& operator ()(int i, int j) const {
            return n[i][j];
        }

        Vector2& operator ()(int i) {
            return *(reinterpret_cast <Vector2*> (n[i]));
        }

        const Vector2& operator ()(int i) const {
            return *(reinterpret_cast <const Vector2*> (n[i]));
        }

        Matrix2 T() {
            return Matrix2(n[0][0], n[1][0],
                           n[0][1], n[1][1]);
        }
    };

    Matrix2 operator *(const Matrix2& A, const Matrix2& B) {
        return Matrix2(
            A.n[0][0] * B.n[0][0] + A.n[1][0] * B.n[0][1],
            A.n[0][1] * B.n[0][0] + A.n[1][1] * B.n[0][1],
            A.n[0][0] * B.n[1][0] + A.n[1][0] * B.n[1][1],
            A.n[0][1] * B.n[1][0] + A.n[1][1] * B.n[1][1]);
    }

    Vector2 operator*(const Matrix2& M, const Vector2& v) {
        return Vector2(M.n[0][0] * v.x + M.n[1][0] * v.y,
                       M.n[0][1] * v.x + M.n[1][1] * v.y);
    }

    Vector2 Dot(const Matrix2& M, const Vector2& v) {
        return Vector2(M.n[0][0] * v.x + M.n[0][1] * v.y,
                       M.n[1][0] * v.x + M.n[1][1] * v.y);
    }

    inline float Dot(const Vector2& a, const Vector2& b) {
        return (a.x * b.x + a.y * b.y);
    }

    inline float Cross(const Vector2& a, const Vector2& b) {
        return a.x * b.y - a.y * b.x;
    }

    inline Vector2 Project(const Vector2& a, const Vector2& b) {
        return (b * (Dot(a, b) / Dot(b, b)));
    }

    inline Vector2 Reject(const Vector2& a, const Vector2& b) {
        return a - (b * (Dot(a, b) / Dot(b, b)));
    }

    inline Vector2 UnitProject(const Vector2& v, const Vector2& u) {
        return (u * Dot(v, u));
    }

    inline Vector2 UnitReject(const Vector2& v, const Vector2& u) {
        return v - (u * Dot(v, u));
    }

    Matrix2 Normalize(const Matrix2& M) {
        return Matrix2(Normalize(M(0)), Normalize(M(1)));
    }

    Matrix2 MakeDirection(const Vector2& dir, bool dirAxis) {
        Matrix2 M;
        M(dirAxis) = dir;
        float scalar = 1 - 2 * dirAxis;
        M(!dirAxis) = Vector2(dir.y * -1.0F * scalar, dir.x * scalar);
        return M;
    }

    Matrix2 MakeRotation(float radians) {
        float sine = sin(radians);
        float cosine = cos(radians);
        return Matrix2(cosine, sine, -sine, cosine);
    }

    ostream& operator << (ostream& os, const Matrix2& M) {
        os << "{{ " << M.n[0][0] << ", " << M.n[1][0] << " }\n";
        os << " { " << M.n[0][1] << ", " << M.n[1][1] << " }}\n";
        return os;
    }

    struct Transform2 {
        Matrix2 R;
        Vector2 t;
        Transform2() {
            R = Matrix2();
            t = Vector2();
        }

        Transform2(const Transform2& T) {
            R = T.R;
            t = T.t;
        }

        Transform2(float n00, float n01,
                   float n10, float n11,
                   float n20, float n21) {
            R = Matrix2(n00, n01, n10, n11);
            t = Vector2(n20, n21);
        }

        Transform2(const Vector2& a, const Vector2& b, const Vector2& c) {
            R = Matrix2(a, b);
            t = c;
        }

        Transform2(const Matrix2& M, const Vector2& v) {
            R = M;
            t = v;
        }

        Transform2 T() {
            Matrix2 RT = R.T();
            return Transform2(RT, RT * -t);
        }
    };

    inline Vector2 operator *(const Transform2& T, const Vector2 p) {
        return T.R * p + T.t;
    }

    inline Transform2 operator *(const Transform2& A, const Transform2& B) {
        return Transform2(A.R * B.R, A.R * B.t + A.t);
    }

    inline Transform2 operator *(const Transform2& T, const Matrix2& M) {
        return Transform2(T.R * M, T.t);
    }

    inline Transform2 operator *(const Matrix2& M, const Transform2& T) {
        return Transform2(M * T.R, M * T.t);
    }

    inline Transform2 operator +(const Transform2& T, const Vector2& v) {
        return Transform2(T.R, T.t + v);
    }

    inline Transform2 operator -(const Transform2& T, const Vector2& v) {
        return Transform2(T.R, T.t - v);
    }

    inline Vector2 GetTranslation(const Transform2& T) {
        return T.t;
    }

    inline Matrix2 GetRotation(const Transform2& T) {
        return T.R;
    }

    struct Line2 {
        Vector2 d, t;
        Line2() = default;

        Line2(float d_x, float d_y, float t_x, float t_y) {
            d = Vector2(d_x, d_y);
            t = Vector2(t_x, t_y);
        }

        Line2(const Vector2& direction, const Vector2& translation) {
            d = direction;
            t = translation;
        }
    };

    inline Line2 Join(const Vector2& a, const Vector2& b) {
        return Line2(Normalize(b - a), a);
    }

    inline Vector2 Meet(const Line2& a, const Line2& b) {
        return a.t + a.d * (Cross(b.d, a.d) / Cross(a.t - b.t, b.d));
    }

    inline Vector2 Project(const Vector2& a, const Line2& l) {
        return l.d * Dot(l.t - a, l.d);
    }

    inline Line2 operator *(const Matrix2& R, const Line2& l) {
        return Line2(R * l.d, R * l.t);
    }

    inline Line2 operator *(const Transform2& T, const Line2& l) {
        return Line2(T.R * l.d, T * l.t);
    }

    inline BaryCoords Barycentric(const Vector2& a, const Vector2& b, const Vector2& c, const Vector2& p) {
        const Vector2 pa = p - a;
        const Vector2 pb = p - b;
        const Vector2 pc = p - c;
        return BaryCoords(Cross(pb, pc), Cross(pc, pa), Cross(pa, pb)) / Cross(b - a, c - a);
    }

    struct AlignedBox2 {
        Vector2 minV, maxV;

        AlignedBox2() = default;

        AlignedBox2(const Vector2& newMin, const Vector2& newMax) {
            minV = newMin;
            maxV = newMax;
        }

        AlignedBox2(const Vector2* includedPoints, const int ptsLen) {
            int ptCt = ptsLen;
            minV = includedPoints[0];
            maxV = includedPoints[0];
            for (int i = 0; i < ptCt; i++) {
                minV.x = min(minV.x, includedPoints[i].x);
                minV.y = min(minV.y, includedPoints[i].y);
                maxV.x = max(maxV.x, includedPoints[i].x);
                maxV.y = max(maxV.y, includedPoints[i].y);
            }
        }

        Vector2& operator [](int extreme) {
            return (&minV)[extreme];
        }

        const Vector2& operator[](int extreme) const {
            return (&minV)[extreme];
        }
    };

    AlignedBox2 Combine(const AlignedBox2& aBox, const Vector2& point) {
        Vector2 points[3] = { aBox.minV, aBox.maxV, point };
        return AlignedBox2(points, 3);
    }

    AlignedBox2 Combine(const AlignedBox2& aBox1, const AlignedBox2& aBox2) {
        Vector2 points[4] = { aBox1.minV, aBox1.maxV, aBox2.minV, aBox2.maxV };
        return AlignedBox2(points, 4);
    }

    AlignedBox2 Combine(const AlignedBox2* aBoxes, const int len) {
        Vector2* points = new Vector2[len * 2];
        int nextPoint = 0;
        for (int i = 0; i < len; i++) {
            points[nextPoint++] = aBoxes[i].minV;
            points[nextPoint++] = aBoxes[i].maxV;
        }
        AlignedBox2 aBox(points, len * 2);
        delete[] points;
        return aBox;
    }

    bool Contains(const AlignedBox2& abox, const Vector2& p) {
        return abox.minV.x < p.x&& abox.maxV.x > p.x &&
            abox.minV.y < p.y&& abox.maxV.y > p.y;
    }

    AlignedBox2 Constrict(const AlignedBox2& constrictee, const AlignedBox2& constrictor) {
        return AlignedBox2(Vector2(min(max(constrictee.minV.x, constrictor.minV.x), constrictor.maxV.x), 
                                   min(max(constrictee.minV.y, constrictor.minV.y), constrictor.maxV.y)), 
                           Vector2(max(min(constrictee.maxV.x, constrictor.maxV.x), constrictor.minV.x), 
                                   max(min(constrictee.maxV.y, constrictor.maxV.y), constrictor.minV.y)));
    }
}