#ifndef UNICODE
#define UNICODE
#endif 

#include <windows.h>
#include <windowsx.h>
#include <string>

extern "C" {
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
}
#include "systems.hpp"

#define PIXEL_DEFAULT_SIZE (sizeof(int)*8)
#define RENDERDIST 0x100
#define RENDERRAYS (false)
#define RENDERCUBES (true)
#define RENDERCONTAINERS (false)
bool RENDERMODE = true;

bool running;
int* mainpxlbuf;
int mainwidth;
int mainheight;
BITMAPINFO mainbmi;
BITMAPINFOHEADER mainbmih;

int mousex, mousey;
float mousesensitivity = 0.0005;

unsigned int* vertexcolors = 0;
Vector3* vertices = 0;
float* certainties = 0;
int vwidth, cwidth, vheight, batch_width;

PlaneContainer* plane_map = 0;
Vector3 pm_disp;
int log2_pmwidth;

Vector3 vpoints[4] = { Vector3(-0.5, -0.5, 4), Vector3(0.5, -0.5, 4), Vector3(0.5, 0.5, 4), Vector3(-0.5, 0.5, 4) };
Vector3 vcolors[4] = { Vector3(1, 0, 0), Vector3(0, 1, 0), Vector3(0, 0, 1), Vector3(1, 1, 1) };

float* mainzbuf;
float farplane = 1000.0F;
float nearplane = 0.125F;
float double_tangent = 1.1547F;
Transform3 camtfrm = Transform3();
Transform3 tfrm = Transform3();
Transform2 screentfrm = Transform2();
AlignedBox2 screenbox;


void ResetPixelBuffer(int newwidth, int newheight, int newcolor) {
    mainwidth = newwidth;
    mainheight = newheight;
    screenbox = AlignedBox2(Vector2(), Vector2(mainwidth - 1, mainheight - 1));
    screentfrm = Transform2(1.0F, 0.0F, 0.0F, -1.0F, mainwidth / 2.0F, mainheight / 2.0F);
    delete[] mainpxlbuf;
    delete[] mainzbuf;
    int bufarea = mainwidth * mainheight;
    mainpxlbuf = new int[bufarea];
    mainzbuf = new float[bufarea];
    for (int i = 0; i < bufarea; i++) {
        mainpxlbuf[i] = newcolor;
        mainzbuf[i] = farplane;
    }

    mainbmih.biWidth = mainwidth;
    mainbmih.biHeight = -mainheight;
    mainbmi.bmiHeader = mainbmih;
}

void ClearPixels() {
    int bufarea = mainwidth * mainheight;
    for (int i = 0; i < bufarea; i++) {
        mainpxlbuf[i] = 0;
        mainzbuf[i] = farplane;
    }
}

void DrawPixels(HWND hwnd) {
    SetDIBitsToDevice(GetDC(hwnd),
        0,
        0,
        mainwidth,
        mainheight,
        0,
        0,
        0,
        mainheight,
        mainpxlbuf,
        &mainbmi,
        DIB_RGB_COLORS);
}

void WriteToPixelBuffer(int* pxlbuf, int newwidth, int newheight) {
    int minwidth = mainwidth + (newwidth < mainwidth) * (newwidth - mainwidth);
    int minheight = mainheight + (newheight < mainheight) * (newheight - mainheight);
    for (int x = 0; x < minwidth; x++) {
        for (int y = 0; y < minheight; y++) {
            mainpxlbuf[x + y * mainwidth] = pxlbuf[x + y * newwidth];
        }
    }
}

void WriteToPixelBuffer(int* pxlbuf, int newwidth, int newheight, int x_disp, int y_disp) {
    if (x_disp >= mainwidth || y_disp >= mainheight || (x_disp + newwidth) <= 0 || (y_disp + newheight) <= 0) return;
    int minx = (x_disp > 0) * x_disp;
    int maxx = mainwidth + ((x_disp + newwidth) < mainwidth) * ((x_disp + newwidth) - mainwidth);
    int miny = (y_disp > 0) * y_disp;
    int maxy = mainheight + ((y_disp + newheight) < mainheight) * ((y_disp + newheight) - mainheight);
    for (int x = minx; x < maxx; x++) {
        for (int y = miny; y < maxy; y++) {
            mainpxlbuf[x + y * mainwidth] = pxlbuf[(x - x_disp) + (y - y_disp) * newwidth];
        }
    }
}

inline void WriteToPixelBuffer(unsigned int color, int x, int y) {
    if (x >= 0 && x < mainwidth && y >= 0 && y < mainheight) mainpxlbuf[x + y * mainwidth] = color;
}

inline void WriteToPixelBuffer(float R, float G, float B, int x, int y, float z) {
    int r, g, b;
    r = 255.0F * R; g = 255.0F * G; b = 255.0F * B;
    r += (r > 255) * (255 - r);
    g += (g > 255) * (255 - g);
    b += (b > 255) * (255 - b);
    if (z < mainzbuf[x + y * mainwidth] && z > nearplane) {
        mainpxlbuf[x + y * mainwidth] = (r << 16) + (g << 8) + b;
        mainzbuf[x + y * mainwidth] = z;
    }
}

inline void WriteToPixelBuffer(const Vector3& color, int x, int y, float z) {
    int r, g, b;
    r = 255.0F * color.x; g = 255.0F * color.y; b = 255.0F * color.z;
    r += (r > 255) * (255 - r);
    g += (g > 255) * (255 - g);
    b += (b > 255) * (255 - b);
    if (z <= mainzbuf[x + y * mainwidth] && z > nearplane) {
        mainpxlbuf[x + y * mainwidth] = (r << 16) + (g << 8) + b;
        mainzbuf[x + y * mainwidth] = z;
    }
}

void WriteToPixelBuffer(float* pxlbuf, int newwidth, int newheight) {
    int minwidth = mainwidth + (newwidth < mainwidth) * (newwidth - mainwidth);
    int minheight = mainheight + (newheight < mainheight) * (newheight - mainheight);
    int r, g, b;
    for (int x = 0; x < minwidth; x++) {
        for (int y = 0; y < minheight; y++) {
            r = 255 * pxlbuf[(x + y * newwidth) * 4];
            r += (r > 255) * (255 - r);
            g = 255 * pxlbuf[(x + y * newwidth) * 4 + 1];
            g += (g > 255) * (255 - g);
            b = 255 * pxlbuf[(x + y * newwidth) * 4 + 2];
            b += (b > 255) * (255 - b);
            mainpxlbuf[x + y * mainwidth] = (r << 16) + (g << 8) + b;
        }
    }
}

inline Vector3 VectorizeColor(unsigned int color) {
    return Vector3((color >> 16) & 0b11111111, (color >> 8) & 0b11111111, color & 0b11111111) / 255.0F;
}

void RBSwap(unsigned int* pxlbuf, int length) {
    for (int i = 0; i < length; i++)
        pxlbuf[i] =
        ((pxlbuf[i] & 0b11111111) << 16) +
        (pxlbuf[i] & (0b11111111 << 8)) +
        ((pxlbuf[i] & (0b11111111 << 16)) >> 16);
}

void RenderTriangle(const Vector3& a, const Vector3& b, const Vector3& c, 
    const Vector3& colora, const Vector3& colorb, const Vector3& colorc) {
    Vector3 tfrma = tfrm * a;
    Vector3 tfrmb = tfrm * b;
    Vector3 tfrmc = tfrm * c;
    bool exclusions[3] = { tfrma.z < nearplane, tfrmb.z < nearplane, tfrmc.z < nearplane };
    
    Plane3 tfrmplane = Join(tfrma, tfrmb, tfrmc);
    if (tfrmplane.w == 0.0F || (exclusions[0] && exclusions[1] && exclusions[2])) return;
    float scale = (mainwidth + (mainheight > mainwidth) * (mainheight - mainwidth)) / double_tangent;
    float invscale = 1.0F / scale;

    Transform2 screentfrm(1.0F, 0.0F, 0.0F, -1.0F, mainwidth / 2.0F, mainheight / 2.0F);
    Transform2 invscreentfrm = screentfrm.T();

    AlignedBox2 scanbox;
    if (exclusions[0] || exclusions[1] || exclusions[2]) scanbox = screenbox;
    else {
        Vector2 proj_points[3] = {
            screentfrm * (Vector2(tfrma.x, tfrma.y) * scale / tfrma.z),
            screentfrm * (Vector2(tfrmb.x, tfrmb.y) * scale / tfrmb.z),
            screentfrm * (Vector2(tfrmc.x, tfrmc.y) * scale / tfrmc.z)
        };
        scanbox = Constrict(AlignedBox2(proj_points, 3), screenbox);
    }
    int minx = scanbox.minV.x;
    int maxx = scanbox.maxV.x + 1;
    int miny = scanbox.minV.y;
    int maxy = scanbox.maxV.y + 1;

    Vector2 point;
    Vector3 colorp;
    BaryCoords bary;
    float zp;
    
    for (int y = miny; y < maxy; y++) {
        for (int x = minx; x < maxx; x++) {
            point = (invscreentfrm * (Vector2(x, y))) * invscale;
            zp = -tfrmplane.w / (point.x * tfrmplane.x + point.y * tfrmplane.y + tfrmplane.z);
            bary = Barycentric(tfrma, tfrmb, tfrmc, tfrmplane, Vector3(point.x * zp, point.y * zp, zp));
            colorp = colora * bary.a + colorb * bary.b + colorc * bary.c;
            if (bary.a >= 0.01 && bary.b >= 0.01 && bary.c >= 0.01) WriteToPixelBuffer(colorp.x, colorp.y, colorp.z, x, y, zp);
            else if (bary.a >= 0 && bary.b >= 0 && bary.c >= 0) WriteToPixelBuffer(1.0F, 1.0F, 1.0F, x, y, zp);
        }
    }
}

void RenderQuad(const Vector3& a, const Vector3& b, const Vector3& c, const Vector3& d,
    const Vector3& colora, const Vector3& colorb, const Vector3& colorc, const Vector3& colord) {
    RenderTriangle(a, b, c, colora, colorb, colorc);
    RenderTriangle(c, d, a, colorc, colord, colora);
}

inline Vector3 faoijf(const Vector3& a, const Vector3& b, int dim) {
    // a is target
    if (a[dim] < screenbox[0][dim] && b[dim] >= screenbox[1][dim])
        return b + (a - b) * ((b[dim] - screenbox[0][dim]) / (b[dim] - a[dim]));
    if (a[dim] > screenbox[1][dim] && b[dim] <= screenbox[1][dim])
        return b + (a - b) * ((screenbox[1][dim] - b[dim]) / (a[dim] - b[dim]));
    return a;
}

void RenderLine2D(const Vector3& outera, const Vector3& outerb) {
    Vector3 a = faoijf(faoijf(outera, outerb, 0), outerb, 1);
    Vector3 b = faoijf(faoijf(outerb, outera, 0), outera, 1);
    float scanlength = Magnitude(MakeVector2(b - a));
    Vector3 direction = (b - a) / scanlength;

    Vector3 point = a;
    for (int l = 0; l < scanlength; l++) {
        point += direction;
        if (Contains(screenbox, MakeVector2(point))) WriteToPixelBuffer(1.0F, 1.0F, 1.0F, point.x, point.y, 1.0F / point.z);
    }
}

void RenderLine(const Vector3& a, const Vector3& b) {
    Vector3 projA = tfrm * a;
    Vector3 projB = tfrm * b;
    if (projA.z < nearplane && projB.z < nearplane) return;
    if (projA.z < nearplane) projA = projB + (projA - projB) * ((projB.z - nearplane) / (projB.z - projA.z));
    if (projB.z < nearplane) projB = projA + (projB - projA) * ((projA.z - nearplane) / (projA.z - projB.z));
    float scale = (mainwidth + (mainheight > mainwidth) * (mainheight - mainwidth)) / double_tangent;
    float invzA = 1.0F / projA.z; float invzB = 1.0F / projB.z;
    RenderLine2D(MakeVector3(screentfrm * (Vector2(projA.x, projA.y) * scale * invzA), invzA),
        MakeVector3(screentfrm * (Vector2(projB.x, projB.y) * scale * invzB), invzB));
}

void RenderAlignedCube(const Vector3& offset, float width) {
    RenderLine(offset + Vector3(0, 0, 0) * width, offset + Vector3(1, 0, 0) * width);
    RenderLine(offset + Vector3(1, 0, 0) * width, offset + Vector3(1, 0, 1) * width);
    RenderLine(offset + Vector3(1, 0, 1) * width, offset + Vector3(0, 0, 1) * width);
    RenderLine(offset + Vector3(0, 0, 1) * width, offset + Vector3(0, 0, 0) * width);

    RenderLine(offset + Vector3(0, 0, 0) * width, offset + Vector3(0, 1, 0) * width);
    RenderLine(offset + Vector3(1, 0, 0) * width, offset + Vector3(1, 1, 0) * width);
    RenderLine(offset + Vector3(1, 0, 1) * width, offset + Vector3(1, 1, 1) * width);
    RenderLine(offset + Vector3(0, 0, 1) * width, offset + Vector3(0, 1, 1) * width);

    RenderLine(offset + Vector3(0, 1, 0) * width, offset + Vector3(1, 1, 0) * width);
    RenderLine(offset + Vector3(1, 1, 0) * width, offset + Vector3(1, 1, 1) * width);
    RenderLine(offset + Vector3(1, 1, 1) * width, offset + Vector3(0, 1, 1) * width);
    RenderLine(offset + Vector3(0, 1, 1) * width, offset + Vector3(0, 1, 0) * width);
}

void RenderPlane(const PlaneCube* pc, int plane_idx, const Vector3& pc_disp) {
    int i = (plane_idx & 1); int k = (plane_idx & 2) >> 1; int j = (plane_idx & 4) >> 2;
    Vector3 plane_disp = pc_disp + Vector3(i, j, k);
    AlignedBox3 planebox = AlignedBox3(plane_disp, plane_disp + Vector3(1, 1, 1));
    CornerBox3 tfrmpts = tfrm * planebox;
    AlignedBox3 tfrmbox = MakeAlignedBox(tfrmpts);

    if ((Dot(tfrm * pc->n[i][j][k], pc->p[i][j][k]) == 0.0F) || (tfrmbox.maxV.z < nearplane)) return;
    float scale = (mainwidth + (mainheight > mainwidth) * (mainheight - mainwidth)) / double_tangent;
    float invscale = 1.0F / scale;

    Transform2 invscreentfrm = screentfrm.T();

    AlignedBox2 scanbox;

    if (tfrmbox.minV.z < nearplane) scanbox = screenbox;
    else {
        Vector2 proj_points[8] = {
            screentfrm * (Vector2(tfrmpts.corners[0].x, tfrmpts.corners[0].y) * scale / tfrmpts.corners[0].z),
            screentfrm * (Vector2(tfrmpts.corners[1].x, tfrmpts.corners[1].y) * scale / tfrmpts.corners[1].z),
            screentfrm * (Vector2(tfrmpts.corners[2].x, tfrmpts.corners[2].y) * scale / tfrmpts.corners[2].z),
            screentfrm * (Vector2(tfrmpts.corners[3].x, tfrmpts.corners[3].y) * scale / tfrmpts.corners[3].z),
            screentfrm * (Vector2(tfrmpts.corners[4].x, tfrmpts.corners[4].y) * scale / tfrmpts.corners[4].z),
            screentfrm * (Vector2(tfrmpts.corners[5].x, tfrmpts.corners[5].y) * scale / tfrmpts.corners[5].z),
            screentfrm * (Vector2(tfrmpts.corners[6].x, tfrmpts.corners[6].y) * scale / tfrmpts.corners[6].z),
            screentfrm * (Vector2(tfrmpts.corners[7].x, tfrmpts.corners[7].y) * scale / tfrmpts.corners[7].z)
        };
        scanbox = Constrict(AlignedBox2(proj_points, 8), screenbox);
    }
    int minx = scanbox.minV.x;
    int maxx = scanbox.maxV.x + 1;
    int miny = scanbox.minV.y;
    int maxy = scanbox.maxV.y + 1;

    Vector3 point;
    Line3 ray;
    Vector3 colorp;
    float zp;

    for (int y = miny; y < maxy; y++) {
        for (int x = minx; x < maxx; x++) {
            ray = Line3(camtfrm.R * MakeVector3((invscreentfrm * (Vector2(x, y))) * invscale, 1), camtfrm.t);
            point = Meet(ray, Plane3(pc->n[i][j][k], pc->p[i][j][k]));
            zp = Dot(point - camtfrm.t, camtfrm.R(2));
            colorp = InterpolateHue(HueCoefficients(point - plane_disp), *pc, i, j, k) * pc->c[i][j][k];
            if (Contains(planebox, point)) WriteToPixelBuffer(colorp.x, colorp.y, colorp.z, x, y, zp);
        }
    }

    if (RENDERCUBES) RenderAlignedCube(plane_disp, 1);
}

void RenderPlaneContainer(const PlaneContainer* PC, int log2_width, const Vector3& PC_disp) {
    int halfwidth = 1 << (log2_width - 1);
    if (RENDERCONTAINERS && log2_width == 1) RenderAlignedCube(PC_disp, 1 << log2_width);

    if (log2_width > 1) {
        for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
        for (int k = 0; k < 2; k++)
        if (PC->sub[i][j][k])
        if (PCHasContents(*(PC->sub[i][j][k])))
        RenderPlaneContainer(PC->sub[i][j][k], log2_width - 1, PC_disp + Vector3(halfwidth * i, halfwidth * j, halfwidth * k));
    }
    else if (PC->planes) {
        for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
        for (int k = 0; k < 2; k++)
        if (PC->planes->c[i][j][k])
        RenderPlane(PC->planes, i + 2 * j + 4 * k, PC_disp);
    }
}

void RenderTheThing(HWND hwnd) {
    if (RENDERMODE) RenderPlaneContainer(plane_map, log2_pmwidth, pm_disp);
    else {
        int vidx1, cidx1, vidx2, cidx2;
        for (int y = 0; y < vheight - 1; y++) {
            for (int x = 0; x < vwidth - 1; x++) {

                std::string stringabingaling = std::to_string(x) + ", " + std::to_string(y);
                SetWindowTextA(hwnd, stringabingaling.c_str());

                float color_scalar = certainties[x + y * vwidth];

                vidx1 = x + y * vwidth;
                cidx1 = x + y * cwidth;
                vidx2 = vidx1 + vwidth;
                cidx2 = cidx1 + cwidth;
                RenderTriangle(vertices[vidx1], vertices[vidx2], vertices[vidx2 + 1],
                    VectorizeColor(vertexcolors[cidx1]) * color_scalar,
                    VectorizeColor(vertexcolors[cidx2]) * color_scalar,
                    VectorizeColor(vertexcolors[cidx2 + 1]) * color_scalar);
                RenderTriangle(vertices[vidx2 + 1], vertices[vidx1 + 1], vertices[vidx1],
                    VectorizeColor(vertexcolors[cidx2 + 1]) * color_scalar,
                    VectorizeColor(vertexcolors[cidx1 + 1]) * color_scalar,
                    VectorizeColor(vertexcolors[cidx1]) * color_scalar);
                if (RENDERRAYS) RenderLine(Vector3(), vertices[vidx1]);
            }
        }
    }
    RenderTriangle(vpoints[0], vpoints[1], vpoints[2], vcolors[0], vcolors[1], vcolors[2]);
    RenderTriangle(vpoints[2], vpoints[3], vpoints[0], vcolors[2], vcolors[3], vcolors[0]);
}

LRESULT CALLBACK WindowProc(HWND hwnd, UINT uMsg, WPARAM wParam, LPARAM lParam)
{
    int newmousex, newmousey;
    switch (uMsg) {
    case WM_DESTROY:
        running = false;
        PostQuitMessage(0);
        return 0;

    case WM_CREATE:
        mainbmih.biSize = sizeof(BITMAPINFOHEADER);
        mainbmih.biPlanes = 1;
        mainbmih.biBitCount = PIXEL_DEFAULT_SIZE;
        mainbmih.biCompression = BI_RGB;
        mainbmih.biSizeImage = NULL;
        mainbmih.biXPelsPerMeter = NULL; // ?
        mainbmih.biYPelsPerMeter = NULL; // ?
        mainbmih.biClrUsed = NULL;
        mainbmih.biClrImportant = NULL;

        RECT wndrect;
        GetWindowRect(hwnd, &wndrect);
        ResetPixelBuffer(wndrect.right - wndrect.left, wndrect.bottom - wndrect.top, 0xFFp16);
        running = true;


        return 0;

    case WM_SIZE:
        ResetPixelBuffer(LOWORD(lParam), HIWORD(lParam), 0x00p16);
        return 0;

    case WM_LBUTTONDOWN:
        mousex = GET_X_LPARAM(lParam);
        mousey = GET_Y_LPARAM(lParam);

        return 0;

    case WM_LBUTTONUP:
        newmousex = GET_X_LPARAM(lParam);
        newmousey = GET_Y_LPARAM(lParam);

        camtfrm.R = MakeRotation(1, (newmousex - mousex) * mousesensitivity) * camtfrm.R;
        camtfrm.R = MakeRotation(0, (newmousey - mousey) * mousesensitivity) * camtfrm.R;

        mousex = newmousex;
        mousey = newmousey;

        return 0;

    case WM_CHAR:
        std::string stringabingaling = std::to_string(wParam);
        //SetWindowTextA(hwnd, stringabingaling.c_str());
        switch (wParam) {
            case 'w':
                camtfrm = camtfrm + camtfrm.R(2);
                break;

            case 'a':
                camtfrm = camtfrm - camtfrm.R(0);
                break;

            case 's':
                camtfrm = camtfrm - camtfrm.R(2);
                break;

            case 'd':
                camtfrm = camtfrm + camtfrm.R(0);
                break;

            case 'c':
                camtfrm = camtfrm + camtfrm.R(1);
                break;

            case 'x':
                camtfrm = camtfrm - camtfrm.R(1);
                break;

            case 't':
                RENDERMODE = !RENDERMODE;
                break;

            case ' ':
                tfrm = camtfrm.T();
                ClearPixels();
                RenderTheThing(hwnd);
                DrawPixels(hwnd);
                break;
        }
        stringabingaling += " { " + 
            std::to_string(camtfrm.t.x) + ", " +
            std::to_string(camtfrm.t.y) + ", " +
            std::to_string(camtfrm.t.z) + " } ";
        SetWindowTextA(hwnd, stringabingaling.c_str());
        return 0;
    }
    return DefWindowProc(hwnd, uMsg, wParam, lParam);
}

int WINAPI wWinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, PWSTR pCmdLine, int nCmdShow)
{
    // Register the window class.
    const wchar_t CLASS_NAME[] = L"Sample Window Class";

    WNDCLASS wc = { };

    wc.lpfnWndProc = WindowProc;
    wc.hInstance = hInstance;
    wc.lpszClassName = CLASS_NAME;

    RegisterClass(&wc);

    // Create the window.

    HWND hwnd = CreateWindowEx(
        0,                              // Optional window styles.
        CLASS_NAME,                     // Window class
        L"Learn to Program Windows",    // Window text
        WS_OVERLAPPEDWINDOW,            // Window style

        // Size and position
        CW_USEDEFAULT, CW_USEDEFAULT, CW_USEDEFAULT, CW_USEDEFAULT,

        NULL,       // Parent window
        NULL,       // Menu
        hInstance,  // Instance handle
        NULL        // Additional application data
    );

    if (hwnd == NULL)
    {
        return 0;
    }

    ShowWindow(hwnd, nCmdShow);

    /*
    */
    int x, y, n;
    unsigned int* image1 = reinterpret_cast<unsigned int*>(stbi_load("MCImg1.png", &x, &y, &n, 4));
    unsigned int* image2 = reinterpret_cast<unsigned int*>(stbi_load("MCImg2.png", &x, &y, &n, 4));

    Gen_State gen = Gen_State();
    Cam_State cam(x, y - 128 * 6, 0.57735F, 2);
    DP_State dps = DP_State(&gen, &cam, 8, 1, 1);
    IFE_State ife(&gen, &dps);
    PMU_State pmu(&ife);

    dps.SetImages(&(image1[x * (256)]), &(image2[x * (256)]));

    dps.Execute();

    ife.Execute();
    for (int i = 0; i < 50; i++) pmu.Execute();

    vertexcolors = GetVColors(&dps);
    vertices = GetVertices(&dps);
    certainties = GetCertainties(&dps);
    RBSwap(vertexcolors, (dps.shifts_dims[0] + 1) * dps.shifts_dims[1]);

    vwidth = dps.shifts_dims[0];
    cwidth = dps.vcolors_dims[0];
    vheight = dps.shifts_dims[1];
    batch_width = 1 << dps.log2_bwidth;

    for (int i = 0; i < vwidth * vheight; i++) { certainties[i] *= 0.1; if (certainties[i] > 1.0F) certainties[i] = 1.0F; }

    plane_map = pmu.plane_map;
    log2_pmwidth = pmu.log2_pmwidth;
    pm_disp = pmu.pm_disp;

    stbi_image_free(image1);
    stbi_image_free(image2);

    RenderTheThing(hwnd);

    DrawPixels(hwnd);

    MSG msg = { };
    while (GetMessage(&msg, NULL, 0, 0))
    {
        TranslateMessage(&msg);
        DispatchMessage(&msg);
    }
    running = false;

    delete[] vertexcolors;
    delete[] vertices;
    delete[] certainties;

    delete[] mainpxlbuf;
    delete[] mainzbuf;

    return 0;
}