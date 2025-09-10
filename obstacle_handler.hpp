#pragma once
#include "simulation_state.hpp"
#include "obstacle_rendering.hpp"
#include "utils.hpp"

struct ObstacleHandler {
    FluidState &fas;
    ObstacleRenderer &obstacle_renderer;

    std::vector<float> phi;

    int pencilRadius = 1;

    bool solidDrawing = false;

    float rightWallPos;
    float floorPos;

    ObstacleHandler(FluidState &fluid_attributes, ObstacleRenderer &ore): fas(fluid_attributes), obstacle_renderer(ore) {
        phi.resize(fas.gridSize);

        rightWallPos = (fas.numX - 1) * fas.cellSpacing;
        floorPos     = (fas.numY - 1) * fas.cellSpacing;
    }

    void SetUpPhi() {
        for (int i = 0; i < fas.gridSize; ++i) {
            if (fas.cellType[i] == fas.SOLID) {
                phi[i] = 1e10f;
            }
            else {
                phi[i] = 0.f;
            }
        }
    }

    void FastSweep() {
        SetUpPhi();

        const int NSweeps = 4;
        // sweep directions { start, end, step }
        const int dirX[NSweeps][3] = { {0, fas.numX - 1,  1}, {fas.numX - 1, 0, -1}, {fas.numX - 1, 0,  -1}, {0, fas.numX - 1,   1} };
        const int dirY[NSweeps][3] = { {0, fas.numY - 1, 1}, {0, fas.numY - 1, 1}, {fas.numY - 1, 0, -1}, {fas.numY - 1, 0, -1} };

        double aa[2];
        double d_new, a, b;
        int s, ix, iy, gridPos;
        const double h = fas.cellSpacing, f = 1.0;

        for (s = 0; s < NSweeps; s++) {
            for (ix = dirX[s][0]; dirX[s][2] * ix <= dirX[s][1]; ix += dirX[s][2]) {
                for (iy = dirY[s][0]; dirY[s][2] * iy <= dirY[s][1]; iy += dirY[s][2]) {
                    gridPos = ix * fas.numY + iy;
                    int left   = gridPos - fas.numY;
                    int right  = gridPos + fas.numY;
                    int top    = gridPos - 1;
                    int bottom = gridPos + 1;

                    if (fas.cellType[gridPos] == fas.SOLID) {
                        if (iy == 0 || iy == (fas.numY - 1)) {
                            if (iy == 0) {
                                aa[1] = phi[gridPos] < phi[bottom] ? phi[gridPos] : phi[bottom];
                            }
                                if (iy == (fas.numY - 1)) {
                                aa[1] = phi[top] < phi[gridPos] ? phi[top] : phi[gridPos];
                            }
                        }
                        else {
                            aa[1] = phi[top] < phi[bottom] ? phi[top] : phi[bottom];
                        }
                    
                        if (ix == 0 || ix == (fas.numX - 1)) {
                            if (ix == 0) {
                                aa[0] = phi[gridPos] < phi[right] ? phi[gridPos] : phi[right];
                            }
                            if (ix == (fas.numX - 1)) {
                                aa[0] = phi[left] < phi[gridPos] ? phi[left] : phi[gridPos];
                            }
                        }
                        else {
                            aa[0] = phi[left] < phi[right] ? phi[left] : phi[right];
                        }
                    
                        a = aa[0]; b = aa[1];
                        d_new = (fabs(a - b) < f * h ? (a + b + sqrt(2.0 * f * f * h * h - (a - b) * (a - b))) * 0.5 : std::fminf(a, b) + f * h);
                    
                        phi[gridPos] = phi[gridPos] < d_new ? phi[gridPos] : d_new;
                    }
                }
            }
        }
    }

    float samplePhi(Vector2 pos) {
        int x = pos.x;
        int y = pos.y;

        float gx = x / fas.cellSpacing - 0.5f;
        float gy = y / fas.cellSpacing - 0.5f;

        int i0 = clamp(int(std::floor(gx)), 0, fas.numX - 1);
        int j0 = clamp(int(std::floor(gy)), 0, fas.numY - 1);
        int i1 = std::min(i0 + 1, fas.numX - 1);
        int j1 = std::min(j0 + 1, fas.numY - 1);
        
        float fx = gx - i0;
        float fy = gy - j0;

        int topLeft     = i0 * fas.numY + j0;
        int topRight    = i1 * fas.numY + j0;
        int bottomLeft  = i0 * fas.numY + j1;
        int bottomRight = i1 * fas.numY + j1;

        float topLeftVal     = phi[topLeft];
        float topRightVal    = phi[topRight];
        float bottomLeftVal  = phi[bottomLeft];
        float bottomRightVal = phi[bottomRight];

        float v0 = (1 - fx) * topLeftVal    + fx * topRightVal;
        float v1 = (1 - fx) * bottomLeftVal + fx * bottomRightVal;
        return (1 - fy) * v0 + fy * v1;
    }

    Vector2 sampleGradient(Vector2 pos) {
        float eps = fas.cellSpacing * 0.25f;
        
        float ddx = (samplePhi(pos + Vector2{eps, 0}) -
                     samplePhi(pos - Vector2{eps, 0})) / (2 * eps);

        float ddy = (samplePhi(pos + Vector2{0, eps}) -
                     samplePhi(pos - Vector2{0, eps})) / (2 * eps);

        Vector2 g{ddx, ddy};
        float len = std::sqrt(g.x * g.x + g.y * g.y);
        return (len > 1e-6f) ? g / len : Vector2{0,0};
    }

    Vector2 ClosestSurfacePoint(Vector2 pos) {
        Vector2 p = {clamp(pos.x, fas.cellSpacing, rightWallPos), clamp(pos.y, fas.cellSpacing, floorPos)};
        auto phi_p = samplePhi(p);
        if (phi_p == 0) {
            return pos;
        }
        auto d = sampleGradient(p);

        float eps = 1e-1;

        for (int i = 0; i < 2; ++i) {
            float alpha = 1.f;

            for (int i = 0; i < 2; ++i) {
                auto q = p - d * alpha * phi_p;
                auto phi_q = samplePhi(q);
                if (abs(phi_q) < abs(phi_p)) {
                    p = q;
                    phi_p = phi_q;
                    d = sampleGradient(q);
                    if (abs(phi_p) < eps) {
                        return p;
                    }
                }
                else {
                    alpha = 0.7 * alpha;
                }
            }
        }

        return p;
    }

    void CollideSurfaces() {
        for (int i = 0; i < fas.num_particles; ++i) {
            int pIdx = 2 * i;
            float vx = fas.velocities[pIdx];
            float vy = fas.velocities[pIdx + 1];

            Vector2 pos = {fas.positions[pIdx], fas.positions[pIdx + 1]};

            Vector2 surface = ClosestSurfacePoint(pos);

            if (surface == pos) continue;

            float dx = surface.x - pos.x;
            float dy = surface.y - pos.y;

            float dist = sqrt(dx * dx + dy * dy);
            
            float nX, nY;
            if (dist > 0) {
                nX = dx / dist;
                nY = dy / dist;
            }
            else {
                nX = 0;
                nY = 0;
            }

            float velocityNormal = vx * -nX + vy * -nY;
            if (velocityNormal < 0) {
                vx -= velocityNormal * -nX;
                vy -= velocityNormal * -nY;
            }

            fas.positions[pIdx]     = surface.x;
            fas.positions[pIdx + 1] = surface.y;

        }
    }

    void drawSolids() {
        int localX = static_cast<int>(fas.frame_context.world_mouse_pos.x / fas.cellSpacing);
        int localY = static_cast<int>(fas.frame_context.world_mouse_pos.y / fas.cellSpacing);

        int numObstacles = fas.obstaclePositions.size();

        int numPotentialAddedObstacles = ((2 * pencilRadius + 1) * (2 * pencilRadius + 1));
        fas.obstaclePositions.resize(numObstacles + numPotentialAddedObstacles);

        int numAddedObstacles = 0;
        for (int i = -pencilRadius; i <= pencilRadius; ++i) {
            for (int j = -pencilRadius; j <= pencilRadius; ++j) {
                int x = localX + i;
                int y = localY + j;

                if (x > 0 && y > 0 && x < fas.numX - 1 && y < fas.numY - 1) {

                    int idx = x * fas.numY + y;

                    if (fas.cellType[idx] != fas.SOLID) {
                        fas.cellType[idx] = fas.SOLID;
                        fas.obstaclePositions[numObstacles + numAddedObstacles] = Vector2i{x, y};
                        numAddedObstacles++;
                    }
                }
            }
        }

        fas.obstaclePositions.resize(numObstacles + numAddedObstacles);

        if (numAddedObstacles > 0) {
            obstacle_renderer.ResizeAndUploadNewObstacleMesh(numObstacles + numAddedObstacles);
        }  

        FastSweep();
    }

    void eraseSolids() {
        int localX = static_cast<int>(fas.frame_context.world_mouse_pos.x / fas.cellSpacing);
        int localY = static_cast<int>(fas.frame_context.world_mouse_pos.y / fas.cellSpacing);

        int numObstacles = fas.obstaclePositions.size();
        int numFreedCells = 0;
        for (int i = -pencilRadius; i <= pencilRadius; ++i) {
            for (int j = -pencilRadius; j <= pencilRadius; ++j) {
                int x = localX + i;
                int y = localY + j;

                int idx = x * fas.n + y;
                if (fas.cellType[idx] == fas.SOLID && x > 0 && y > 0 && x < fas.numX - 1 && y < fas.numY - 1) {
                    fas.cellType[idx] = fas.AIR;
                    // find the index of obstacleSet with the position we want to remove
                    // replace that index with the [numObstacles - (numFreedCells + 1)] index of obstacleSet
                    // then just resize the obstacleSet and obstacleVa after this double loop
                    numFreedCells++;

                    int obstacleIdx = find(fas.obstaclePositions, Vector2i{x, y});
                    int endIdx = numObstacles - numFreedCells;

                    fas.obstaclePositions[obstacleIdx] = fas.obstaclePositions[endIdx];

                    obstacle_renderer.replaceObstacleWithEnd(obstacleIdx, endIdx);
                }
            }
        }
        obstacle_renderer.ResizeAndUploadNewObstacleMesh(numObstacles - numFreedCells);
        fas.obstaclePositions.resize(numObstacles - numFreedCells);

        FastSweep();
    }

    void render_objects() {
        if (getPencilActive()) {
            obstacle_renderer.drawPencil();
        }
    }

    bool getPencilActive() {
        return this->solidDrawing;
    }

    int getPencilRadius() {
        return this->pencilRadius;
    }

    void addToPencilRadius(int add) {
        this->pencilRadius += add;
        obstacle_renderer.changePencilRadius(pencilRadius);
    }

    void setPencilRadius(int set) {
        this->pencilRadius = set;
        obstacle_renderer.changePencilRadius(pencilRadius);
    }

    void setSolidDrawer(bool set) {
        this->solidDrawing = set;
    }

};