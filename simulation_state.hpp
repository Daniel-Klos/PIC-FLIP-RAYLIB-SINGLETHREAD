#pragma once
#include <raylib.h>
#include <vector>
#include <cmath>

#include "collision_grid.hpp"
#include "frame_context.hpp"
#include "utils.hpp"

#include <iostream>

struct FluidState {
    int num_particles;
    //std::vector<float> masses;
    std::vector<float> positions;
    std::vector<float> renderPositions;
    std::vector<float> velocities;
    std::vector<int> cellType;
    std::vector<float> u;
    std::vector<float> v;
    std::vector<float> sumUGridWeights;
    std::vector<float> sumVGridWeights;
    std::vector<float> prevU;
    std::vector<float> prevV;
    std::vector<float> cellDensities;
    //std::vector<float> cellMasses
    std::vector<int> fluid_cells;
    std::vector<Vector2i> obstaclePositions;
    std::vector<float> phi;
    std::vector<int> particleAges;

    float cellSpacing;
    float halfSpacing;
    float invSpacing;
    float radius;
    int numX;
    int numY;
    int n; // same as numY
    int gridSize;
    int num_fluid_cells;
    int age_constant = 5;
    float rightWallPos;
    float floorPos;

    float particleRestDensity = 0.f; // INITIALIZE TO 0.f

    float vorticityStrength;
    float flipRatio;

    float gravityX;
    float gravityY;

    bool stop = false;
    bool step = false; 

    bool fireActive = false;
    std::vector<float> temperatures;
    float groundConductivity = 20000.f;  // how quickly the ground transfers heat to the particles
    float interConductivity = 10000.f;    // how quickly particles transfer heat between themselves
    float fireStrength = 30.f;         // how quickly particles accelerate upwards due to heat
    float tempDiffusion = 50.f;        // how quickly particles lose heat

    //std::vector<bool> debug; // for any debug condition you want, just set pIdx = 1

    int FLUID = 0;
    int AIR = 1;
    int SOLID = 2;

    FrameContext frame_context;

    CollisionGrid cellOccupants;

    FluidState(int num_particles_, int numX_, float vorticityStrength_, float flipRatio_, float gravityX_, float gravityY_, float simWIDTH, float simHEIGHT): num_particles(num_particles_), numX(numX_), vorticityStrength(vorticityStrength_), flipRatio(flipRatio_), gravityX(gravityX_), gravityY(gravityY_), frame_context(simWIDTH, simHEIGHT) {

        cellSpacing = frame_context.WIDTH / numX;
        halfSpacing = cellSpacing * 0.5;
        invSpacing = 1.f / cellSpacing;
        radius = 0.3 * cellSpacing;

        numY = std::floor(frame_context.HEIGHT / cellSpacing);
        n = numY;
        gridSize = numX * numY;

        rightWallPos = (numX - 1) * cellSpacing;
        floorPos     = (numY - 1) * cellSpacing;

        particleAges.resize(num_particles);
        std::fill(begin(particleAges), end(particleAges), 0);

        cellOccupants = CollisionGrid(numX, numY);

        positions.resize(2 * num_particles);
        renderPositions.resize(2 * num_particles);
        velocities.resize(2 * num_particles);
        temperatures.resize(num_particles);
        cellType.resize(gridSize);
        u.resize(gridSize);
        v.resize(gridSize);
        sumUGridWeights.resize(gridSize);
        sumVGridWeights.resize(gridSize);
        prevU.resize(gridSize);
        prevV.resize(gridSize);
        cellDensities.resize(gridSize);
        fluid_cells.resize(gridSize);

        //debug.resize(num_particles);
        
        // initializing particle positions
        float separation = 2.1; 

        int wideNum = std::floor((frame_context.WIDTH - 2 * cellSpacing - 2) / (radius * 2.1));
        int highNum = num_particles / wideNum;

        float starting_px = radius + cellSpacing + 2;
        float starting_py = (frame_context.HEIGHT - (radius * separation * highNum)) / 2 + radius;

        float px = starting_px;
        float py = starting_py;

        int addTo = num_particles - (wideNum * highNum);

        bool offset = true;
        for (int i = 0; i < wideNum * highNum + addTo; ++i) {
            this->positions[i * 2] = px;
            this->positions[i * 2 + 1] = py;

            px += this->radius * separation;

            if ((i + 1) % wideNum == 0) {
                px = starting_px;
                if (offset) {
                    px += this->radius;
                }
                py += this->radius * separation;
                offset = !offset;
            }
        }
    }

    float Weight(float px, float py, float gx, float gy) {
        gx *= cellSpacing;
        gy *= cellSpacing;

        float dx = 1.f - std::abs((px - gx)) * invSpacing;
        float dy = 1.f - std::abs((py - gy)) * invSpacing;

        return dx * dy;
    }

    void WeightGradientFD(float px, float py, float gx, float gy, float &gradX, float &gradY) {
        float eps = cellSpacing * 0.001f;

        gradX = (Weight(px + eps, py, gx, gy) - Weight(px - eps, py, gx, gy)) / (2 * eps);
        gradY = (Weight(px, py + eps, gx, gy) - Weight(px, py - eps, gx, gy)) / (2 * eps);
    }

    void WeightGradientTest(float px, float py, float gx, float gy, float &gradX, float &gradY) {
        
    }

    void WeightGradient(float px, float py, float gx, float gy, float &gradx, float &grady) {
        gx *= cellSpacing;
        gy *= cellSpacing;

        float dx = (px - gx) * invSpacing;
        float dy = (py - gy) * invSpacing;

        gradx = -sign(dx) * (1.f - abs(dy)) * invSpacing;
        grady = -sign(dy) * (1.f - abs(dx)) * invSpacing;
    }

    float samplePhi(Vector2 pos) {
        int x = pos.x;
        int y = pos.y;

        float gx = x / cellSpacing - 0.5f;
        float gy = y / cellSpacing - 0.5f;

        int i0 = clamp(int(std::floor(gx)), 0, numX - 1);
        int j0 = clamp(int(std::floor(gy)), 0, numY - 1);
        int i1 = std::min(i0 + 1, numX - 1);
        int j1 = std::min(j0 + 1, numY - 1);
        
        float fx = gx - i0;
        float fy = gy - j0;

        int topLeft     = i0 * numY + j0;
        int topRight    = i1 * numY + j0;
        int bottomLeft  = i0 * numY + j1;
        int bottomRight = i1 * numY + j1;

        float topLeftVal     = phi[topLeft];
        float topRightVal    = phi[topRight];
        float bottomLeftVal  = phi[bottomLeft];
        float bottomRightVal = phi[bottomRight];

        float v0 = (1 - fx) * topLeftVal    + fx * topRightVal;
        float v1 = (1 - fx) * bottomLeftVal + fx * bottomRightVal;
        return (1 - fy) * v0 + fy * v1;
    }
    
    Vector2 sampleGradient(Vector2 pos) {
        float eps = cellSpacing * 0.25f;
        
        float ddx = (samplePhi(pos + Vector2{eps, 0}) -
                     samplePhi(pos - Vector2{eps, 0})) / (2 * eps);

        float ddy = (samplePhi(pos + Vector2{0, eps}) -
                     samplePhi(pos - Vector2{0, eps})) / (2 * eps);

        Vector2 g{ddx, ddy};
        float len = std::sqrt(g.x * g.x + g.y * g.y);
        return (len > 1e-6f) ? g / len : Vector2{0,0};
    }

    Vector2 ClosestSurfacePoint(Vector2 pos) {
        Vector2 p = {clamp(pos.x, cellSpacing, rightWallPos), clamp(pos.y, cellSpacing, floorPos)};
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
        for (int i = 0; i < num_particles; ++i) {
            int pIdx = 2 * i;
            float vx = velocities[pIdx];
            float vy = velocities[pIdx + 1];

            Vector2 pos = {positions[pIdx], positions[pIdx + 1]};

            Vector2 surface = ClosestSurfacePoint(pos);

            if (surface == pos) continue;

            float dx = surface.x - pos.x;
            float dy = surface.y - pos.y;

            float dist = sqrt(dx * dx + dy * dy);
            
            float nX, nY;
            
            nX = dx / dist;
            nY = dy / dist;

            float velocityNormal = -(vx * nX + vy * nY);
            if (velocityNormal < 0) {
                vx += velocityNormal * nX;
                vy += velocityNormal * nY;
            }

            positions[pIdx]     = surface.x;
            positions[pIdx + 1] = surface.y;
        }
    }

    void FillCellOccupants() {
        cellOccupants.clear();

        const float minX = cellSpacing;
        const float maxX = frame_context.WIDTH - cellSpacing;
        const float minY = cellSpacing;
        const float maxY = frame_context.HEIGHT - cellSpacing;

        uint32_t i{0};

        for (int32_t index = 0; index < num_particles; ++index) {
            float x = positions[2 * index];
            float y = positions[2 * index + 1];
            if (x > minX && x < maxX && y > minY && y < maxY) {
                
                int32_t cellOccupantsX = x / cellSpacing;
                int32_t cellOccupantsY = y / cellSpacing;
                cellOccupants.addAtom(cellOccupantsX, cellOccupantsY, i);

            }
            ++i;
        }
    }

    void GetStaggeredGridPositions(float px, float py, int& x0, int& x1, int& y0, int& y1, int component) {

        // if you ever implement extrapolation, see if this bug goes away
        x0 = clamp(static_cast<int>(std::floor(px * invSpacing)), 1, numX - 2);
        x1 = std::min(x0 + 1, numX - 2); // <- see if you can make this 1
    
        y0 = clamp(static_cast<int>(std::floor(py * invSpacing)), 1, numY - 2);
        y1 = std::min(y0 + 1, numY - 1);

        // this stops the fluid touching the left wall from going down
        if (component == 0 && x0 == 1) { // ceiling
            x0 = x1;
        }
        if (component == 1 && x0 == 1) { // left wall
            x1 = x0;
        }
    }

    float curl(int i, int j) {
        int idx = i * n + j;
        const float denom = 1.f / (2.f * cellSpacing);
        const int leftType = cellType[idx - n] == 0;
        const int rightType = cellType[idx + n] == 0;
        const int topType = cellType[idx - 1] == 0;
        const int bottomType = cellType[idx + 1] == 0;
        if (!leftType || !rightType || !topType || !bottomType) {
            return 0.f;
        }
        return ((v[(i + 1) * n + j] * bottomType - v[(i - 1) * n + j] * topType) - (u[i * n + j + 1] * rightType - u[i * n + j - 1] * leftType)) * denom;
    }

    float calcVorticity(int i, int j) {
        float curl = this->curl(i, j);
        return curl * curl; // std::abs for just curl
    }

    // Returns world coordinates to the top left of the cell
    Vector2 gridCellToPos(int idx) {
        int i = idx % numY;
        int j = idx / numY;
        float x = j * cellSpacing;
        float y = i * cellSpacing;
        return Vector2{x, y};
    }

    int getNumParticles() {
        return this->num_particles;
    }

    bool getStop() {
        return this->stop;
    }

    void setStop(bool set) {
        this->stop = set;
    }

    bool getStep() {
        return this->step;
    }

    void setStep(bool set) {
        this->step = set;
    }

    void addToVorticityStrength(float add) {
        this->vorticityStrength += add;
    }

    float getVorticityStrength() {
        return this->vorticityStrength;
    }

    float getFlipRatio() {
        return this->flipRatio;
    }

    void addToFlipRatio(const float add) {
        this->flipRatio += add;
    }

    void setFlipRatio(float set) {
        this->flipRatio = set;
    }

    void addToGravityX(float add) {
        this->gravityX += add;
    }

    float getGravityX() {
        return this->gravityX;
    }

    void addToGravityY(float add) {
        this->gravityY += add;
    }

    float getGravityY() {
        return this->gravityY;
    }

    void setFireActive(bool active) {
        this->fireActive = active;
    }

    float getFireActive() {
        return this->fireActive;
    }

    int getNumX() {
        return this->numX;
    }

    int getNumY() {
        return this->numY;
    }
};