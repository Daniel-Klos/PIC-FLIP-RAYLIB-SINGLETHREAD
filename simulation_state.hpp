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
    std::vector<int> fluid_cells;
    std::vector<Vector2i> obstaclePositions;
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

    float particleRestDensity;

    float vorticityStrength;
    float flipRatio; // move this to TransferGrid struct

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