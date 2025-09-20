#pragma once
#include "simulation_state.hpp"
#include <random>
#include <algorithm>

struct IDPSolver {
    FluidState &fas;

    int n;
    int gridSize;
    int numDensityIters = 40;

    float overRelaxation = 1.9f;

    std::vector<float> δx;
    std::vector<float> δy;

    std::mt19937 gen;

    std::normal_distribution<float> distribution;

    float restDensity = 1.f; // 0.8f

    IDPSolver(FluidState &fas): fas(fas), 
                                gen(std::random_device{}()),
                                distribution(fas.halfSpacing * 0.5, fas.halfSpacing)

    {                          
        n = fas.numY;
        gridSize = fas.gridSize;

        δx.resize(gridSize);
        δy.resize(gridSize);
    }

    void ProjectDensity() {
        ComputeGridDensity();

        // HandleDegenerateConditions();

        ApplyBoundaryConditions();

        ProjectDensityRedBlackSOR();

        ApplyPositionChanges();
    }

    void ComputeGridDensity() {
        UpdateCellDensities();

        CalculateRealRestDensity();

        RefineCellDensities();
    }

    void ComputeGridDensityOnly() {
        UpdateCellDensities();

        CalculateRealRestDensity();
    }

    void UpdateCellDensities() {

        std::fill(begin(fas.cellDensities), end(fas.cellDensities), 0.f);

        for (int i = 0; i < fas.num_particles; ++i) {
            float x = fas.positions[2 * i];
            float y = fas.positions[2 * i + 1];

            x -= fas.halfSpacing;
            y -= fas.halfSpacing;

            x = clamp(x, fas.cellSpacing, (fas.numX - 1) * fas.cellSpacing);
            y = clamp(y, fas.cellSpacing, (fas.numY - 1) * fas.cellSpacing);

            // can't use fas.GetStaggeredGridPositions() with this because cellDensities is not a staggered grid
            int x0 = std::max(1, std::min((int)(std::floor(x * fas.invSpacing)), fas.numX - 2));
            int x1 = std::min(x0 + 1, fas.numX - 1);

            int y0 = std::max(1, std::min((int)(std::floor(y * fas.invSpacing)), fas.numY - 2));
            int y1 = std::min(y0 + 1, fas.numY - 2);
           
            float dxLeft = (x - x0 * fas.cellSpacing) * fas.invSpacing;
            float dyTop = (y - y0 * fas.cellSpacing) * fas.invSpacing;
            float dxRight = 1.f - dxLeft;
            float dyBottom = 1.f - dyTop;

            int topLeft = x0 * n + y0;
            int topRight = x1 * n + y0;
            int bottomLeft = x0 * n + y1;
            int bottomRight = x1 * n + y1;

            // if you ever implement multiphase, know that you should interpolate masses and divide by cell volume when computing cell densities
            fas.cellDensities[topLeft]     += dxRight * dyBottom;
            
            fas.cellDensities[topRight]    += dxLeft * dyBottom;

            fas.cellDensities[bottomLeft]  += dxRight * dyTop;

            fas.cellDensities[bottomRight] += dxLeft * dyTop;
        }
    }

    void CalculateRealRestDensity() {
        if (fas.particleRestDensity == 0.f) {
            float sum = 0.f;

            for (int i = 0; i < fas.gridSize; ++i) {
                if (fas.cellType[i] == fas.FLUID) {
                    sum += fas.cellDensities[i];
                }
            }

            if (fas.num_fluid_cells > 0) {
                fas.particleRestDensity = sum / fas.num_fluid_cells;
            }
        }
    }

    // clamp cell densities next to air cells
    // apply ghost pressure
    // normalize by cell volume
    void RefineCellDensities() {

        for (int i = 0; i < fas.numX; ++i) {
            for (int j = 0; j < fas.numY; ++j) {
                int gIdx = i * fas.numY + j;
                if (fas.cellType[gIdx] != fas.FLUID) continue;

                int leftType   = fas.cellType[gIdx - fas.numY];
                int rightType  = fas.cellType[gIdx + fas.numY];
                int topType    = fas.cellType[gIdx - 1];
                int bottomType = fas.cellType[gIdx + 1];


                float ghostPressure = CalculateGhostPressures(gIdx);
                fas.cellDensities[gIdx] += ghostPressure;

                fas.cellDensities[gIdx] *= 1.f / fas.particleRestDensity;

                bool surrounding_air_cell = leftType   == fas.AIR ||
                                            rightType  == fas.AIR ||
                                            topType    == fas.AIR ||
                                            bottomType == fas.AIR;
                if (surrounding_air_cell) {
                    fas.cellDensities[gIdx] = std::max(fas.cellDensities[gIdx], restDensity);
                }
            }
        }
    }

    float CalculateGhostPressures(int gIdx) { 
        int solidLeft   = fas.cellType[gIdx - fas.numY]  == fas.SOLID;
        int solidRight  = fas.cellType[gIdx + fas.numY]  == fas.SOLID;
        int solidTop    = fas.cellType[gIdx - 1]         == fas.SOLID;
        int solidBottom = fas.cellType[gIdx + 1]         == fas.SOLID;

        int solidTopLeftCorner     = fas.cellType[gIdx - fas.numY - 1] == fas.SOLID;
        int solidTopRightCorner    = fas.cellType[gIdx + fas.numY - 1] == fas.SOLID;
        int solidBottomLeftCorner  = fas.cellType[gIdx - fas.numY + 1] == fas.SOLID;
        int solidBottomRightCorner = fas.cellType[gIdx + fas.numY + 1] == fas.SOLID;

        float edgeContribution   = 0.5 * fas.particleRestDensity;
        float cornerContribution = 0.25 * fas.particleRestDensity;

        float ghostPressure = edgeContribution   * (solidLeft + solidRight + solidTop + solidBottom) +
                              cornerContribution * (solidTopLeftCorner + solidTopRightCorner + solidBottomLeftCorner + solidBottomRightCorner);

        return ghostPressure;
    }

    // Redistribute cells that contain too many particles
    // If cellDensity > 1.5 * particleRestDensity, redistribute
    // Redistribution:
        // split cell into uniform subcells
        // place each particle randomly into a regoin close ot the center of a subcell to ensure the particles cover the whole volume of a cell
        // interpolate new particle velocities from the grid
    void HandleDegenerateConditions() {
        fas.FillCellOccupants();
        for (int i = 0; i < fas.num_fluid_cells; ++i) {
            int cellIdx = fas.fluid_cells[i];
            if (fas.cellDensities[cellIdx] > 1.5 * fas.particleRestDensity) { //1.5
                RedistributeParticles(cellIdx);
            }
        }
    }

    void RedistributeParticles(int cellIdx) {
        Vector2 cellPos = fas.gridCellToPos(cellIdx);
        auto cell = fas.cellOccupants.data[cellIdx];

        float particleCount = static_cast<float>(cell.objects_count);
        int particlesPerQuadrant = std::ceil(particleCount / 4.f);

        // Split cell into 4 quadrants and redistribute particles evenly to them
        float offsetsX[4] = {0.f,  fas.halfSpacing,        0,         fas.halfSpacing};
        float offsetsY[4] = {0.f,          0,       fas.halfSpacing,  fas.halfSpacing};
        for (int i = 0; i < 4; ++i) {
            float Qx = cellPos.x + offsetsX[i];
            float Qy = cellPos.y + offsetsY[i];

            for (int j = 0; j < particlesPerQuadrant; ++j) {
                int dataIdx = i * particlesPerQuadrant + j;
                if (dataIdx > particleCount - 1) {
                    break;
                }

                int pIdx = cell.objects[dataIdx];

                float offsetX = distribution(gen);
                float offsetY = distribution(gen);

                fas.positions[2 * pIdx] = Qx + offsetX;
                fas.positions[2 * pIdx + 1] = Qy + offsetY;
            }
        }
    }

    // similar boundary conditions as normal:
        // Dirichlet for air boundaries
        // Neumann for solid boundaries
            // This time though, we use δx * n = 0, not u * n = 0
    // Push-Out Boundary Condition
        // for every solid cell with particle(s) in it:
            // query the level set SDF for each particle to see where it should be pushed out
            // if there are multiple particles in the cell, use the SDF query of the particle with the deepest penetration on all particles in the cell
        // avoiding very drastic position corrections:
            // clamp displacement δx at the solid boundary to half of the cell width
    // Ghost pressures for solid obstacle cells ("ghost", because pressures usually only defined for fluid and air cells)
        // lots of derivations later...
        // we see that the push-out Neumann BC can be implemented by subtracting  δx_i * particleRestDensity/(∆x * ∆t) from the RHS

    // Mark cells that need to ignore Neumann boundary conditions
        // when particles are about to be separated from SDF in the next timestep (we only care about particles actually being separated in the next timestep, remember that obstacle separation is clamped), δx at the boundary between the separation cell and the solid needs to be the displacement of deepest particle
            // set this displacement at ALL MAC faces of the cell
    void ApplyBoundaryConditions() {
        fas.CollideSurfaces();
    }

    // In pressure projection we use the formula:
        // (∆t / ρ0) * ∇^2 p1 = ∇ · u*                      (9)
    // Here we want to solve:
        // (∆t / ρ0) * ∇^2 p2 = (1 / ∆t)(1 − ρ∗(t)/ρ0)      (10)
    // Note: large density errors can lead to drastic correction displacements in a single timestep, which can cause oscillatoins and popping artifacts. 
        // To avoid this, just clamp the right hand side of Eq (10).
        // clamp ρ* /ρ0 to the interval [0.5, 1.5]
    void passRedBlackSOR(bool red) {
        float invDt = 1.f / fas.frame_context.dt;
        for (int i = 1; i < fas.numX - 1; ++i) {
            for (int j = (i + red) % 2; j < fas.numY - 1; j += 2) {
                int idx = i * n + j;
                if (fas.cellType[idx] != fas.FLUID) continue;

                float leftType   = fas.cellType[idx - n] <= fas.AIR ? 1 : 0;
                float rightType  = fas.cellType[idx + n] <= fas.AIR ? 1 : 0;
                float topType    = fas.cellType[idx - 1] <= fas.AIR ? 1 : 0;
                float bottomType = fas.cellType[idx + 1] <= fas.AIR ? 1 : 0;

                float divideBy = leftType + rightType + topType + bottomType;
            
                if (divideBy == 0.f) continue;

                float rhoRatio = fas.cellDensities[idx] / restDensity;

                //rhoRatio = clamp(rhoRatio, 0.5f, 1.5f);

                float RHS = (1.f - rhoRatio) * invDt;

                float solution = RHS / divideBy;
                solution *= overRelaxation;

                δx[idx]     += leftType   * solution;
                δx[idx + n] -= rightType  * solution;
                δy[idx]     += topType    * solution;
                δy[idx + 1] -= bottomType * solution;
            }
        }
    }

    void ProjectDensityRedBlackSOR() {
        std::fill(begin(δx), end(δx), 0.f);
        std::fill(begin(δy), end(δy), 0.f);
        
        for (int iter = 0; iter < numDensityIters; ++iter) {
            passRedBlackSOR(0);
            passRedBlackSOR(1);
        }
    }

    // To get actual position changes on the grid:
        // plug p2 into the equation u(t + ∆t) = u* - ∆t * (1 / ρ0) * ∇p    (5)
        // So,
            // δx = -(∆t^2 / ρ0) * ∇p2        (11)
    void ApplyPositionChanges() {
        float dt = fas.frame_context.dt;
        float mul = (dt * dt) / restDensity;
        for (int component = 0; component < 2; ++component) {

            auto δdir = component == 0 ? δx : δy;

            for (int32_t pIdx = 0; pIdx < fas.num_particles; ++pIdx) {
                int i = 2 * pIdx;
                float px = fas.positions[i];
                float py = fas.positions[i + 1];

                px = clamp(px, fas.cellSpacing, (fas.numX - 1) * fas.cellSpacing);
                py = clamp(py, fas.cellSpacing, (fas.numY - 1) * fas.cellSpacing);

                px -= (component == 1) * fas.halfSpacing;
                py -= (component == 0) * fas.halfSpacing;

                int x0, x1, y0, y1;
                x0 = x1 = y0 = y1 = 0;

                fas.GetStaggeredGridPositions(px, py, x0, x1, y0, y1, component);

                int32_t topLeftCell     = x0 * n + y0;
                int32_t topRightCell    = x1 * n + y0;
                int32_t bottomLeftCell  = x0 * n + y1;
                int32_t bottomRightCell = x1 * n + y1;

                float topLeftDelta    = δdir[topLeftCell];
                float topRightDelta   = δdir[topRightCell];
                float bottomLeftDelta = δdir[bottomLeftCell];
                float bottomRighDelta = δdir[bottomRightCell];

                float topLeftWeight     = fas.Weight(px, py, x0, y0);
                float topRightWeight    = fas.Weight(px, py, x1, y0);
                float bottomLeftWeight  = fas.Weight(px, py, x0, y1);
                float bottomRightWeight = fas.Weight(px, py, x1, y1);

                int direction = component == 0 ? n : 1;
                const bool validTopLeftWeight     = fas.cellType[topLeftCell]                     != fas.AIR   || 
                                                    fas.cellType[topLeftCell - direction]         != fas.AIR;
                const bool validTopRightWeight    = fas.cellType[topRightCell]                    != fas.AIR   || 
                                                    fas.cellType[topRightCell - direction]        != fas.AIR;
                const bool validBottomLeftWeight  = fas.cellType[bottomLeftCell]                  != fas.AIR   || 
                                                    fas.cellType[bottomLeftCell - direction]      != fas.AIR;
                const bool validBottomRightWeight = fas.cellType[bottomRightCell]                 != fas.AIR   || 
                                                    fas.cellType[bottomRightCell - direction]     != fas.AIR;

                const float div = validTopLeftWeight     * topLeftWeight + 
                                  validTopRightWeight    * topRightWeight + 
                                  validBottomLeftWeight  * bottomLeftWeight + 
                                  validBottomRightWeight * bottomRightWeight;

                if (div > 0.f) {
                    float deltaGrad = (validTopLeftWeight     * topLeftWeight     * topLeftDelta    + 
                                       validTopRightWeight    * topRightWeight    * topRightDelta   +
                                       validBottomLeftWeight  * bottomLeftWeight  * bottomLeftDelta +
                                       validBottomRightWeight * bottomRightWeight * bottomRighDelta
                                       ) / div;
                    
                    fas.positions[i + component] += mul * deltaGrad;
                }
            }
        }
    }

    void AddToRestDensity(float add) {
        if (restDensity + add > 0.f) {
            restDensity += add;
        }
        else {
            restDensity = 0.1f;
        }
    }

    float getRestDensity() {
        return restDensity;
    }
};