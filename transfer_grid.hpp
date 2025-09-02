#pragma once
#include <vector>
#include <cstdint>

#include "simulation_state.hpp"
#include "utils.hpp"
#include "collision_grid.hpp"

class TransferGrid {

    float invSpacing;
    float halfSpacing;
    int n;

    int SOLID_CELL;
    int FLUID_CELL;
    int AIR_CELL;

    std::vector<float> Cu;
    std::vector<float> Cv;

    std::vector<float> dxLeft;  // x_p - x_cell for u and v grids respectively
    std::vector<float> dxRight;
    std::vector<float> dyBottom;
    std::vector<float> dyTop;
    
    std::vector<int32_t> topLeftCells;
    std::vector<int32_t> topRightCells;
    std::vector<int32_t> bottomRightCells;
    std::vector<int32_t> bottomLeftCells;

    std::vector<float> topLeftWeights;
    std::vector<float> topRightWeights;
    std::vector<float> bottomLeftWeights;
    std::vector<float> bottomRightWeights;

    std::vector<float> topLeftWeightDerivs;
    std::vector<float> topRightWeightDerivs;
    std::vector<float> bottomLeftWeightDerivs;
    std::vector<float> bottomRightWeightDerivs;

    FluidState &fluid_attributes;

    CollisionGrid cellOccupantsGrid;

public:

    TransferGrid(FluidState &fas): fluid_attributes(fas) {
        invSpacing =  1.f / fluid_attributes.cellSpacing;
        halfSpacing = 0.5f * fluid_attributes.cellSpacing;
        n = fluid_attributes.numY;

        topLeftCells.resize(2 * fluid_attributes.num_particles);
        topRightCells.resize(2 * fluid_attributes.num_particles);
        bottomRightCells.resize(2 * fluid_attributes.num_particles);
        bottomLeftCells.resize(2 * fluid_attributes.num_particles);
    
        topLeftWeights.resize(2 * fluid_attributes.num_particles);
        topRightWeights.resize(2 * fluid_attributes.num_particles);
        bottomLeftWeights.resize(2 * fluid_attributes.num_particles);
        bottomRightWeights.resize(2 * fluid_attributes.num_particles);

        topLeftWeightDerivs.resize(2 * fluid_attributes.num_particles);
        topRightWeightDerivs.resize(2 * fluid_attributes.num_particles);
        bottomLeftWeightDerivs.resize(2 * fluid_attributes.num_particles);
        bottomRightWeightDerivs.resize(2 * fluid_attributes.num_particles);

        dxLeft.resize(2 * fluid_attributes.num_particles);
        dxRight.resize(2 * fluid_attributes.num_particles);
        dyBottom.resize(2 * fluid_attributes.num_particles);
        dyTop.resize(2 * fluid_attributes.num_particles);

        Cu.resize(2 * fluid_attributes.num_particles);
        Cv.resize(2 * fluid_attributes.num_particles);
        std::fill(begin(Cu), end(Cu), 0.f);
        std::fill(begin(Cv), end(Cv), 0.f);

        SOLID_CELL = fluid_attributes.SOLID_CELL;
        FLUID_CELL = fluid_attributes.FLUID_CELL;
        AIR_CELL = fluid_attributes.AIR_CELL;

        cellOccupantsGrid = CollisionGrid(fluid_attributes.numX, fluid_attributes.numY);
    }

    void TransferToGrid() {
        cacheTransferNodesMulti();
        setUpTransferGrids();
        transferParticleVelocitiesToGridMulti();
    }

    void TransferToParticles() {
        transferUGridToParticles(0, fluid_attributes.num_particles);
        transferVGridToParticles(0, fluid_attributes.num_particles);
    }

    void updateCellDensitiesMulti() {
        std::fill(begin(fluid_attributes.cellDensities), end(fluid_attributes.cellDensities), 0.f);
        
        updateCellDensities();

        calculateRestDensity();
    }

    void Resize(int newSize) {
        topLeftCells.resize(2 * newSize);
        topRightCells.resize(2 * newSize);
        bottomRightCells.resize(2 * newSize);
        bottomLeftCells.resize(2 * newSize);
    
        topLeftWeights.resize(2 * newSize);
        topRightWeights.resize(2 * newSize);
        bottomLeftWeights.resize(2 * newSize);
        bottomRightWeights.resize(2 * newSize);

        topLeftWeightDerivs.resize(2 * newSize);
        topRightWeightDerivs.resize(2 * newSize);
        bottomLeftWeightDerivs.resize(2 * newSize);
        bottomRightWeightDerivs.resize(2 * newSize);

        dxLeft.resize(2 * newSize);
        dxRight.resize(2 * newSize);
        dyBottom.resize(2 * newSize);
        dyTop.resize(2 * newSize);

        Cu.resize(2 * newSize);
        Cv.resize(2 * newSize);
    }

private:

    void GetGridPositions(float px, float py, int& x0, int& x1, int& y0, int& y1, int component) {
        x0 = std::max(1, std::min(static_cast<int>(std::floor(px * invSpacing)), fluid_attributes.numX - 2));
        x1 = std::min(x0 + 1, fluid_attributes.numX - 2);
    
        y0 = std::max(0, std::min(static_cast<int>(std::floor(py * invSpacing)), fluid_attributes.numY - 2));
        y1 = std::min(y0 + 1, fluid_attributes.numY - 1);

        // this stops the fluid touching the left wall from going down
        if (component == 0 && x0 == 1) { // ceiling
            x0 = x1;
        }
        if (component == 1 && x0 == 1) { // left wall
            x1 = x0;
        }
    }

    float Weight(float px, float py, float gx, float gy) {
        gx *= fluid_attributes.cellSpacing;
        gy *= fluid_attributes.cellSpacing;

        float dx = 1.f - std::abs((px - gx)) * invSpacing;
        float dy = 1.f - std::abs((py - gy)) * invSpacing;

        return dx * dy;
    }

    float WeightDerivativeFD(float px, float py, float gx, float gy, int component) {
        float eps = fluid_attributes.cellSpacing * 0.001f;

        if (component == 0) { // v
            return (Weight(px, py + eps, gx, gy) - Weight(px, py - eps, gx, gy)) / (2 * eps);
        }
        else { // u
            return (Weight(px + eps, py, gx, gy) - Weight(px - eps, py, gx, gy)) / (2 * eps);
        }
    }

    void WeightDerivative(float px, float py, float gx, float gy, float &gradx, float &grady) {
        gx *= fluid_attributes.cellSpacing;
        gy *= fluid_attributes.cellSpacing;

        float dx = (px - gx) * invSpacing;
        float dy = (py - gy) * invSpacing;

        gradx = -sign(dx) * (1.f - abs(dy));
        grady = -sign(dy) * (1.f - abs(dx));
    }

    // add RKstep as a parameter
    void cacheGridInfoAt(int pIdx, int component) {
        const float dx = (component == 1) * halfSpacing;
        const float dy = (component == 0) * halfSpacing;

        float px = fluid_attributes.positions[2 * pIdx];
        float py = fluid_attributes.positions[2 * pIdx + 1];

        px = clamp(px, fluid_attributes.cellSpacing, (fluid_attributes.numX - 1) * fluid_attributes.cellSpacing);
        py = clamp(py, fluid_attributes.cellSpacing, (fluid_attributes.numY - 1) * fluid_attributes.cellSpacing);

        px -= dx;
        py -= dy;

        int x0, x1, y0, y1;
        x0 = x1 = y0 = y1 = 0;

        GetGridPositions(px, py, x0, x1, y0, y1, component);

        int gIdx = 2 * pIdx + component;

        dxLeft[gIdx]   = (px - x0 * fluid_attributes.cellSpacing) * invSpacing;
        dyTop[gIdx]    = (py - y0 * fluid_attributes.cellSpacing) * invSpacing;
        dxRight[gIdx]  = 1.f - dxLeft[gIdx];
        dyBottom[gIdx] = 1.f - dyTop[gIdx];
    
        topLeftWeights[gIdx]     = Weight(px, py, x0, y0); // * fluid_attributes.masses[pIdx];
        topRightWeights[gIdx]    = Weight(px, py, x1, y0); // * fluid_attributes.masses[pIdx];
        bottomLeftWeights[gIdx]  = Weight(px, py, x0, y1); // * fluid_attributes.masses[pIdx];
        bottomRightWeights[gIdx] = Weight(px, py, x1, y1); // * fluid_attributes.masses[pIdx];

        /*topLeftWeightDerivs[gIdx]     = WeightDerivative(px, py, x0, y0);
        topRightWeightDerivs[gIdx]    = WeightDerivative(px, py, x1, y0);
        bottomLeftWeightDerivs[gIdx]  = WeightDerivative(px, py, x0, y1);
        bottomRightWeightDerivs[gIdx] = WeightDerivative(px, py, x1, y1);*/

        topLeftCells[gIdx]     = x0 * n + y0;
        topRightCells[gIdx]    = x1 * n + y0;
        bottomLeftCells[gIdx]  = x0 * n + y1;
        bottomRightCells[gIdx] = x1 * n + y1;
    }

    // make topLeftCell a vector of arrays of size whatever RK you need
    void cacheTransferNodes(int32_t start, int32_t end, int32_t component) {
        for (int32_t i = start; i < end; ++i) {
            cacheGridInfoAt(i, component);
        }
    }

    void cacheTransferNodesMulti() {
        cacheTransferNodes(0, fluid_attributes.num_particles, 0);
        cacheTransferNodes(0, fluid_attributes.num_particles, 1);
    }

    void setUpTransferGrids() {
        std::copy(std::begin(fluid_attributes.u), std::end(fluid_attributes.u), std::begin(fluid_attributes.prevU));
        std::copy(std::begin(fluid_attributes.v), std::end(fluid_attributes.v), std::begin(fluid_attributes.prevV));
        std::fill(begin(fluid_attributes.sumUGridWeights), end(fluid_attributes.sumUGridWeights), 0.f);
        std::fill(begin(fluid_attributes.sumVGridWeights), end(fluid_attributes.sumVGridWeights), 0.f);
        std::fill(begin(fluid_attributes.u), end(fluid_attributes.u), 0.f);
        std::fill(begin(fluid_attributes.v), end(fluid_attributes.v), 0.f);

        // initialize every inside cell to air
        for (int i = 0; i < fluid_attributes.numX; ++i) {
            for (int j = 0; j < fluid_attributes.numY; ++j) {
                if (fluid_attributes.cellType[i * n + j] != SOLID_CELL) {
                    fluid_attributes.cellType[i * n + j] = AIR_CELL;
                }
            }
        }

        // initialize all cells that particles are in to fluid and store fluid cells
        fluid_attributes.num_fluid_cells = 0;
        for (int i = 0; i < fluid_attributes.num_particles; ++i) {
            float x = fluid_attributes.positions[2 * i];
            float y = fluid_attributes.positions[2 * i + 1];

            int xi = clamp(std::floor(x * invSpacing), 0, fluid_attributes.numX - 1);
            int yi = clamp(std::floor(y * invSpacing), 0, fluid_attributes.numY - 1);

            int cellNr = xi * n + yi;
            if (fluid_attributes.cellType[cellNr] == AIR_CELL) {
                fluid_attributes.cellType[cellNr] = FLUID_CELL;
                fluid_attributes.fluid_cells[fluid_attributes.num_fluid_cells] = cellNr;
                fluid_attributes.num_fluid_cells++;
            }
        }

        // add particles to occupantsGrid for multithreading
        cellOccupantsGrid.clear();

        const float minX = fluid_attributes.cellSpacing;
        const float maxX = fluid_attributes.frame_context.WIDTH - fluid_attributes.cellSpacing;
        const float minY = fluid_attributes.cellSpacing;
        const float maxY = fluid_attributes.frame_context.HEIGHT - fluid_attributes.cellSpacing;

        uint32_t i{0};

        for (int32_t index = 0; index < fluid_attributes.num_particles; ++index) {
            float x = fluid_attributes.positions[2 * index];
            float y = fluid_attributes.positions[2 * index + 1];
            if (x > minX && x < maxX && y > minY && y < maxY) {
                
                int32_t cellOccupantsX = x / fluid_attributes.cellSpacing;
                int32_t cellOccupantsY = y / fluid_attributes.cellSpacing;
                cellOccupantsGrid.addAtom(cellOccupantsX, cellOccupantsY, i);

            }
            ++i;
        }
    }

    void transferToUGridCells(int idx) {
        const auto cell = cellOccupantsGrid.data[idx];
    
        for (uint32_t id{0}; id < cell.objects_count; ++id) {
            const uint32_t pIdx = cell.objects[id];

            const int32_t ui = 2 * pIdx;

            int32_t topLeftCell_u     = topLeftCells[ui];
            int32_t topRightCell_u    = topRightCells[ui];
            int32_t bottomRightCell_u = bottomRightCells[ui];
            int32_t bottomLeftCell_u  = bottomLeftCells[ui];

            float topLeftWeight_u     = topLeftWeights[ui];
            float topRightWeight_u    = topRightWeights[ui];
            float bottomLeftWeight_u  = bottomLeftWeights[ui];
            float bottomRightWeight_u = bottomRightWeights[ui];

            float pvx = fluid_attributes.velocities[ui];

            float leftDist   = dxLeft[ui];
            float rightDist  = dxRight[ui];
            float bottomDist = dyBottom[ui];
            float topDist    = dyTop[ui];

            int matIdx = 2 * pIdx;

            float C11 = Cu[matIdx];
            float C21 = Cu[matIdx + 1];

            float affinepvxBottomLeft = pvx + 
                                    C11 * leftDist + 
                                    C21 * bottomDist;

            float affinepvxBottomRight = pvx + 
                                    C11 * rightDist + 
                                    C21 * bottomDist;

            float affinepvxTopLeft = pvx + 
                                    C11 * leftDist + 
                                    C21 * topDist;

            float affinepvTopRight = pvx + 
                                    C11 * rightDist + 
                                    C21 * topDist;
           
            /*fluid_attributes.u[topLeftCell_u]     += affinepvxTopLeft * topLeftWeight_u;  
            fluid_attributes.u[topRightCell_u]    += affinepvTopRight * topRightWeight_u;
            fluid_attributes.u[bottomLeftCell_u]  += affinepvxBottomLeft * bottomLeftWeight_u;
            fluid_attributes.u[bottomRightCell_u] += affinepvxBottomRight * bottomRightWeight_u;*/

            fluid_attributes.u[topLeftCell_u] += pvx * topLeftWeight_u;  
            fluid_attributes.u[topRightCell_u] += pvx * topRightWeight_u;
            fluid_attributes.u[bottomLeftCell_u] += pvx * bottomLeftWeight_u;
            fluid_attributes.u[bottomRightCell_u] += pvx * bottomRightWeight_u;

            fluid_attributes.sumUGridWeights[topLeftCell_u]     += topLeftWeight_u;
            fluid_attributes.sumUGridWeights[topRightCell_u]    += topRightWeight_u;
            fluid_attributes.sumUGridWeights[bottomLeftCell_u]  += bottomLeftWeight_u;
            fluid_attributes.sumUGridWeights[bottomRightCell_u] += bottomRightWeight_u;
        }
    }

    void transferToUGrid(int start, int end) {
        for (int i = start; i < end; ++i) {
            transferToUGridCells(fluid_attributes.fluid_cells[i]);
        }
    }

    void transferToVGridCells(int idx) {
        const auto cell = cellOccupantsGrid.data[idx];
    
        for (uint32_t id{0}; id < cell.objects_count; ++id) {
            const uint32_t pIdx = cell.objects[id];

            const int32_t vi = 2 * pIdx + 1;

            const int32_t topLeftCell_v     = topLeftCells[vi];
            const int32_t topRightCell_v    = topRightCells[vi];
            const int32_t bottomRightCell_v = bottomRightCells[vi];
            const int32_t bottomLeftCell_v  = bottomLeftCells[vi];

            const float topLeftWeight_v     = topLeftWeights[vi];
            const float topRightWeight_v    = topRightWeights[vi];
            const float bottomLeftWeight_v  = bottomLeftWeights[vi];
            const float bottomRightWeight_v = bottomRightWeights[vi];

            const float pvy = fluid_attributes.velocities[vi];

            float leftDist   = dxLeft[vi];
            float rightDist  = dxRight[vi];
            float bottomDist = dyBottom[vi];
            float topDist    = dyTop[vi];

            int matIdx = 2 * pIdx;

            float C12 = Cv[matIdx];
            float C22 = Cv[matIdx + 1];

            float affinepvyBottomLeft = pvy + 
                                    C12 * leftDist + 
                                    C22 * bottomDist;

            float affinepvyBottomRight = pvy + 
                                    C12 * rightDist + 
                                    C22 * bottomDist;

            float affinepvyTopLeft = pvy + 
                                    C12 * leftDist + 
                                    C22 * topDist;

            float affinepvyTopRight = pvy + 
                                    C12 * rightDist + 
                                    C22 * topDist;
    
            /*fluid_attributes.v[topLeftCell_v]     += affinepvyTopLeft * topLeftWeight_v;  
            fluid_attributes.v[topRightCell_v]    += affinepvyTopRight * topRightWeight_v;
            fluid_attributes.v[bottomLeftCell_v]  += affinepvyBottomLeft * bottomLeftWeight_v;
            fluid_attributes.v[bottomRightCell_v] += affinepvyBottomRight * bottomRightWeight_v;*/

            fluid_attributes.v[topLeftCell_v] += pvy * topLeftWeight_v;  
            fluid_attributes.v[topRightCell_v] += pvy * topRightWeight_v;
            fluid_attributes.v[bottomLeftCell_v] += pvy * bottomLeftWeight_v;
            fluid_attributes.v[bottomRightCell_v] += pvy * bottomRightWeight_v;

            fluid_attributes.sumVGridWeights[topLeftCell_v]     += topLeftWeight_v;
            fluid_attributes.sumVGridWeights[topRightCell_v]    += topRightWeight_v;
            fluid_attributes.sumVGridWeights[bottomLeftCell_v]  += bottomLeftWeight_v;
            fluid_attributes.sumVGridWeights[bottomRightCell_v] += bottomRightWeight_v;
        }
    }

    void transferToVGrid(int start, int end) {
        for (int i = start; i < end; ++i) {
            transferToVGridCells(fluid_attributes.fluid_cells[i]);
        }
    }

    void momentumToVelocity(int start, int end) {
        for (int i = start; i < end; ++i) {
            float prevNode = fluid_attributes.sumUGridWeights[i];
            if (prevNode > 0.f) {
                fluid_attributes.u[i] /= prevNode;
            }
            prevNode = fluid_attributes.sumVGridWeights[i];
            if (prevNode > 0.f) {
                fluid_attributes.v[i] /= prevNode;
            }
        }
    }

    void momentumToVelocityMulti() {
        momentumToVelocity(0, fluid_attributes.gridSize);
    }

    void enforceNoSlip(int start, int end) {
        for (int i = start; i < end; ++i) {
            for (int j = 0; j < fluid_attributes.numY; ++j) {
                int idx = i * n + j;
                bool solid = fluid_attributes.cellType[idx] == SOLID_CELL;
                if (solid || (i > 0 && fluid_attributes.cellType[idx - n] == SOLID_CELL)) {
                    fluid_attributes.u[idx] = 0;
                }
                if (solid || (j > 0 && fluid_attributes.cellType[idx - 1] == SOLID_CELL)) {
                    fluid_attributes.v[idx] = 0;
                }
            }
        }
    }

    void enforceNoSlipMulti() {
        enforceNoSlip(0, fluid_attributes.numX - 1);
    }

    void transferParticleVelocitiesToGridMulti() {
        transferToUGrid(0, fluid_attributes.num_fluid_cells);

        transferToVGrid(0, fluid_attributes.num_fluid_cells);

        momentumToVelocityMulti();

        enforceNoSlipMulti();

        // store prev grid for FLIP
        std::copy(std::begin(fluid_attributes.u), std::end(fluid_attributes.u), std::begin(fluid_attributes.prevU));
        std::copy(std::begin(fluid_attributes.v), std::end(fluid_attributes.v), std::begin(fluid_attributes.prevV));
    }

    void transferUGridToParticles(int32_t startIndex, int32_t endIndex) {
        for (int32_t i = startIndex; i < endIndex; ++i) {
            int ui = 2 * i;
            float pvx = fluid_attributes.velocities[ui];

            int topLeftCell     = topLeftCells[ui];
            int topRightCell    = topRightCells[ui];
            int bottomRightCell = bottomRightCells[ui];
            int bottomLeftCell  = bottomLeftCells[ui];

            float topLeftWeight     = topLeftWeights[ui];
            float topRightWeight    = topRightWeights[ui];
            float bottomLeftWeight  = bottomLeftWeights[ui];
            float bottomRightWeight = bottomRightWeights[ui];

            float gridVelX_topLeft     = fluid_attributes.u[topLeftCell];
            float gridVelX_topRight    = fluid_attributes.u[topRightCell];
            float gridVelX_bottomRight = fluid_attributes.u[bottomRightCell];
            float gridVelX_bottomLeft  = fluid_attributes.u[bottomLeftCell];
           
            const float validTopLeftWeight     = fluid_attributes.cellType[topLeftCell]         != AIR_CELL || 
                                                 fluid_attributes.cellType[topLeftCell - n]     != AIR_CELL;
            const float validTopRightWeight    = fluid_attributes.cellType[topRightCell]        != AIR_CELL || 
                                                 fluid_attributes.cellType[topRightCell - n]    != AIR_CELL;
            const float validBottomLeftWeight  = fluid_attributes.cellType[bottomLeftCell]      != AIR_CELL || 
                                                 fluid_attributes.cellType[bottomLeftCell - n]  != AIR_CELL;
            const float validBottomRightWeight = fluid_attributes.cellType[bottomRightCell]     != AIR_CELL || 
                                                 fluid_attributes.cellType[bottomRightCell - n] != AIR_CELL;

            const float divX = validTopLeftWeight * topLeftWeight + 
                               validTopRightWeight * topRightWeight + 
                               validBottomLeftWeight * bottomLeftWeight + 
                               validBottomRightWeight * bottomRightWeight;

            float picV;
            float corr;
            float flipV;

            if (divX > 0.f) {
                picV = (validTopLeftWeight     * topLeftWeight     * gridVelX_topLeft + 
                        validTopRightWeight    * topRightWeight    * gridVelX_topRight + 
                        validBottomLeftWeight  * bottomLeftWeight  * gridVelX_bottomLeft + 
                        validBottomRightWeight * bottomRightWeight * gridVelX_bottomRight
                        ) / divX;

                corr = (validTopLeftWeight     * topLeftWeight     * (gridVelX_topLeft - fluid_attributes.prevU[topLeftCell]) + 
                        validTopRightWeight    * topRightWeight    * (gridVelX_topRight - fluid_attributes.prevU[topRightCell]) + 
                        validBottomLeftWeight  * bottomLeftWeight  * (gridVelX_bottomLeft - fluid_attributes.prevU[bottomLeftCell]) + 
                        validBottomRightWeight * bottomRightWeight * (gridVelX_bottomRight - fluid_attributes.prevU[bottomRightCell])
                        ) / divX;

                flipV = pvx + corr;
                fluid_attributes.velocities[ui] = (1.f - fluid_attributes.flipRatio) * picV + fluid_attributes.flipRatio * flipV;
            }

            /*float fy = dyTop[i].vu[1];
            float fx = dxLeft[i].vu[1];

            int gIdx = 2 * i;
            Cu[gIdx] = ((fx - 1.f) * gridVelX_topLeft +
                       (1.f - fx)  * gridVelX_bottomLeft  +
                       (-fx)       * gridVelX_topRight  +
                       (fx)        * gridVelX_bottomRight) * invSpacing;
            
            Cu[gIdx + 1] = ((fy - 1.f) * gridVelX_topLeft +
                           (1.f - fy)  * gridVelX_bottomLeft  +
                           (-fy)       * gridVelX_topRight  +
                           (fy)        * gridVelX_bottomRight) * invSpacing;*/
        }
    }

    void transferVGridToParticles(int32_t startIndex, int32_t endIndex) {
        for (int32_t i = startIndex; i < endIndex; ++i) {
            int32_t vi = 2 * i + 1;
            float pvy = fluid_attributes.velocities[vi];

            int32_t topLeftCell     = topLeftCells[vi];
            int32_t topRightCell    = topRightCells[vi];
            int32_t bottomRightCell = bottomRightCells[vi];
            int32_t bottomLeftCell  = bottomLeftCells[vi];

            float topLeftWeight     = topLeftWeights[vi];
            float topRightWeight    = topRightWeights[vi];
            float bottomLeftWeight  = bottomLeftWeights[vi];
            float bottomRightWeight = bottomRightWeights[vi];

            float gridVelY_topLeft     = fluid_attributes.v[topLeftCell];
            float gridVelY_topRight    = fluid_attributes.v[topRightCell];
            float gridVelY_bottomRight = fluid_attributes.v[bottomRightCell];
            float gridVelY_bottomLeft  = fluid_attributes.v[bottomLeftCell];
           
            const float validTopLeftWeight     = fluid_attributes.cellType[topLeftCell]         != AIR_CELL || 
                                                 fluid_attributes.cellType[topLeftCell - 1]     != AIR_CELL;
            const float validTopRightWeight    = fluid_attributes.cellType[topRightCell]        != AIR_CELL || 
                                                 fluid_attributes.cellType[topRightCell - 1]    != AIR_CELL;
            const float validBottomLeftWeight  = fluid_attributes.cellType[bottomLeftCell]      != AIR_CELL || 
                                                 fluid_attributes.cellType[bottomLeftCell - 1]  != AIR_CELL;
            const float validBottomRightWeight = fluid_attributes.cellType[bottomRightCell]     != AIR_CELL || 
                                                 fluid_attributes.cellType[bottomRightCell - 1] != AIR_CELL;

            const float divY = validTopLeftWeight * topLeftWeight + 
                               validTopRightWeight * topRightWeight + 
                               validBottomLeftWeight * bottomLeftWeight + 
                               validBottomRightWeight * bottomRightWeight;

            float picV;
            float corr;
            float flipV;

            if (divY > 0.f) {
                picV = (validTopLeftWeight     * topLeftWeight     * gridVelY_topLeft + 
                        validTopRightWeight    * topRightWeight    * gridVelY_topRight + 
                        validBottomLeftWeight  * bottomLeftWeight  * gridVelY_bottomLeft + 
                        validBottomRightWeight * bottomRightWeight * gridVelY_bottomRight
                        ) / divY;

                corr = (validTopLeftWeight     * topLeftWeight     * (gridVelY_topLeft - fluid_attributes.prevV[topLeftCell]) + 
                        validTopRightWeight    * topRightWeight    * (gridVelY_topRight - fluid_attributes.prevV[topRightCell]) + 
                        validBottomLeftWeight  * bottomLeftWeight  * (gridVelY_bottomLeft - fluid_attributes.prevV[bottomLeftCell]) + 
                        validBottomRightWeight * bottomRightWeight * (gridVelY_bottomRight - fluid_attributes.prevV[bottomRightCell])
                        ) / divY;

                flipV = pvy + corr;
                fluid_attributes.velocities[vi] = (1.f - fluid_attributes.flipRatio) * picV + fluid_attributes.flipRatio * flipV;
            }

            /*float fy = dyTop[i].vu[0];
            float fx = dxLeft[i].vu[0];

            int gIdx = 2 * i;
            Cv[gIdx] = ((fx - 1.f) * gridVelY_topLeft +
                       (1.f - fx)  * gridVelY_bottomLeft  +
                       (-fx)       * gridVelY_topRight  +
                       (fx)        * gridVelY_bottomRight) * invSpacing;
            
            Cv[gIdx + 1] = ((fy - 1.f) * gridVelY_topLeft +
                           (1.f - fy)  * gridVelY_bottomLeft  +
                           (-fy)       * gridVelY_topRight  +
                           (fy)        * gridVelY_bottomRight) * invSpacing;*/
        }
    }

    void calculateRestDensity() {
        if (fluid_attributes.particleRestDensity == 0.f) {
            float sum = 0.f;
            int numFluidCells = 0;

            for (int i = 0; i < fluid_attributes.gridSize; ++i) {
                if (fluid_attributes.cellType[i] == FLUID_CELL) {
                    sum += fluid_attributes.cellDensities[i];
                    numFluidCells++;
                }
            }

            if (numFluidCells > 0) {
                fluid_attributes.particleRestDensity = sum / numFluidCells;
            }
        }
    }

    void updateCellDensities() {
        for (int i = 0; i < fluid_attributes.num_particles; ++i) {
            float x = fluid_attributes.positions[2 * i];
            float y = fluid_attributes.positions[2 * i + 1];

            x = clamp(x, fluid_attributes.cellSpacing, (fluid_attributes.numX - 1) * fluid_attributes.cellSpacing);
            y = clamp(y, fluid_attributes.cellSpacing, (fluid_attributes.numY - 1) * fluid_attributes.cellSpacing);

            int x0 = std::max(1, std::min((int)(std::floor((x - halfSpacing) * invSpacing)), fluid_attributes.numX - 2));
            int x1 = std::min(x0 + 1, fluid_attributes.numX - 1);

            int y0 = std::max(1, std::min((int)(std::floor((y - halfSpacing) * invSpacing)), fluid_attributes.numY - 2));
            int y1 = std::min(y0 + 1, fluid_attributes.numY - 2);
           
            float dxLeft = ((x - halfSpacing) - x0 * fluid_attributes.cellSpacing) * invSpacing;
            float dyTop = ((y - halfSpacing) - y0 * fluid_attributes.cellSpacing) * invSpacing;
            float dxRight = 1.f - dxLeft;
            float dyBottom = 1.f - dyTop;

            int topLeft = x0 * n + y0;
            int topRight = x1 * n + y0;
            int bottomLeft = x0 * n + y1;
            int bottomRight = x1 * n + y1;

            fluid_attributes.cellDensities[topLeft]     += dxRight * dyBottom;
            
            fluid_attributes.cellDensities[topRight]    += dxLeft * dyBottom;

            fluid_attributes.cellDensities[bottomLeft]  += dxRight * dyTop;

            fluid_attributes.cellDensities[bottomRight] += dxLeft * dyTop;
        }
    }

};