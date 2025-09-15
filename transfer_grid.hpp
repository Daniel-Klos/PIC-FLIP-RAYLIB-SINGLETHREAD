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
    float volume;
    float invVolume;

    int SOLID_CELL;
    int FLUID_CELL;
    int AIR_CELL;

    std::vector<float> Cu;
    std::vector<float> Cv;

    FluidState &fluid_attributes;

public:

    TransferGrid(FluidState &fas): fluid_attributes(fas) {
        invSpacing =  1.f / fluid_attributes.cellSpacing;
        halfSpacing = 0.5f * fluid_attributes.cellSpacing;
        n = fluid_attributes.numY;

        Cu.resize(2 * fluid_attributes.num_particles);
        Cv.resize(2 * fluid_attributes.num_particles);
        std::fill(begin(Cu), end(Cu), 0.f);
        std::fill(begin(Cv), end(Cv), 0.f);

        SOLID_CELL = fluid_attributes.SOLID;
        FLUID_CELL = fluid_attributes.FLUID;
        AIR_CELL = fluid_attributes.AIR;

        volume = fluid_attributes.cellSpacing * fluid_attributes.cellSpacing;
        invVolume = 1.f / volume;
    }

    void TransferToGrid() {
        SetUpTransferGrid();

        P2G();

        momentumToVelocity();

        NeumannBC();

        // store prev grid for FLIP
        std::copy(std::begin(fluid_attributes.u), std::end(fluid_attributes.u), std::begin(fluid_attributes.prevU));
        std::copy(std::begin(fluid_attributes.v), std::end(fluid_attributes.v), std::begin(fluid_attributes.prevV));
    }

    void TransferToParticles() {
        G2P();
    }

    void Resize(int newSize) {
        Cu.resize(2 * newSize);
        Cv.resize(2 * newSize);
    }

private:

    void SetUpTransferGrid() {
        std::copy(std::begin(fluid_attributes.u), std::end(fluid_attributes.u), std::begin(fluid_attributes.prevU));
        std::copy(std::begin(fluid_attributes.v), std::end(fluid_attributes.v), std::begin(fluid_attributes.prevV));
        std::fill(begin(fluid_attributes.sumUGridWeights), end(fluid_attributes.sumUGridWeights), 0.f);
        std::fill(begin(fluid_attributes.sumVGridWeights), end(fluid_attributes.sumVGridWeights), 0.f);
        std::fill(begin(fluid_attributes.u), end(fluid_attributes.u), 0.f);
        std::fill(begin(fluid_attributes.v), end(fluid_attributes.v), 0.f);
    }

    void P2G() {
        for (int component = 0; component < 2; ++component) {
            auto &gridVel = component == 0 ? fluid_attributes.u : fluid_attributes.v;
            auto &gridWeights = component == 0 ? fluid_attributes.sumUGridWeights : fluid_attributes.sumVGridWeights;
            auto &C = component == 0 ? Cu : Cv;

            for (int pIdx = 0; pIdx < fluid_attributes.num_particles; ++pIdx) {
                const int32_t i = 2 * pIdx;

                float px = fluid_attributes.positions[i];
                float py = fluid_attributes.positions[i + 1];

                px -= (component == 1) * halfSpacing;
                py -= (component == 0) * halfSpacing;

                px = clamp(px, fluid_attributes.cellSpacing, (fluid_attributes.numX - 1) * fluid_attributes.cellSpacing);
                py = clamp(py, fluid_attributes.cellSpacing, (fluid_attributes.numY - 1) * fluid_attributes.cellSpacing);

                int x0, x1, y0, y1;
                x0 = x1 = y0 = y1 = 0;

                fluid_attributes.GetStaggeredGridPositions(px, py, x0, x1, y0, y1, component);

                int32_t topLeftCell     = x0 * n + y0;
                int32_t topRightCell    = x1 * n + y0;
                int32_t bottomLeftCell  = x0 * n + y1;
                int32_t bottomRightCell = x1 * n + y1;

                float topLeftWeight     = fluid_attributes.Weight(px, py, x0, y0);
                float topRightWeight    = fluid_attributes.Weight(px, py, x1, y0);
                float bottomLeftWeight  = fluid_attributes.Weight(px, py, x0, y1);
                float bottomRightWeight = fluid_attributes.Weight(px, py, x1, y1);

                float dxLeft   = -(px - x0 * fluid_attributes.cellSpacing) * invSpacing;
                float dyTop    = -(py - y0 * fluid_attributes.cellSpacing) * invSpacing;
                float dxRight  = 1.f - abs(dxLeft);
                float dyBottom = 1.f - abs(dyTop);

                float pv = fluid_attributes.velocities[i + component];

                int matIdx = 2 * pIdx;

                float C1 = C[matIdx];
                float C2 = C[matIdx + 1];

                float affinepvxBottomLeft = pv + 
                                        C1 * dxLeft + 
                                        C2 * dyBottom;

                float affinepvxBottomRight = pv + 
                                        C1 * dxRight + 
                                        C2 * dyBottom;

                float affinepvxTopLeft = pv + 
                                        C1 * dxLeft + 
                                        C2 * dyTop;

                float affinepvTopRight = pv + 
                                        C1 * dxRight + 
                                        C2 * dyTop;
            
                gridVel[topLeftCell]     += affinepvxTopLeft     * topLeftWeight;
                gridVel[topRightCell]    += affinepvTopRight     * topRightWeight;
                gridVel[bottomLeftCell]  += affinepvxBottomLeft  * bottomLeftWeight;
                gridVel[bottomRightCell] += affinepvxBottomRight * bottomRightWeight;

                //grid[topLeftCell]     += pv * topLeftWeight;
                //grid[topRightCell]    += pv * topRightWeight;
                //grid[bottomLeftCell]  += pv * bottomLeftWeight;
                //grid[bottomRightCell] += pv * bottomRightWeight;

                gridWeights[topLeftCell]     += topLeftWeight;
                gridWeights[topRightCell]    += topRightWeight;
                gridWeights[bottomLeftCell]  += bottomLeftWeight;
                gridWeights[bottomRightCell] += bottomRightWeight;
            }
        }
    }

    void momentumToVelocity() {
        for (int i = 0; i < fluid_attributes.gridSize; ++i) {
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

    void NeumannBC() {
        for (int i = 0; i < fluid_attributes.numX - 1; ++i) {
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

    void G2P() {
        for (int component = 0; component < 2; ++component) {
            auto &grid = component == 0 ? fluid_attributes.u : fluid_attributes.v;
            auto &prevGrid = component == 0 ? fluid_attributes.prevU : fluid_attributes.prevV;
            auto &C = component == 0 ? Cu : Cv;

            for (int32_t pIdx = 0; pIdx < fluid_attributes.num_particles; ++pIdx) {
                int i = 2 * pIdx;

                float px = fluid_attributes.positions[i];
                float py = fluid_attributes.positions[i + 1];

                px = clamp(px, fluid_attributes.cellSpacing, (fluid_attributes.numX - 1) * fluid_attributes.cellSpacing);
                py = clamp(py, fluid_attributes.cellSpacing, (fluid_attributes.numY - 1) * fluid_attributes.cellSpacing);

                px -= (component == 1) * halfSpacing;
                py -= (component == 0) * halfSpacing;

                int x0, x1, y0, y1;
                x0 = x1 = y0 = y1 = 0;

                fluid_attributes.GetStaggeredGridPositions(px, py, x0, x1, y0, y1, component);

                int32_t topLeftCell     = x0 * n + y0;
                int32_t topRightCell    = x1 * n + y0;
                int32_t bottomLeftCell  = x0 * n + y1;
                int32_t bottomRightCell = x1 * n + y1;

                float topLeftWeight     = fluid_attributes.Weight(px, py, x0, y0);
                float topRightWeight    = fluid_attributes.Weight(px, py, x1, y0);
                float bottomLeftWeight  = fluid_attributes.Weight(px, py, x0, y1);
                float bottomRightWeight = fluid_attributes.Weight(px, py, x1, y1);

                float pv = fluid_attributes.velocities[i + component];

                float gridVel_topLeft     = grid[topLeftCell];
                float gridVel_topRight    = grid[topRightCell];
                float gridVel_bottomRight = grid[bottomRightCell];
                float gridVel_bottomLeft  = grid[bottomLeftCell];
            
                int direction = component == 0 ? n : 1;
                const bool validTopLeftWeight     = fluid_attributes.cellType[topLeftCell]                 != AIR_CELL || 
                                                    fluid_attributes.cellType[topLeftCell - direction]     != AIR_CELL;
                const bool validTopRightWeight    = fluid_attributes.cellType[topRightCell]                != AIR_CELL || 
                                                    fluid_attributes.cellType[topRightCell - direction]    != AIR_CELL;
                const bool validBottomLeftWeight  = fluid_attributes.cellType[bottomLeftCell]              != AIR_CELL || 
                                                    fluid_attributes.cellType[bottomLeftCell - direction]  != AIR_CELL;
                const bool validBottomRightWeight = fluid_attributes.cellType[bottomRightCell]             != AIR_CELL || 
                                                    fluid_attributes.cellType[bottomRightCell - direction] != AIR_CELL;

                const float div =  validTopLeftWeight     * topLeftWeight + 
                                   validTopRightWeight    * topRightWeight + 
                                   validBottomLeftWeight  * bottomLeftWeight + 
                                   validBottomRightWeight * bottomRightWeight;

                float picV;
                float corr;
                float flipV;

                if (div > 0.f) {
                    picV = (validTopLeftWeight     * topLeftWeight     * gridVel_topLeft + 
                            validTopRightWeight    * topRightWeight    * gridVel_topRight + 
                            validBottomLeftWeight  * bottomLeftWeight  * gridVel_bottomLeft + 
                            validBottomRightWeight * bottomRightWeight * gridVel_bottomRight
                            ) / div;

                    corr = (validTopLeftWeight     * topLeftWeight     * (gridVel_topLeft     - prevGrid[topLeftCell]) + 
                            validTopRightWeight    * topRightWeight    * (gridVel_topRight    - prevGrid[topRightCell]) + 
                            validBottomLeftWeight  * bottomLeftWeight  * (gridVel_bottomLeft  - prevGrid[bottomLeftCell]) + 
                            validBottomRightWeight * bottomRightWeight * (gridVel_bottomRight - prevGrid[bottomRightCell])
                            ) / div;

                    flipV = pv + corr;
                    fluid_attributes.velocities[i + component] = (1.f - fluid_attributes.flipRatio) * picV + fluid_attributes.flipRatio * flipV;
                }

                const bool validTopLeftWeightX     = fluid_attributes.cellType[topLeftCell]         != AIR_CELL || 
                                                     fluid_attributes.cellType[topLeftCell - n]     != AIR_CELL;
                const bool validTopRightWeightX    = fluid_attributes.cellType[topRightCell]        != AIR_CELL || 
                                                     fluid_attributes.cellType[topRightCell - n]    != AIR_CELL;
                const bool validBottomLeftWeightX  = fluid_attributes.cellType[bottomLeftCell]      != AIR_CELL || 
                                                     fluid_attributes.cellType[bottomLeftCell - n]  != AIR_CELL;
                const bool validBottomRightWeightX = fluid_attributes.cellType[bottomRightCell]      != AIR_CELL || 
                                                     fluid_attributes.cellType[bottomRightCell - n] != AIR_CELL;

                const bool validTopLeftWeightY     = fluid_attributes.cellType[topLeftCell]         != AIR_CELL || 
                                                     fluid_attributes.cellType[topLeftCell - 1]     != AIR_CELL;
                const bool validTopRightWeightY    = fluid_attributes.cellType[topRightCell]        != AIR_CELL || 
                                                     fluid_attributes.cellType[topRightCell - 1]    != AIR_CELL;
                const bool validBottomLeftWeightY  = fluid_attributes.cellType[bottomLeftCell]      != AIR_CELL || 
                                                     fluid_attributes.cellType[bottomLeftCell - 1]  != AIR_CELL;
                const bool validBottomRightWeightY = fluid_attributes.cellType[bottomRightCell]     != AIR_CELL || 
                                                     fluid_attributes.cellType[bottomRightCell - 1] != AIR_CELL;

                float gridVelX_topLeft     = validTopLeftWeightX     * fluid_attributes.u[topLeftCell];
                float gridVelX_topRight    = validTopRightWeightX    * fluid_attributes.u[topRightCell];
                float gridVelX_bottomLeft  = validBottomLeftWeightX  * fluid_attributes.u[bottomLeftCell];
                float gridVelX_bottomRight = validBottomRightWeightX * fluid_attributes.u[bottomRightCell];

                float gridVelY_topLeft     = validTopLeftWeightY     * fluid_attributes.v[topLeftCell];
                float gridVelY_topRight    = validTopRightWeightY    * fluid_attributes.v[topRightCell];
                float gridVelY_bottomLeft  = validBottomLeftWeightY  * fluid_attributes.v[bottomLeftCell];
                float gridVelY_bottomRight = validBottomRightWeightY * fluid_attributes.v[bottomRightCell];

                float topLeftWeightDerivativeX, topRightWeightDerivativeX, bottomLeftWeightDerivativeX, bottomRightWeightDerivativeX;
                float topLeftWeightDerivativeY, topRightWeightDerivativeY, bottomLeftWeightDerivativeY, bottomRightWeightDerivativeY;

                fluid_attributes.WeightGradient(px, py, x0, y0, topLeftWeightDerivativeX, topLeftWeightDerivativeY);
                fluid_attributes.WeightGradient(px, py, x1, y0, topRightWeightDerivativeX, topRightWeightDerivativeY);
                fluid_attributes.WeightGradient(px, py, x0, y1, bottomLeftWeightDerivativeX, bottomLeftWeightDerivativeY);
                fluid_attributes.WeightGradient(px, py, x1, y1, bottomRightWeightDerivativeX, bottomRightWeightDerivativeY);

                C[i] = (topLeftWeightDerivativeX    * gridVelX_topLeft +
                       topRightWeightDerivativeX    * gridVelX_topRight +
                       bottomLeftWeightDerivativeX  * gridVelX_bottomLeft +
                       bottomRightWeightDerivativeX * gridVelX_bottomRight);

                C[i + 1] = (topLeftWeightDerivativeY    * gridVelY_topLeft +
                           topRightWeightDerivativeY    * gridVelY_topRight +
                           bottomLeftWeightDerivativeY  * gridVelY_bottomLeft +
                           bottomRightWeightDerivativeY * gridVelY_bottomRight);
            }
        }
    }
};