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

    FluidState &fas;

public:

    TransferGrid(FluidState &fluid_attributes): fas(fluid_attributes) {
        invSpacing =  1.f / fas.cellSpacing;
        halfSpacing = 0.5f * fas.cellSpacing;
        n = fas.numY;

        Cu.resize(2 * fas.num_particles);
        Cv.resize(2 * fas.num_particles);
        std::fill(begin(Cu), end(Cu), 0.f);
        std::fill(begin(Cv), end(Cv), 0.f);

        SOLID_CELL = fas.SOLID;
        FLUID_CELL = fas.FLUID;
        AIR_CELL = fas.AIR;

        volume = fas.cellSpacing * fas.cellSpacing;
        invVolume = 1.f / volume;
    }

    void TransferToGrid() {
        SetUpTransferGrid();

        P2G();

        MomentumToVelocity();

        NeumannBC();

        // store prev grid for FLIP
        std::copy(std::begin(fas.u), std::end(fas.u), std::begin(fas.prevU));
        std::copy(std::begin(fas.v), std::end(fas.v), std::begin(fas.prevV));
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
        std::copy(std::begin(fas.u), std::end(fas.u), std::begin(fas.prevU));
        std::copy(std::begin(fas.v), std::end(fas.v), std::begin(fas.prevV));
        std::fill(begin(fas.sumUGridWeights), end(fas.sumUGridWeights), 0.f);
        std::fill(begin(fas.sumVGridWeights), end(fas.sumVGridWeights), 0.f);
        std::fill(begin(fas.u), end(fas.u), 0.f);
        std::fill(begin(fas.v), end(fas.v), 0.f);
    }

    void P2G() {
        for (int component = 0; component < 2; ++component) {
            auto &gridVel = component == 0 ? fas.u : fas.v;
            auto &gridWeights = component == 0 ? fas.sumUGridWeights : fas.sumVGridWeights;
            auto &C = component == 0 ? Cu : Cv;

            for (int pIdx = 0; pIdx < fas.num_particles; ++pIdx) {
                const int32_t i = 2 * pIdx;

                float px = fas.positions[i];
                float py = fas.positions[i + 1];

                px -= (component == 1) * halfSpacing;
                py -= (component == 0) * halfSpacing;

                px = clamp(px, fas.cellSpacing, (fas.numX - 1) * fas.cellSpacing);
                py = clamp(py, fas.cellSpacing, (fas.numY - 1) * fas.cellSpacing);

                int x0, x1, y0, y1;
                x0 = x1 = y0 = y1 = 0;

                fas.GetStaggeredGridPositions(px, py, x0, x1, y0, y1, component);

                int32_t topLeftCell     = x0 * n + y0;
                int32_t topRightCell    = x1 * n + y0;
                int32_t bottomLeftCell  = x0 * n + y1;
                int32_t bottomRightCell = x1 * n + y1;

                float topLeftWeight     = fas.Weight(px, py, x0, y0);
                float topRightWeight    = fas.Weight(px, py, x1, y0);
                float bottomLeftWeight  = fas.Weight(px, py, x0, y1);
                float bottomRightWeight = fas.Weight(px, py, x1, y1);

                float dxLeft   = (x0 * fas.cellSpacing - px) * invSpacing;
                float dyTop    = (y0 * fas.cellSpacing - py) * invSpacing;
                float dxRight  = (x1 * fas.cellSpacing - px) * invSpacing;
                float dyBottom = (y1 * fas.cellSpacing - py) * invSpacing;

                float pv = fas.velocities[i + component];

                int matIdx = 2 * pIdx;

                float C1 = C[matIdx];
                float C2 = C[matIdx + 1];

                float affinepvBottomLeft = pv + 
                                        C1 * dxLeft + 
                                        C2 * dyBottom;

                float affinepvBottomRight = pv + 
                                        C1 * dxRight + 
                                        C2 * dyBottom;

                float affinepvTopLeft = pv + 
                                        C1 * dxLeft + 
                                        C2 * dyTop;

                float affinepvTopRight = pv + 
                                        C1 * dxRight + 
                                        C2 * dyTop;
            
                gridVel[topLeftCell]     += affinepvTopLeft     * topLeftWeight;
                gridVel[topRightCell]    += affinepvTopRight    * topRightWeight;
                gridVel[bottomLeftCell]  += affinepvBottomLeft  * bottomLeftWeight;
                gridVel[bottomRightCell] += affinepvBottomRight * bottomRightWeight;

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

    void MomentumToVelocity() {
        for (int i = 0; i < fas.gridSize; ++i) {
            float prevNode = fas.sumUGridWeights[i];
            if (prevNode > 0.f) {
                fas.u[i] /= prevNode;
            }
            prevNode = fas.sumVGridWeights[i];
            if (prevNode > 0.f) {
                fas.v[i] /= prevNode;
            }
        }
    }

    void NeumannBC() {
        for (int i = 0; i < fas.numX; ++i) {
            for (int j = 0; j < fas.numY; ++j) {
                int idx = i * n + j;
                bool solid = fas.cellType[idx] == SOLID_CELL;
                if (solid || (i > 0 && fas.cellType[idx - n] == SOLID_CELL)) {
                    fas.u[idx] = 0;
                }
                if (solid || (j > 0 && fas.cellType[idx - 1] == SOLID_CELL)) {
                    fas.v[idx] = 0;
                }
            }
        }
    }

    void G2P() {
        for (int component = 0; component < 2; ++component) {
            auto &grid = component == 0 ? fas.u : fas.v;
            auto &prevGrid = component == 0 ? fas.prevU : fas.prevV;
            auto &C = component == 0 ? Cu : Cv;

            for (int32_t pIdx = 0; pIdx < fas.num_particles; ++pIdx) {
                int i = 2 * pIdx;

                float px = fas.positions[i];
                float py = fas.positions[i + 1];

                px = clamp(px, fas.cellSpacing, (fas.numX - 1) * fas.cellSpacing);
                py = clamp(py, fas.cellSpacing, (fas.numY - 1) * fas.cellSpacing);

                px -= (component == 1) * halfSpacing;
                py -= (component == 0) * halfSpacing;

                int x0, x1, y0, y1;
                x0 = x1 = y0 = y1 = 0;

                fas.GetStaggeredGridPositions(px, py, x0, x1, y0, y1, component);

                int32_t topLeftCell     = x0 * n + y0;
                int32_t topRightCell    = x1 * n + y0;
                int32_t bottomLeftCell  = x0 * n + y1;
                int32_t bottomRightCell = x1 * n + y1;

                float topLeftWeight     = fas.Weight(px, py, x0, y0);
                float topRightWeight    = fas.Weight(px, py, x1, y0);
                float bottomLeftWeight  = fas.Weight(px, py, x0, y1);
                float bottomRightWeight = fas.Weight(px, py, x1, y1);

                float pv = fas.velocities[i + component];

                float gridVel_topLeft     = grid[topLeftCell];
                float gridVel_topRight    = grid[topRightCell];
                float gridVel_bottomRight = grid[bottomRightCell];
                float gridVel_bottomLeft  = grid[bottomLeftCell];
            
                int direction = component == 0 ? n : 1;
                const bool validTopLeftWeight     = fas.cellType[topLeftCell]                 != AIR_CELL || 
                                                    fas.cellType[topLeftCell - direction]     != AIR_CELL;
                const bool validTopRightWeight    = fas.cellType[topRightCell]                != AIR_CELL || 
                                                    fas.cellType[topRightCell - direction]    != AIR_CELL;
                const bool validBottomLeftWeight  = fas.cellType[bottomLeftCell]              != AIR_CELL || 
                                                    fas.cellType[bottomLeftCell - direction]  != AIR_CELL;
                const bool validBottomRightWeight = fas.cellType[bottomRightCell]             != AIR_CELL || 
                                                    fas.cellType[bottomRightCell - direction] != AIR_CELL;

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
                    fas.velocities[i + component] = (1.f - fas.flipRatio) * picV + fas.flipRatio * flipV;
                }

                const bool validTopLeftWeightX     = fas.cellType[topLeftCell]         != AIR_CELL || 
                                                     fas.cellType[topLeftCell - n]     != AIR_CELL;
                const bool validTopRightWeightX    = fas.cellType[topRightCell]        != AIR_CELL || 
                                                     fas.cellType[topRightCell - n]    != AIR_CELL;
                const bool validBottomLeftWeightX  = fas.cellType[bottomLeftCell]      != AIR_CELL || 
                                                     fas.cellType[bottomLeftCell - n]  != AIR_CELL;
                const bool validBottomRightWeightX = fas.cellType[bottomRightCell]      != AIR_CELL || 
                                                     fas.cellType[bottomRightCell - n] != AIR_CELL;

                const bool validTopLeftWeightY     = fas.cellType[topLeftCell]         != AIR_CELL || 
                                                     fas.cellType[topLeftCell - 1]     != AIR_CELL;
                const bool validTopRightWeightY    = fas.cellType[topRightCell]        != AIR_CELL || 
                                                     fas.cellType[topRightCell - 1]    != AIR_CELL;
                const bool validBottomLeftWeightY  = fas.cellType[bottomLeftCell]      != AIR_CELL || 
                                                     fas.cellType[bottomLeftCell - 1]  != AIR_CELL;
                const bool validBottomRightWeightY = fas.cellType[bottomRightCell]     != AIR_CELL || 
                                                     fas.cellType[bottomRightCell - 1] != AIR_CELL;

                float gridVelX_topLeft     = validTopLeftWeightX     * fas.u[topLeftCell];
                float gridVelX_topRight    = validTopRightWeightX    * fas.u[topRightCell];
                float gridVelX_bottomLeft  = validBottomLeftWeightX  * fas.u[bottomLeftCell];
                float gridVelX_bottomRight = validBottomRightWeightX * fas.u[bottomRightCell];

                float gridVelY_topLeft     = validTopLeftWeightY     * fas.v[topLeftCell];
                float gridVelY_topRight    = validTopRightWeightY    * fas.v[topRightCell];
                float gridVelY_bottomLeft  = validBottomLeftWeightY  * fas.v[bottomLeftCell];
                float gridVelY_bottomRight = validBottomRightWeightY * fas.v[bottomRightCell];

                float topLeftWeightDerivativeX, topRightWeightDerivativeX, bottomLeftWeightDerivativeX, bottomRightWeightDerivativeX;
                float topLeftWeightDerivativeY, topRightWeightDerivativeY, bottomLeftWeightDerivativeY, bottomRightWeightDerivativeY;

                fas.WeightGradientFD(px, py, x0, y0, topLeftWeightDerivativeX, topLeftWeightDerivativeY);
                fas.WeightGradientFD(px, py, x1, y0, topRightWeightDerivativeX, topRightWeightDerivativeY);
                fas.WeightGradientFD(px, py, x0, y1, bottomLeftWeightDerivativeX, bottomLeftWeightDerivativeY);
                fas.WeightGradientFD(px, py, x1, y1, bottomRightWeightDerivativeX, bottomRightWeightDerivativeY);

                C[i] = (topLeftWeightDerivativeX    * gridVelX_topLeft +
                       topRightWeightDerivativeX    * gridVelX_topRight +
                       bottomLeftWeightDerivativeX  * gridVelX_bottomLeft +
                       bottomRightWeightDerivativeX * gridVelX_bottomRight);

                C[i + 1] = (topLeftWeightDerivativeY    * gridVelY_topLeft +
                           topRightWeightDerivativeY    * gridVelY_topRight +
                           bottomLeftWeightDerivativeY  * gridVelY_bottomLeft +
                           bottomRightWeightDerivativeY * gridVelY_bottomRight);

                /*if (pIdx == 35) {
                    std::cout << "Top Left: ";

                    float gradXTopLeft;
                    float gradYTopLeft;

                    WeightGradientFD(px, py, x0, y0, gradXTopLeft, gradYTopLeft);

                    std::cout << topLeftWeightDerivativeX << ", " << gradXTopLeft << ", " << topLeftWeightDerivativeY << ", " << gradYTopLeft << "\n";


                    std::cout << "Top Right: ";

                    float gradXTopRight;
                    float gradYTopRight;

                    WeightGradientFD(px, py, x0, y0, gradXTopRight, gradYTopRight);

                    std::cout << topRightWeightDerivativeX << ", " << gradXTopRight << ", " << topRightWeightDerivativeY << ", " << gradYTopRight << "\n";


                    std::cout << "Bottom Left: ";

                    float gradXBottomLeft;
                    float gradYBottomLeft;

                    WeightGradientFD(px, py, x0, y0, gradXBottomLeft, gradYBottomLeft);

                    std::cout << bottomLeftWeightDerivativeX << ", " << gradXBottomLeft << ", " << bottomLeftWeightDerivativeY << ", " << gradYBottomLeft << "\n";

                    std::cout << "Bottom Right: ";

                    float gradXBottomRight;
                    float gradYBottomRight;

                    WeightGradientFD(px, py, x0, y0, gradXBottomRight, gradYBottomRight);

                    std::cout << bottomRightWeightDerivativeX << ", " << gradXBottomRight << ", " << bottomRightWeightDerivativeY << ", " << gradYBottomRight << "\n";
                }*/
            }
        }
    }
};