#pragma once
#include <functional>
#include "simulation_state.hpp"

struct PressureSolver {
    FluidState &fluid_attributes;
    int n;
    int gridSize;
    int numPressureIters = 40;
    int FLUID = fluid_attributes.FLUID;
    int AIR = fluid_attributes.AIR;
    int SOLID = fluid_attributes.SOLID;

    float k = 10.f;
    float overRelaxation = 1.9f;

    /*std::vector<double> residual;
    std::vector<double> Adiag;
    std::vector<double> precon;
    std::vector<double> direction;
    std::vector<double> dotProducts;
    std::vector<double> pressure;

    std::vector<uint8_t> neighbors;
    static constexpr uint8_t LEFT = 64;
    static constexpr uint8_t RIGHT = 32;
    static constexpr uint8_t BOTTOM = 16;
    static constexpr uint8_t TOP = 8;
    static constexpr uint8_t CENTER = 7;
    static constexpr std::array<uint8_t, 4> directions{LEFT, RIGHT, BOTTOM, TOP};*/

    std::vector<float> coarseU;
    int max_multigrid_coarseness = 4;

    PressureSolver(FluidState &fas): fluid_attributes(fas) {
        n = fluid_attributes.numY;

        gridSize = fluid_attributes.gridSize;
        /*this->Adiag.resize(gridSize);
        this->neighbors.resize(gridSize);
        this->precon.resize(gridSize);
        this->direction.resize(gridSize);
        this->residual.resize(gridSize);
        this->pressure.resize(gridSize);*/
    }


    // What I'm using currently
    // --------------------------------------------------------------------------------------------------------------------------------------------

    void passRedBlackSOR(bool red) {
        for (int i = 1; i < fluid_attributes.numX - 1; ++i) {
            for (int j = (i + red) % 2; j < fluid_attributes.numY - 1; j += 2) {
                int idx = i * n + j;
                if (fluid_attributes.cellType[idx] != FLUID) continue;

                float leftType = fluid_attributes.cellType[idx - n] <= AIR ? 1 : 0;
                float rightType = fluid_attributes.cellType[idx + n] <= AIR ? 1 : 0;
                float topType = fluid_attributes.cellType[idx - 1] <= AIR ? 1 : 0;
                float bottomType = fluid_attributes.cellType[idx + 1] <= AIR ? 1 : 0;

                float divideBy = leftType + rightType + topType + bottomType;
            
                if (divideBy == 0.f) continue;

                float divergence = fluid_attributes.u[idx + n] - fluid_attributes.u[idx] + fluid_attributes.v[idx + 1] - fluid_attributes.v[idx];
                if (fluid_attributes.particleRestDensity > 0.f) {
                    float compression = fluid_attributes.cellDensities[idx] - fluid_attributes.particleRestDensity;
                    if (compression > 0.f) {
                        divergence -= k * compression;
                   }
                }

                float p = divergence / divideBy;
                p *= overRelaxation;

                fluid_attributes.u[idx]     += leftType * p;
                fluid_attributes.u[idx + n] -= rightType * p;
                fluid_attributes.v[idx]     += topType * p;
                fluid_attributes.v[idx + 1] -= bottomType * p;
            }
        }
    }

    void SolvePressure() {
        for (int iter = 0; iter < numPressureIters; ++iter) {
            passRedBlackSOR(0);
            passRedBlackSOR(1);
        }
    }
    // --------------------------------------------------------------------------------------------------------------------------------------------


    // Soon to be Multigrid code
    // --------------------------------------------------------------------------------------------------------------------------------------------
    void restrict() {
        
    }

    void prolongate() {

    }

    void projectVCycle() {

    }
    // --------------------------------------------------------------------------------------------------------------------------------------------



    // conjugate gradient code
    // --------------------------------------------------------------------------------------------------------------------------------------------

    /*void setUpResidual() {
        for (int32_t i = 1; i < fluid_attributes.numX - 1; ++i) {
            for (int32_t j = 1; j < fluid_attributes.numY - 1; ++j) {
                int32_t idx = i * n + j;
                if (fluid_attributes.cellType[idx] != FLUID)  {
                    continue;
                }

                float divergence = -(fluid_attributes.u[idx + n] - fluid_attributes.u[idx] + fluid_attributes.v[idx + 1] - fluid_attributes.v[idx]);

                if (fluid_attributes.particleRestDensity > 0.f) {
                    float compression = fluid_attributes.cellDensities[idx] - fluid_attributes.particleRestDensity;
                    divergence += (compression > 0.f) * this->k * compression;
                }

                residual[idx] = divergence;

            }
        }
    }

    void ScaledAdd(std::vector<double> &a, std::vector<double> &b, double c) {
        for (int i = 1; i < fluid_attributes.numX - 1; ++i) {
            if (fluid_attributes.cellType[i] == FLUID) {
                a[i] += b[i] * c;
            }
        }
    }

    double Dot(std::vector<double> &a, std::vector<double> &b) {
        double res = 0.0;
        for (int i = 1; i < fluid_attributes.numX - 1; ++i) {
            if (fluid_attributes.cellType[i] == FLUID) {
                res += a[i] * b[i];
            }
        }
        return res;
    }

    void EqualsPlusTimes(std::vector<double> &a, std::vector<double> &b, double c) {
        for (int i = 1; i < fluid_attributes.numX - 1; ++i) {
            if (fluid_attributes.cellType[i] == FLUID) {
                a[i] = b[i] + a[i] * c;
            }
        }
    }

    void applyPressure() {
        //const float density = 1000.f;
        //const float scale = dt / (density * cellSpacing);
    
        for (int i = 1; i < fluid_attributes.numX - 1; ++i) {
            for (int j = 1; j < fluid_attributes.numY - 1; ++j) {
                int idx = i * n + j;
                int leftIdx =  idx - n;
                int upIdx = idx - 1;

                if ((fluid_attributes.cellType[idx] == FLUID || fluid_attributes.cellType[leftIdx] == FLUID) &&
                    fluid_attributes.cellType[idx] != SOLID && fluid_attributes.cellType[leftIdx] != SOLID) {
                    float p = pressure[idx];
                    float pLeft = pressure[leftIdx];
                    fluid_attributes.u[idx] -= 1 * (p - pLeft);
                }
    
                if ((fluid_attributes.cellType[idx] == FLUID || fluid_attributes.cellType[upIdx] == FLUID) &&
                    fluid_attributes.cellType[idx] != SOLID && fluid_attributes.cellType[upIdx] != SOLID) {
                    float p = pressure[idx];
                    float pTop = pressure[upIdx];
                    fluid_attributes.v[idx] -= 1 * (p - pTop);
                }
            }
        }
    }

    uint8_t getMaterialFromDirection(uint8_t direction, int idx) {
        switch (direction) {
            case LEFT:
                return fluid_attributes.cellType[idx - n];
            case RIGHT:
                return fluid_attributes.cellType[idx + n];
            case BOTTOM:
                return fluid_attributes.cellType[idx + 1];
            case TOP:
                return fluid_attributes.cellType[idx - 1];
        }

        return SOLID; // should never get executed
    }

    uint8_t updateNbrFromNeighbor(uint8_t material, uint8_t nbr_info, uint8_t dir) {

        // if neighbor cell is fluid or air, add to the first 3 bits
        if (material != SOLID) {
            nbr_info++;
        }

        // if neighbor cell is not fluid, end here
        if (material != FLUID) {
            return nbr_info;
        }

        // else, add the direction onto the last 4 bits of nbr_info. Using OR here because we don't wanna touch the first 3 bits
        return nbr_info | dir;
    }

    void setUpNeighbors() {
        for (int i = 1; i < fluid_attributes.numX - 1; ++i) {
            for (int j = 1; j < fluid_attributes.numY - 1; ++j) {
                int idx = i * n + j;
                if (fluid_attributes.cellType[idx] != FLUID) continue;

                uint8_t nbr_info = 0u;
                for (uint8_t dir : directions) {
                    uint8_t material = getMaterialFromDirection(dir, idx);
                    nbr_info = updateNbrFromNeighbor(material, nbr_info, dir);
                }

                neighbors[idx] = nbr_info;
            }
        }
    }

    void MatVec() {
        for (int i = 1; i < fluid_attributes.numX - 1; ++i) {
            for (int j = 1; j < fluid_attributes.numY - 1; ++j) {
                int idx = i * n + j;

                uint8_t nbrs = neighbors[idx];
                if (!nbrs) { // if a cell has no neighbors continue
                    Adiag[idx] = 0.0;
                    continue;
                }

                // all of these & statements are just looking at different bits in nbrs
                // nbrs & CENTER just extracts the first 3 bits of nbrs; if nbrs is 1001 010, then nbrs & center = 1001 010 & 0000 111 = 0000 010
                Adiag[idx] = 
                    ((nbrs & CENTER) * direction[idx]) -
                    ((nbrs & LEFT) ? direction[idx - n] : 0) -
                    ((nbrs & RIGHT) ? direction[idx + n] : 0) -
                    ((nbrs & BOTTOM) ? direction[idx + 1] : 0) -
                    ((nbrs & TOP) ? direction[idx - 1] : 0);
            }
        }
    }

    // conjugate gradient
    void projectCG() {

        setUpNeighbors();

        std::fill(begin(pressure), end(pressure), 0.0);
        std::fill(begin(Adiag), end(Adiag), 0.0);

        setUpResidual();

        std::copy(begin(residual), end(residual), begin(direction));

        double sigma = DotMulti(residual, residual);
        
        for (int iter = 0; iter < numPressureIters && sigma > 0; ++iter) {
            MatVecMulti();
            
            double alpha = sigma / DotMulti(direction, Adiag);
            ScaledAddMulti(pressure, direction, alpha);
            ScaledAddMulti(residual, Adiag, -alpha);
            
            double sigmaOld = sigma;
            sigma = DotMulti(residual, residual);
            EqualsPlusTimesMulti(direction, residual, sigma / sigmaOld);
        }

        applyPressureMulti();
    }

    // general preconditioned conjugate gradient (using matrixless GS/Jacobi methods, so no z needed)
    void projectPCG(std::function<void()> precondition) {

        setUpNeighbors();

        std::fill(begin(pressure), end(pressure), 0.0);
        std::fill(begin(Adiag), end(Adiag), 0.0);

        setUpResidual();

        precondition();

        std::copy(begin(residual), end(residual), begin(direction));

        double sigma = DotMulti(residual, residual);
        
        for (int iter = 0; iter < numPressureIters && sigma > 0; ++iter) {
            MatVecMulti();
            
            double alpha = sigma / DotMulti(direction, Adiag);
            ScaledAddMulti(pressure, direction, alpha);
            ScaledAddMulti(residual, Adiag, -alpha);

            precondition();
            
            double sigmaOld = sigma;
            sigma = DotMulti(residual, residual);
            EqualsPlusTimesMulti(direction, residual, sigma / sigmaOld);
        }

        applyPressureMulti();
    }*/
    // --------------------------------------------------------------------------------------------------------------------------------------------
};