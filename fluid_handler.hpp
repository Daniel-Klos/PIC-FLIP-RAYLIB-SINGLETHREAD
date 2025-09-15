#pragma once
#include <raylib.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <random>
#include <array>

#include "collision_grid.hpp"
#include "linear_solver.hpp"
#include "density_projection.hpp"
#include "transfer_grid.hpp"
#include "fluid_rendering.hpp"
#include "scene_renderer.hpp"
#include "utils.hpp"

class FluidHandler {
    float moveDist;
    float checkSeparationDist;

public:

    float scalingFactor;
    int32_t scaledWIDTH;
    int32_t scaledHEIGHT;

    CollisionGrid collisionGrid;

    FluidState &fluid_attributes;

    float objectSimRadius = 250;
    float objectRenderRadius = objectSimRadius;

    bool forceObjectActive = true;
    bool dragObjectActive = false;
    bool generatorActive = false;

    std::vector<uint32_t> collisions;

    FluidRenderer &fluid_renderer;
    LinearSolver pressure_solver;
    IDPSolver density_solver;
    TransferGrid transfer_grid;

    Color objectColor = {255, 255, 255, 255};

    FluidHandler(FluidState& fas, FluidRenderer &fr): fluid_attributes(fas),
                                                      fluid_renderer(fr),
                                                      pressure_solver(fas),
                                                      density_solver(fas),
                                                      transfer_grid(fas)

    {
        
        this->moveDist = 2 * fluid_attributes.radius;
        this->checkSeparationDist = moveDist * moveDist;

        this->collisions.resize(fluid_attributes.num_particles);

        this->scalingFactor = 2 * fluid_attributes.radius;

        this->scaledWIDTH = std::ceil(static_cast<float>(fluid_attributes.frame_context.WIDTH) / scalingFactor);
        this->scaledHEIGHT = std::ceil(static_cast<float>(fluid_attributes.frame_context.HEIGHT) / scalingFactor);

        collisionGrid = CollisionGrid(scaledWIDTH, scaledHEIGHT);

        objectRenderRadius = objectSimRadius / fluid_attributes.frame_context.zoom_amount;
    }

    void ApplyBodyForces() {
        if (dragObjectActive) {                                                          // dragging
            includeDragObject();
        } if (forceObjectActive && fluid_attributes.frame_context.leftMouseDown) {  // force pull
            includeForceObject(-250);
        } if (forceObjectActive && fluid_attributes.frame_context.rightMouseDown) { // force push
            includeForceObject(1000); // pushing, 1000
        }
        if (fluid_renderer.getRenderPattern() == 3 && fluid_attributes.fireActive) {     // applying buoyant force to high temperature partiles
            makeFire();
        }
        if (fluid_renderer.getRenderPattern() == 3) {                              // ground->particle thermal conducting. I just lump it in here but it's not an actual force
            heatGround();                                                          // note that particle->particle thermal conducting is handled in solveCollisions()
        }
        if (fluid_attributes.vorticityStrength != 0) {           // vorticity confinement
            applyVorticityConfinementRedBlack();
        }
    }

    void UpdateEnvironment() {
        Advect();

        SolveCollisions();

        fluid_attributes.CollideSurfaces();

        MarkAirAndFluidCells();

        density_solver.ComputeGridDensity();

        //density_solver.ProjectDensity();
        //MarkAirAndFluidCells();

        transfer_grid.TransferToGrid();

        ApplyBodyForces();
        
        pressure_solver.SolvePressure();

        transfer_grid.TransferToParticles();
    }

    void Advect() {
        for (int i = 0; i < fluid_attributes.num_particles; ++i) {
            fluid_attributes.positions[2 * i] += fluid_attributes.velocities[2 * i] * fluid_attributes.frame_context.dt;
            fluid_attributes.positions[2 * i + 1] += fluid_attributes.velocities[2 * i + 1] * fluid_attributes.frame_context.dt;
            fluid_attributes.velocities[2 * i] += fluid_attributes.gravityX * fluid_attributes.frame_context.dt;
            fluid_attributes.velocities[2 * i + 1] += fluid_attributes.gravityY * fluid_attributes.frame_context.dt;

            fluid_attributes.particleAges[i] = std::min(fluid_attributes.age_constant, fluid_attributes.particleAges[i] + 1);
        }
    }

    void MarkAirAndFluidCells() {
        for (int i = 0; i < fluid_attributes.numX; ++i) {
            for (int j = 0; j < fluid_attributes.numY; ++j) {
                int idx = i * fluid_attributes.numY + j;
                if (fluid_attributes.cellType[idx] != fluid_attributes.SOLID) {
                    fluid_attributes.cellType[idx] = fluid_attributes.AIR;
                }
            }
        }

        fluid_attributes.num_fluid_cells = 0;
        for (int i = 0; i < fluid_attributes.num_particles; ++i) {
            float x = fluid_attributes.positions[2 * i];
            float y = fluid_attributes.positions[2 * i + 1];

            int xi = clamp(std::floor(x * fluid_attributes.invSpacing), 0, fluid_attributes.numX - 1);
            int yi = clamp(std::floor(y * fluid_attributes.invSpacing), 0, fluid_attributes.numY - 1);

            int cellNr = xi * fluid_attributes.numY + yi;
            if (fluid_attributes.cellType[cellNr] == fluid_attributes.AIR) {
                fluid_attributes.cellType[cellNr] = fluid_attributes.FLUID;
                fluid_attributes.fluid_cells[fluid_attributes.num_fluid_cells] = cellNr;
                fluid_attributes.num_fluid_cells++;
            }
        }
    }

    void heatGround() {
        const float percentRemoved = 0.1f;
        const float leftEdge = fluid_attributes.frame_context.WIDTH * percentRemoved;
        const float rightEdge = fluid_attributes.frame_context.WIDTH * (1.0f - percentRemoved);
        for (int i = 0; i < fluid_attributes.num_particles; ++i) {
            if (fluid_attributes.positions[2 * i + 1] + fluid_attributes.radius > fluid_attributes.frame_context.HEIGHT - 2 * fluid_attributes.cellSpacing) {
                if ((fluid_attributes.temperatures[i] < fluid_renderer.tempGradient.size()) || (fluid_attributes.positions[2 * i] > leftEdge && fluid_attributes.positions[2 * i] < rightEdge)) {
                    fluid_attributes.temperatures[i] += fluid_attributes.groundConductivity * fluid_attributes.frame_context.dt;
                }
            }
        }
    }

    // loop over all cells
        // calculate an aggregate temperature for each cell (for the cell center)
        // by looping over all particles and calculating mean particle temperature or something
    // loop over all cells again
        // for each particle, lerp a value for surround_temp from surrounding cells
            // remember to offset:
                // int x = pos.x;
                // int y = pos.y;

                // float gx = x / fas.cellSpacing - 0.5f;
                // float gy = y / fas.cellSpacing - 0.5f;
            // temperatures[i] += surround_temp * fluid_attributes.interConductivity * dt
    /*void DiffuseHeat() {
        for (int i = 0; i < fluid_attributes.num_fluid_cells; ++i) {
            auto &cell = fluid_attributes.cellOccupants.data[i];

            // loop over all cells and calculate an aggregate temperature (for the cell center)
            
            for (int id = 0; id < cell.objects_count; ++id) {
                int pIdx = cell.objects[id];

                const float transfer = (fluid_attributes.temperatures[index] - fluid_attributes.temperatures[otherIndex]) * fluid_attributes.interConductivity * fluid_attributes.frame_context.dt;
                
                fluid_attributes.temperatures[index] -= transfer * fluid_attributes.frame_context.dt;
                fluid_attributes.temperatures[otherIndex] += transfer * fluid_attributes.frame_context.dt;

                collisions[index]++;
                collisions[otherIndex]++;
        }
    }*/

    void makeFire() {
        for (int i = 0; i < fluid_attributes.num_particles; ++i) {
            if (fluid_attributes.positions[2 * i + 1] < fluid_attributes.frame_context.HEIGHT - fluid_attributes.cellSpacing - 10) {
                fluid_attributes.velocities[2 * i + 1] -= fluid_attributes.fireStrength * fluid_attributes.temperatures[i] * fluid_attributes.frame_context.dt;
                if (fluid_attributes.temperatures[i] > 0) {
                    fluid_attributes.temperatures[i] -= fluid_attributes.tempDiffusion * fluid_attributes.frame_context.dt;
                }
                if (collisions[i] == 0) {
                    fluid_attributes.temperatures[i] = 0;
                }
            }
        }
    }

    void FillCollisionGrid()
    {

        collisionGrid.clear();

        const float minX = fluid_attributes.cellSpacing;
        const float maxX = fluid_attributes.frame_context.WIDTH - fluid_attributes.cellSpacing;
        const float minY = fluid_attributes.cellSpacing;
        const float maxY = fluid_attributes.frame_context.HEIGHT - fluid_attributes.cellSpacing;

        uint32_t i = 0;

        for (int32_t index = 0; index < fluid_attributes.num_particles; ++index) {
            float x = fluid_attributes.positions[2 * index];
            float y = fluid_attributes.positions[2 * index + 1];
            if (x > minX && x < maxX && y > minY && y < maxY) {
                int32_t gridX = x / scalingFactor;
                int32_t gridY = y / scalingFactor;
                collisionGrid.addAtom(gridX, gridY, i);
            }
            ++i;
        }
    }

    void solveContact(uint32_t index, uint32_t otherIndex)
    {
        constexpr float eps = 0.0001f;
        const float o2_o1X = fluid_attributes.positions[2 * index] - fluid_attributes.positions[2 * otherIndex];
        const float o2_o1Y = fluid_attributes.positions[2 * index + 1] - fluid_attributes.positions[2 * otherIndex + 1];

        const float dist2 = o2_o1X * o2_o1X + o2_o1Y * o2_o1Y;

        if (dist2 < checkSeparationDist && dist2 > eps) {
            const float dist          = sqrt(dist2);
            const float delta = 0.5f * (moveDist - dist) / dist;
            const float col_vecX = o2_o1X * delta;
            const float col_vecY = o2_o1Y * delta;

            fluid_attributes.positions[2 * index] += col_vecX;
            fluid_attributes.positions[2 * index + 1] += col_vecY;

            fluid_attributes.positions[2 * otherIndex] -= col_vecX;
            fluid_attributes.positions[2 * otherIndex + 1] -= col_vecY;

            const float transfer = (fluid_attributes.temperatures[index] - fluid_attributes.temperatures[otherIndex]) * fluid_attributes.interConductivity * fluid_attributes.frame_context.dt;

            fluid_attributes.temperatures[index] -= transfer * fluid_attributes.frame_context.dt;
            fluid_attributes.temperatures[otherIndex] += transfer * fluid_attributes.frame_context.dt;

            collisions[index]++;
            collisions[otherIndex]++;
        }
    }

    void checkCellCollisions(uint32_t atom_idx, const CollisionCell& c)
    {
        for (uint32_t i = 0; i < c.objects_count; ++i) {
            solveContact(atom_idx, c.objects[i]);
        }
    }

    void processCell(const CollisionCell& c, uint32_t index)
    {
        for (uint32_t i = 0; i < c.objects_count; ++i) {
            const uint32_t idx = c.objects[i];
            for (int32_t side = 0; side < 2; ++side) {
                checkCellCollisions(idx, collisionGrid.data[index + collisionGrid.height + side]);
            }
            for (int32_t side = 0; side < 2; ++side) {
                checkCellCollisions(idx, collisionGrid.data[index + side]);   
            }
            //checkCellCollisions(idx, grid.data[index - grid.height]);
            checkCellCollisions(idx, collisionGrid.data[index - collisionGrid.height + 1]);
        }
    }

    void SolveCollisions()
    {
        if (fluid_attributes.fireActive) {
            std::fill(begin(collisions), end(collisions), 0);
        }

        FillCollisionGrid();

        for (int32_t idx = 0; idx < collisionGrid.width * collisionGrid.height; ++idx) {
            processCell(collisionGrid.data[idx], idx);
        }
    }

    void calcVorticityConfinement(bool red, int32_t startColumn, int32_t endColumn) {
        for (int32_t i = startColumn; i < endColumn; ++i) {
            for (int32_t j = 1; j < fluid_attributes.numY - 1; ++j) {
                if (red) {
                    if ((i + j) % 2 != 0) continue; 
                }
                else {
                    if ((i + j) % 2 == 0) continue; 
                }
                float dx = abs(fluid_attributes.curl(i, j + 1)) - abs(fluid_attributes.curl(i, j - 1));
                float dy = abs(fluid_attributes.curl(i + 1, j)) - abs(fluid_attributes.curl(i - 1, j));

                const float len = std::sqrt(dx * dx + dy * dy);
                const float invLen = 1.f / (len + (len == 0.f)) - (len == 0.f);

                dx *= invLen;
                dy *= invLen;

                const float c = fluid_attributes.curl(i, j);

                fluid_attributes.v[i * fluid_attributes.n + j] += c * dx * fluid_attributes.frame_context.dt * fluid_attributes.vorticityStrength;
                fluid_attributes.u[i * fluid_attributes.n + j] += c * dy * fluid_attributes.frame_context.dt * fluid_attributes.vorticityStrength;
            }
        }
    }

    void applyVorticityConfinementRedBlack() {
        this->calcVorticityConfinement(true, 1, fluid_attributes.numX - 1);
        this->calcVorticityConfinement(false, 1, fluid_attributes.numX - 1);
    }

    void includeForceObject(const int32_t strength) {
        const float extend = (20 * fluid_attributes.cellSpacing) / fluid_attributes.frame_context.zoom_amount;

        int forceObjectX = fluid_attributes.frame_context.world_mouse_pos.x;
        int forceObjectY = fluid_attributes.frame_context.world_mouse_pos.y;

        int objectCellX = forceObjectX / fluid_attributes.cellSpacing;
        int objectCellY = forceObjectY / fluid_attributes.cellSpacing;
        int objectCellRadius = std::ceil(objectSimRadius / fluid_attributes.cellSpacing);

        for (int i = objectCellX - objectCellRadius; i < objectCellX + objectCellRadius; i++) {
            for (int j = objectCellY - objectCellRadius; j < objectCellY + objectCellRadius; j++) {
                int cellNr = i * fluid_attributes.n + j;
                bool notInBounds = i < 0 || i >= fluid_attributes.numX || j < 0 || j >= fluid_attributes.numY;
                if (notInBounds || (fluid_attributes.cellType[i * fluid_attributes.n + j] == fluid_attributes.SOLID)) continue;
                float dx = (i + 0.5) * fluid_attributes.cellSpacing - forceObjectX;
                float dy = (j + 0.5) * fluid_attributes.cellSpacing - forceObjectY;
                float d2 = dx * dx + dy * dy;

                if (d2 > objectSimRadius * objectSimRadius + extend || d2 == 0.0f) continue;

                float d = std::sqrt(d2);

                float edgeT = d / objectSimRadius;
            
                float centerT = 1 - edgeT;
                    
                if (fluid_attributes.cellType[cellNr - fluid_attributes.n] != fluid_attributes.SOLID) {
                    fluid_attributes.u[cellNr] += (dx * strength - fluid_attributes.u[cellNr]) * centerT * fluid_attributes.frame_context.dt;
                }
                if (fluid_attributes.cellType[cellNr + fluid_attributes.n] != fluid_attributes.SOLID) {
                    fluid_attributes.u[cellNr + fluid_attributes.n] += (dx * strength - fluid_attributes.u[cellNr + fluid_attributes.n]) * centerT * fluid_attributes.frame_context.dt;
                }
                if (fluid_attributes.cellType[cellNr - 1] != fluid_attributes.SOLID) {
                    fluid_attributes.v[cellNr] += (dy * strength - fluid_attributes.v[cellNr]) * centerT * fluid_attributes.frame_context.dt;
                }
                if (fluid_attributes.cellType[cellNr + 1] != fluid_attributes.SOLID) {
                    fluid_attributes.v[cellNr + 1] += (dy * strength - fluid_attributes.v[cellNr + 1]) * centerT * fluid_attributes.frame_context.dt;
                }
            }
        }
    }

    void includeDragObject() {
        // Bug when you hold force object attract then switch to drag object while holding mouse down
        const float extend = (20 * fluid_attributes.cellSpacing) / fluid_attributes.frame_context.zoom_amount;
        if (fluid_attributes.frame_context.leftMouseDown || fluid_attributes.frame_context.rightMouseDown) {
            float objectX = fluid_attributes.frame_context.world_mouse_pos.x;
            float objectY = fluid_attributes.frame_context.world_mouse_pos.y;
            float objectPrevX = fluid_attributes.frame_context.prev_world_mouse_pos.x;
            float objectPrevY = fluid_attributes.frame_context.prev_world_mouse_pos.y;
            float vx = (objectX - objectPrevX) * fluid_attributes.frame_context.maxFps * !fluid_attributes.frame_context.justPressed;
            float vy = (objectY - objectPrevY) * fluid_attributes.frame_context.maxFps * !fluid_attributes.frame_context.justPressed;
            int objectCellX = objectX / fluid_attributes.cellSpacing;
            int objectCellY = objectY / fluid_attributes.cellSpacing;
            int objectCellRadius = std::ceil(objectSimRadius / fluid_attributes.cellSpacing);
            for (int i = objectCellX - objectCellRadius; i < objectCellX + objectCellRadius; i++) {
                for (int j = objectCellY - objectCellRadius; j < objectCellY + objectCellRadius; j++) {
                    int cellNr = i * fluid_attributes.n + j;
                    bool notInBounds = i < 0 || i >= fluid_attributes.numX || j < 0 || j >= fluid_attributes.numY;
                    if (notInBounds || (fluid_attributes.cellType[i * fluid_attributes.n + j] == fluid_attributes.SOLID)) continue;
                    float dx = (i + 0.5) * fluid_attributes.cellSpacing - objectX;
                    float dy = (j + 0.5) * fluid_attributes.cellSpacing - objectY;

                    if (dx * dx + dy * dy < objectSimRadius * objectSimRadius + extend) {
                        
                        if (fluid_attributes.cellType[cellNr - fluid_attributes.n] != fluid_attributes.SOLID) {
                            fluid_attributes.u[cellNr] = vx;
                        }
                        if (fluid_attributes.cellType[cellNr + fluid_attributes.n] != fluid_attributes.SOLID) {
                            fluid_attributes.u[cellNr + fluid_attributes.n] = vx;
                        }
                        if (fluid_attributes.cellType[cellNr - 1] != fluid_attributes.SOLID) {
                            fluid_attributes.v[cellNr] = vy;
                        }
                        if (fluid_attributes.cellType[cellNr + 1] != fluid_attributes.SOLID) {
                            fluid_attributes.v[cellNr + 1] = vy;
                        }
                    }
                }
            }
        }
    }

    void generate() {
        float separation = fluid_attributes.radius * 2.f;
        int wideNum = std::floor((2 * objectSimRadius) / (separation));
        int highNum = wideNum;
        int numPotentiallyAdded = wideNum * highNum;

        float generatorRightBound = fluid_attributes.frame_context.world_mouse_pos.x + objectSimRadius;
        float generatorBottomBound = fluid_attributes.frame_context.world_mouse_pos.y + objectSimRadius;

        float starting_px = std::max(fluid_attributes.frame_context.world_mouse_pos.x - objectSimRadius + fluid_attributes.radius, fluid_attributes.cellSpacing + fluid_attributes.radius);
        float starting_py = std::max(fluid_attributes.frame_context.world_mouse_pos.y - objectSimRadius + fluid_attributes.radius, fluid_attributes.cellSpacing + fluid_attributes.radius);
        float px = starting_px;
        float py = starting_py;
        bool offset = true;

        int originalParticleCount = fluid_attributes.num_particles;
        int maxPossibleParticles = originalParticleCount + numPotentiallyAdded;

        fluid_attributes.positions.resize(2 * maxPossibleParticles);
        fluid_attributes.renderPositions.resize(2 * maxPossibleParticles);
        int addedTo = 0;

        for (int i = 0; i < numPotentiallyAdded; ++i) {
            float prevPx = px;
            float prevPy = py;
            if (prevPy > generatorBottomBound || prevPy + fluid_attributes.radius > fluid_attributes.frame_context.HEIGHT - fluid_attributes.cellSpacing) {
                break;
            }

            int cellX = px / fluid_attributes.cellSpacing;
            int cellY = py / fluid_attributes.cellSpacing;
            int cellNr = cellX * fluid_attributes.n + cellY;

            px += separation;
            if (px > generatorRightBound) {
                px = starting_px;
                /*if (offset) {
                    px += 0.5 * separation;
                }*/
                py += separation;
                offset = !offset;
            }

            if (cellX > fluid_attributes.numX || (fluid_attributes.cellType[cellNr] != fluid_attributes.AIR || prevPx - fluid_attributes.radius < fluid_attributes.cellSpacing || prevPx + fluid_attributes.radius > fluid_attributes.frame_context.WIDTH - fluid_attributes.cellSpacing || prevPy - fluid_attributes.radius < fluid_attributes.cellSpacing)) {
                continue;
            }


            fluid_attributes.positions[2 * (fluid_attributes.num_particles + addedTo)] = prevPx;
            fluid_attributes.positions[2 * (fluid_attributes.num_particles + addedTo) + 1] = prevPy;
            
            fluid_attributes.renderPositions[2 * (fluid_attributes.num_particles + addedTo)] = prevPx;
            fluid_attributes.renderPositions[2 * (fluid_attributes.num_particles + addedTo) + 1] = prevPy;

            addedTo++;
        }

        fluid_attributes.num_particles += addedTo;

        fluid_attributes.positions.resize(2 * fluid_attributes.num_particles);
        fluid_attributes.renderPositions.resize(2 * fluid_attributes.num_particles);
        fluid_attributes.velocities.resize(2 * fluid_attributes.num_particles);
        fluid_renderer.particleDiffusionColors.resize(3 * fluid_attributes.num_particles);

        fluid_attributes.particleAges.resize(fluid_attributes.num_particles);
        //fluid_attributes.debug.resize(fluid_attributes.num_particles);

        int start = fluid_attributes.num_particles - addedTo;

        for (int i = start; i < fluid_attributes.num_particles; ++i) {
            int idx1 = 2 * i;
            int idx2 = 3 * i;

            fluid_attributes.velocities[idx1] = 0.f;
            fluid_attributes.velocities[idx1 + 1] = 0.f;

            fluid_renderer.particleDiffusionColors[idx2] = 0;
            fluid_renderer.particleDiffusionColors[idx2 + 1] = 0;
            fluid_renderer.particleDiffusionColors[idx2 + 2] = 255;

            fluid_attributes.particleAges[i] = 0;
        }

        // make Resize(num_particles) methods within each struct. this is too messy

        this->collisions.resize(fluid_attributes.num_particles);
        fluid_attributes.temperatures.resize(fluid_attributes.num_particles);

        transfer_grid.Resize(fluid_attributes.num_particles);

        fluid_renderer.ResizeAndUpdateMesh(fluid_attributes.num_particles);
    }

    void remove() {
        float spacing = fluid_attributes.cellSpacing;
        const int32_t numCovered = std::ceil(objectSimRadius / spacing);
        const int32_t mouseColumn = std::floor(fluid_attributes.frame_context.world_mouse_pos.x / spacing);
        const int32_t mouseRow = std::floor(fluid_attributes.frame_context.world_mouse_pos.y / spacing);
        const size_t len = fluid_attributes.num_particles;
        const size_t double_len = 2 * len;
        const size_t triple_len = 3 * len;

        size_t remove = 0;
        
        for (int32_t i = -numCovered; i < numCovered + 1; ++i) {
            for (int32_t j = -numCovered; j < numCovered + 1; ++j) {
                if (mouseRow + j <= 1 || mouseRow + j >= fluid_attributes.numY - 1 || mouseColumn + i <= 1 || mouseColumn + i >= fluid_attributes.numX - 1)
                    continue;
    
                const auto &cell = fluid_attributes.cellOccupants.data[mouseRow + j + fluid_attributes.cellOccupants.height * (mouseColumn + i)];
    
                for (uint32_t id{0}; id < cell.objects_count; ++id) {
                    const uint32_t particleIndex = cell.objects[id];

                    size_t doubleid = 2 * particleIndex;
                    size_t tripleid = 3 * particleIndex;

                    size_t doubleremove = 2 * remove;
                    size_t tripleremove = 3 * remove;
    
                    if (doubleid + 2 < double_len) {
                        /*fluid_attributes.particle_densities[particleIndex] = fluid_attributes.particle_densities[particleIndex - 1 - remove];
                        fluid_attributes.densities[particleIndex] = fluid_attributes.particle_densities[particleIndex - 1 - remove];*/

                        fluid_attributes.temperatures[particleIndex] = fluid_attributes.temperatures[len - 1 - remove];

                        fluid_renderer.replaceParticleWithEnd(particleIndex, len - 1 - remove);

                        fluid_attributes.positions[doubleid] = fluid_attributes.positions[double_len - 2 - doubleremove];
                        fluid_attributes.positions[doubleid + 1] = fluid_attributes.positions[double_len - 1 - doubleremove];

                        fluid_attributes.renderPositions[doubleid] = fluid_attributes.renderPositions[double_len - 2 - doubleremove];
                        fluid_attributes.renderPositions[doubleid + 1] = fluid_attributes.renderPositions[double_len - 1 - doubleremove];

                        fluid_attributes.velocities[doubleid] = fluid_attributes.velocities[double_len - 2 - doubleremove];
                        fluid_attributes.velocities[doubleid + 1] = fluid_attributes.velocities[double_len - 1 - doubleremove];

                        fluid_renderer.particleDiffusionColors[tripleid] = fluid_renderer.particleDiffusionColors[triple_len - 3 - tripleremove];
                        fluid_renderer.particleDiffusionColors[tripleid + 1] = fluid_renderer.particleDiffusionColors[triple_len - 2 - tripleremove];
                        fluid_renderer.particleDiffusionColors[tripleid + 2] = fluid_renderer.particleDiffusionColors[triple_len - 1 - tripleremove];
                    }
                    remove++;
                }
            }
        }
        fluid_attributes.num_particles -= remove;
    
        fluid_attributes.positions.resize(2 * fluid_attributes.num_particles);
        fluid_attributes.renderPositions.resize(2 * fluid_attributes.num_particles);
        fluid_attributes.velocities.resize(2 * fluid_attributes.num_particles);
        //fluid_attributes.debug.resize(fluid_attributes.num_particles);
    
        this->collisions.resize(fluid_attributes.num_particles);
        fluid_attributes.temperatures.resize(fluid_attributes.num_particles);

        transfer_grid.Resize(fluid_attributes.num_particles);

        fluid_attributes.particleAges.resize(fluid_attributes.num_particles);

        fluid_renderer.ResizeAndUpdateMesh(fluid_attributes.num_particles);
    }

    void drawCircleObject() {
        DrawCircleLinesV(fluid_attributes.frame_context.screen_mouse_pos, objectRenderRadius, objectColor);
    }

    void drawSquareObject() {
        DrawRectangleLines(fluid_attributes.frame_context.screen_mouse_pos.x - objectRenderRadius, fluid_attributes.frame_context.screen_mouse_pos.y - objectRenderRadius, 2 * objectRenderRadius, 2 * objectRenderRadius, objectColor);
    }

    void render_objects() {
        if (forceObjectActive || dragObjectActive) {
            drawCircleObject();
        }
        if (generatorActive) {
            drawSquareObject();
        }
    }

    void addToObjectRadius(float add) {
        if (objectRenderRadius + add > 0) {
            objectRenderRadius += add;

            objectSimRadius = objectRenderRadius / fluid_attributes.frame_context.zoom_amount;
        }
    }

    bool getObjectsActive() {
        return this->dragObjectActive || this->forceObjectActive || this->generatorActive;
    }

    bool getDragObjectActive() {
        return this->dragObjectActive;
    }

    void setDragObjectActive(bool set) {
        this->dragObjectActive = set;
    }

    bool getForceObjectActive() {
        return this->forceObjectActive;
    }

    void setForceObjectActive(bool active) {
        this->forceObjectActive = active;
    }

    bool getGeneratorActive() {
        return this->generatorActive;
    }

    void setGeneratorActive(bool active) {
        this->generatorActive = active;
    }

};