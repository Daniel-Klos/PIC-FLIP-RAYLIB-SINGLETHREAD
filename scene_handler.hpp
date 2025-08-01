#pragma once
#include <sstream>
#include <iomanip>
#include "simulation_state.hpp"
#include "fluid_handler.hpp"
#include "obstacle_handler.hpp"
#include "scene_renderer.hpp"

#include <rlgl.h>

struct SceneHandler {
    FluidState &fluid_attributes;

    SceneRenderer scene_renderer;
    ObstacleHandler obstacle_handler;
    FluidHandler fluid_handler;

    float mouseX;
    float mouseY;

    int frame = 0;
    float trueDT;
    float setDT;
    int fps;

    float steps = 0.f;
    float CollisionTime = 0.f;
    float ObstacleCollisionTime = 0.f;
    float ToGridTime = 0.f;
    float DensityUpdateTime = 0.f;
    float ProjectionTime = 0.f;
    float ToParticlesTime = 0.f;
    float RenderingTime = 0.f;
    float FillGridTime = 0.f;
    float SimStepTime = 0.f;
    float miscellaneousTime = 0.f;

    std::string showControlsStr = "Show Controls";
    std::string hideControlsStr = "Hide Controls";

    int showControlsX;
    int showControlsWIDTH;
    int showControlsY;
    int showControlsHEIGHT;

    int hideControlsX;  
    int hideControlsWIDTH;
    int hideControlsY;
    int hideControlsHEIGHT;

    bool showControls = false;

    float fontSize;

    SceneHandler(FluidState &fas, int maxFps)
        : fluid_attributes(fas),
          scene_renderer(fas), 
          obstacle_handler(fas, scene_renderer.obstacle_renderer), 
          fluid_handler(fas, scene_renderer.fluid_renderer)
    {
        fluid_attributes.frame_context.maxFps = maxFps;
        setDT = 1.f / maxFps;
        SetTargetFPS(maxFps);

        // initialize obstacle positions
        // ----------------------------------------------------------------------------------------------------------------------------
        size_t numInitObstacles = 2 * fluid_attributes.numX + 2 * (fluid_attributes.numY - 2);
        fluid_attributes.obstaclePositions.resize(numInitObstacles);
        int idx = 0;
        for (int i = 0; i < fluid_attributes.numX; ++i) {
            fluid_attributes.cellType[i * fluid_attributes.n] = fluid_attributes.SOLID_CELL;
            fluid_attributes.obstaclePositions[idx] = Vector2i{i, 0};
            ++idx;
        }
        for (int i = 0; i < fluid_attributes.numX; ++i) {
            fluid_attributes.cellType[i * fluid_attributes.n + fluid_attributes.numY - 1] = fluid_attributes.SOLID_CELL;
            fluid_attributes.obstaclePositions[idx] = Vector2i{i, fluid_attributes.numY - 1};
            ++idx;
        }
        for (int j = 1; j < fluid_attributes.numY - 1; ++j) {
            fluid_attributes.cellType[j] = fluid_attributes.SOLID_CELL;
            fluid_attributes.obstaclePositions[idx] = Vector2i{0, j};
            ++idx;
        }
        for (int j = 1; j < fluid_attributes.numY - 1; ++j) {
            fluid_attributes.cellType[fluid_attributes.numX * fluid_attributes.numY + j - fluid_attributes.n] = fluid_attributes.SOLID_CELL;
            fluid_attributes.obstaclePositions[idx] = Vector2i{fluid_attributes.numX - 1, j};
            ++idx;
        }
        // ----------------------------------------------------------------------------------------------------------------------------


        
        // initialize obstacle_renderer
        // ----------------------------------------------------------------------------------------------------------------------------
        scene_renderer.obstacle_renderer.ResizeAndUploadNewObstacleMesh(numInitObstacles);
        // ----------------------------------------------------------------------------------------------------------------------------

        float WIDTH = fluid_attributes.frame_context.WIDTH;

        fontSize = GetFractionOf(WIDTH, 0.02f);

        showControlsWIDTH = MeasureText(showControlsStr.c_str(), fontSize);
        showControlsHEIGHT = fontSize;
        showControlsX = WIDTH - showControlsWIDTH - 2 * fluid_attributes.cellSpacing;
        showControlsY = showControlsHEIGHT;

        hideControlsWIDTH = MeasureText(hideControlsStr.c_str(), fontSize);
        hideControlsHEIGHT = fontSize;
        hideControlsX = WIDTH - hideControlsWIDTH - 2 * fluid_attributes.cellSpacing;
        hideControlsY = 14 * hideControlsHEIGHT;
    }

    void simulate() {
        trueDT = GetFrameTime();
        trueDT = std::max(trueDT, setDT);

        Vector2 mouse_position = GetMousePosition();

        fluid_attributes.frame_context.screen_mouse_pos = Vector2{static_cast<float>(mouse_position.x), static_cast<float>(mouse_position.y)};
        fluid_attributes.frame_context.world_mouse_pos = scene_renderer.screenToWorld(fluid_attributes.frame_context.screen_mouse_pos);

        fluid_attributes.frame_context.dt = setDT;

        track_key_events();
        track_mouse_events();
            
        if (!fluid_attributes.stop || fluid_attributes.step) {
            update_environment();
            fluid_attributes.step = false;
        }

        handle_zoom();

        BeginDrawing();

        ClearBackground(BLACK);

        scene_renderer.render_scene();

        render_interaction_objects();

        displayGUI();

        EndDrawing();

        fluid_attributes.frame_context.prev_screen_mouse_pos = fluid_attributes.frame_context.screen_mouse_pos;
        fluid_attributes.frame_context.prev_world_mouse_pos = fluid_attributes.frame_context.world_mouse_pos;

        fluid_attributes.frame_context.zooming = false;
        
    }

    void update_environment() {
        ++steps;

        // order of need of implementation/optimization:
            // 1) move event handling into scene_handler
            // 2) add a way to switch between particle view and liquid glass view
            // 3) Thread buffers EVERYWHERE
            // 4) implement Multi Grid, then MGPCG
            // 5) make a sampleVelocity(point) function that can easily incorporate different interpolation functions
                // 6) so that you can do RK advection easier
            // 7) DDA raycasting & making sure that particle-particle collisions dont push particles into obstacles
            // 8) implement implicit density projection
            // 9) finish cleaning up all this code
            // 10) make it so that you pass in static arrays instead of just numbers of particles in the main file

        fluid_handler.integrate();

        // fluid_handler.pressure_solver.project_density_implicit();

        if (fluid_handler.fluid_renderer.getRenderPattern() == 4 && fluid_attributes.fireActive) {
            fluid_handler.makeFire();
        }

        if (fluid_attributes.fireActive) {
            std::fill(begin(fluid_handler.collisions), end(fluid_handler.collisions), 0);
        }

        fluid_handler.addObjectsToGrids();


        if (obstacle_handler.solidDrawing && fluid_attributes.frame_context.leftMouseDown) {
            obstacle_handler.drawSolids();
        } else if (obstacle_handler.solidDrawing && fluid_attributes.frame_context.rightMouseDown) {
            obstacle_handler.eraseSolids();
        }

        if (fluid_handler.generatorActive && fluid_attributes.frame_context.leftMouseDown) {
            fluid_handler.generate();
        } else if (fluid_handler.generatorActive && fluid_attributes.frame_context.rightMouseDown) {
            fluid_handler.remove();
        }


        fluid_handler.solveCollisions();





        if (fluid_handler.fluid_renderer.getRenderPattern() == 4) {
            fluid_handler.heatGround();
        }

        obstacle_handler.collideSurfacesMulti();

        obstacle_handler.constrainWallsMulti();



        fluid_handler.transfer_grid.TransferToGrid();
        
        fluid_handler.transfer_grid.updateCellDensitiesMulti();

        // apply body forces to grid
        if (fluid_handler.dragObjectActive) {
            fluid_handler.includeDragObject();
        } else if (fluid_handler.forceObjectActive && fluid_attributes.frame_context.leftMouseDown) {
            fluid_handler.includeForceObject(-250); // pulling, -250
        } else if (fluid_handler.forceObjectActive && fluid_attributes.frame_context.rightMouseDown) { 
            fluid_handler.includeForceObject(1000); // pushing, 1000
        }

        if (fluid_attributes.vorticityStrength != 0) {
            fluid_handler.applyVorticityConfinementRedBlack();
        }


        //start = std::chrono::high_resolution_clock::now();
        
        fluid_handler.pressure_solver.projectRedBlackSORMulti(fluid_handler.pressure_solver.numPressureIters);

        // PCG with RBSOR as a preconditioner makes super cool high frequency patterns
        /*fluid_handler.pressure_solver.projectPCG([&]() {
            fluid_handler.pressure_solver.projectRedBlackSORMulti(3); // around 2-6 iterations for preconditioner
        });*/
        //fluid_handler.pressure_solver.projectCG();


        fluid_handler.transfer_grid.TransferToParticles();
    }

    void track_key_events() {
        //auto start = std::chrono::high_resolution_clock::now();
        if (IsKeyDown(KEY_B)) {
            if (lteEpsPlus(fluid_attributes.getFlipRatio(), 0.99f)) {
                bool smallIncrement = gteEpsMinus(fluid_attributes.getFlipRatio(), 0.9f);
                bool largeIncrement = ltEpsMinus(fluid_attributes.getFlipRatio(), 0.9f);
                fluid_attributes.addToFlipRatio(0.01 * smallIncrement + 0.1 * largeIncrement);

                if (fluid_attributes.getFlipRatio() > 1.f) {
                    fluid_attributes.setFlipRatio(1.f);
                }
            }
        }
        if (IsKeyDown(KEY_S)) {
            if (gteEpsMinus(fluid_attributes.getFlipRatio(), 0.1f)) {
                bool smallIncrement = gtEpsPlus(fluid_attributes.getFlipRatio(), 0.9f);
                bool largeIncrement = lteEpsPlus(fluid_attributes.getFlipRatio(), 0.9f);
                fluid_attributes.addToFlipRatio(-0.01 * smallIncrement - 0.1 * largeIncrement);

                if (fluid_attributes.getFlipRatio() < 0.f) {
                    fluid_attributes.setFlipRatio(0.f);
                }
            }
        }
        if (IsKeyDown(KEY_E)) {
            fluid_attributes.addToVorticityStrength(10);
        }
        if (IsKeyDown(KEY_W)) {
            if (fluid_attributes.getVorticityStrength() - 10 >= 0.f) {
                fluid_attributes.addToVorticityStrength(-10);
            }
        }
        if (IsKeyDown(KEY_P)) {
            fluid_handler.pressure_solver.addToNumPressureIters(1);
        }
        if (IsKeyDown(KEY_O)) {
            if (fluid_handler.pressure_solver.getNumPressureIters() > 0) {
                fluid_handler.pressure_solver.addToNumPressureIters(-1);
            }
        }
        if (IsKeyDown(KEY_R)) {
            scene_renderer.reset_zoom();
        }
        if (IsKeyDown(KEY_ONE)) {
            fluid_handler.setDragObjectActive(true);
            fluid_handler.setForceObjectActive(false);
            fluid_handler.setGeneratorActive(false);
            obstacle_handler.setSolidDrawer(false);
            scene_renderer.setZoomObjectActive(false);
        }
        if (IsKeyDown(KEY_TWO)) {
            fluid_handler.setDragObjectActive(false);
            fluid_handler.setForceObjectActive(true);
            fluid_handler.setGeneratorActive(false);
            obstacle_handler.setSolidDrawer(false);
            scene_renderer.setZoomObjectActive(false);
        }
        if (IsKeyDown(KEY_THREE)) {
            fluid_handler.setDragObjectActive(false);
            fluid_handler.setForceObjectActive(false);
            fluid_handler.setGeneratorActive(true);
            obstacle_handler.setSolidDrawer(false);
            scene_renderer.setZoomObjectActive(false);
        }
        if (IsKeyDown(KEY_FOUR)) {
            fluid_handler.setDragObjectActive(false);
            fluid_handler.setForceObjectActive(false);
            fluid_handler.setGeneratorActive(false);
            obstacle_handler.setSolidDrawer(true);
            scene_renderer.setZoomObjectActive(false);
        }
        if (IsKeyDown(KEY_FIVE)) {
            fluid_handler.setDragObjectActive(false);
            fluid_handler.setForceObjectActive(false);
            fluid_handler.setGeneratorActive(false);
            obstacle_handler.setSolidDrawer(false);
            scene_renderer.setZoomObjectActive(true);
        }
        if (IsKeyDown(KEY_G)) {
            fluid_attributes.addToGravityY(100);
        }
        if (IsKeyDown(KEY_N)) {
            fluid_attributes.addToGravityY(-100);
        }
        if (IsKeyDown(KEY_M)) {
            fluid_attributes.addToGravityX(100);
        }
        if (IsKeyDown(KEY_H)) {
            fluid_attributes.addToGravityX(-100);
        }
        if (IsKeyDown(KEY_C)) {
            if (fluid_handler.pressure_solver.getDivergenceModifier() > 0) {
                fluid_handler.pressure_solver.addToDivergenceModifier(-1);
            }
        }
        if (IsKeyDown(KEY_D)) {
            fluid_handler.pressure_solver.addToDivergenceModifier(1);
        }
        if (IsKeyPressed(KEY_A)) {
            scene_renderer.fluid_renderer.setNextRenderPattern();
        }
        if (IsKeyPressed(KEY_F)) {
            fluid_attributes.setFireActive(!fluid_attributes.getFireActive());
        }
        if (IsKeyPressed(KEY_Y)) {
            fluid_attributes.setStop(!fluid_attributes.getStop());
        }
        if (IsKeyPressed(KEY_U)) {
            fluid_attributes.setStep(true);
        }
        if (IsKeyPressed(KEY_Q)) {
            std::cout << 
            "Fill Grid: " << getFillGridTime() << "\n" <<
            "Miscellaneous: " << getMiscellaneousTime() << "\n" <<
            "Collision: " << getCollisionTime() << "\n" <<
            "Obstacle Collision: " << getObstacleCollisionTime() << "\n" <<
            "To Grid: " << getToGridTime() << "\n" <<
            "Density Update: " << getDensityUpdateTime() << "\n" <<
            "Projection: " << getProjectionTime() << "\n" <<
            "To Particles: " << getToParticlesTime() << "\n" <<
            "Rendering: " << getRenderingTime() << "\n"; // <<
            //"Whole Step: " << fluid_handler.getSimStepTime() << "\n" <<
            //"Combined: " << fluid_handler.getCombinedTime() << "\n" <<
            //"Before Sim Step: " << beforeSimStep << "\n" <<
            //"After Sim Step: " << afterSimStep << "\n";
            CloseWindow();
        }
        /*end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        addValueToAverage(afterSimStep, duration.count(), numDT);*/
    }

    void track_mouse_events() {

        fluid_attributes.frame_context.leftMouseDown = IsMouseButtonDown(MOUSE_BUTTON_LEFT);
        fluid_attributes.frame_context.rightMouseDown = IsMouseButtonDown(MOUSE_BUTTON_RIGHT);
        fluid_attributes.frame_context.justPressed = IsMouseButtonPressed(MOUSE_BUTTON_LEFT) || IsMouseButtonPressed(MOUSE_BUTTON_RIGHT);
        
        float mouseWheelDelta = GetMouseWheelMove();

        if (IsMouseButtonPressed(MOUSE_BUTTON_LEFT) || IsMouseButtonPressed(MOUSE_BUTTON_RIGHT)) {
            if (!showControls && PointIsInsideRectangle(fluid_attributes.frame_context.screen_mouse_pos.x, fluid_attributes.frame_context.screen_mouse_pos.y, showControlsX, showControlsY, showControlsWIDTH, showControlsHEIGHT)) {
                showControls = true;
            }
            else if (showControls && PointIsInsideRectangle(fluid_attributes.frame_context.screen_mouse_pos.x, fluid_attributes.frame_context.screen_mouse_pos.y, hideControlsX, hideControlsY, hideControlsWIDTH, hideControlsHEIGHT)) {
                showControls = false;
            }
        }

        if (mouseWheelDelta) {
            if (fluid_handler.getObjectsActive()) {
                fluid_handler.addToObjectRadius(20 * mouseWheelDelta);
            }
            
            else if (obstacle_handler.getPencilActive()) {
                int pencilRadius = obstacle_handler.getPencilRadius();
                bool limit_size = (pencilRadius > fluid_attributes.getNumX() / 2 && mouseWheelDelta > 0);
                if (!limit_size) {
                    if (pencilRadius + mouseWheelDelta < 0) {
                        obstacle_handler.setPencilRadius(0);
                    }
                    else {
                        obstacle_handler.addToPencilRadius(mouseWheelDelta);
                    }
                }
            }

            else if (scene_renderer.getZoomObjectActive()) {
                scene_renderer.wheelZoom(mouseWheelDelta);
            }
        }
    }

    int DrawTextNearRightWall(std::string &text, float yPos) {
        int strWIDTH = MeasureText(text.c_str(), fontSize);
        int strHEIGHT = MeasureText(text.c_str(), fontSize);
        DrawText(text.c_str(), fluid_attributes.frame_context.WIDTH - strWIDTH - 2 * fluid_attributes.cellSpacing, yPos, fontSize, WHITE);
        return strHEIGHT;
    }

    void displayGUI() {
        // reset fps tracker every 20 frames
        if (frame == 0) {
            fps = 1 / trueDT;
            frame = 20;
        }
        frame--;

        float WIDTH = fluid_attributes.frame_context.WIDTH;
        float HEIGHT = fluid_attributes.frame_context.HEIGHT;

        int fpsX = GetFractionOf(WIDTH, 0.0077f);
        int fpsY = GetFractionOf(HEIGHT, 0.0087f);

        std::string fps_str = std::to_string(fps);
        if (!fluid_attributes.getStop()) {
            DrawText(fps_str.c_str(), fpsX, fpsY, fontSize, WHITE);
        }
        else {
            DrawText(fps_str.c_str(), fpsX, fpsY, fontSize, RED);
        }

        // controls
        if (showControls) {
            std::string flip_ratio_str = "PIC/FLIP Ratio (S/B):    " + std::to_string(fluid_attributes.getFlipRatio());
            DrawTextNearRightWall(flip_ratio_str, showControlsHEIGHT);

            std::string vorticity_confinement_str = "Vorticity Confinement (W/E):    " + std::to_string(fluid_attributes.getVorticityStrength());
            DrawTextNearRightWall(vorticity_confinement_str, showControlsHEIGHT * 2.5);

            std::string gravityX_str = "Gravity X (H/M):    " + std::to_string(static_cast<int>(fluid_attributes.getGravityX()));
            DrawTextNearRightWall(gravityX_str, showControlsHEIGHT * 4);

            std::string gravityY_str = "Gravity Y (N/G):    " + std::to_string(static_cast<int>(fluid_attributes.getGravityY()));
            DrawTextNearRightWall(gravityY_str, showControlsHEIGHT * 5.5);

            std::string render_str = "Switch Rendering Modes (A)";
            DrawTextNearRightWall(render_str, showControlsHEIGHT * 7.5);

            std::string fire_str = "Toggle Fire During Temperature Mode (F)";
            DrawTextNearRightWall(fire_str, showControlsHEIGHT * 9);

            std::string objects_str = "Switch Object (1-5), Use Mouse Wheel To Resize";
            DrawTextNearRightWall(objects_str, showControlsHEIGHT * 10.5);

            std::string num_particles_str = "Number of Particles:     " + std::to_string(fluid_attributes.num_particles);
            DrawTextNearRightWall(num_particles_str, showControlsHEIGHT * 12);

            DrawTextNearRightWall(hideControlsStr, hideControlsY);

            /*DrawCircle(hideControlsX, hideControlsY, 10, RED);
            DrawCircle(hideControlsX + hideControlsWIDTH, hideControlsY + hideControlsHEIGHT, 10, RED);*/
        }
        else { 
            DrawTextNearRightWall(showControlsStr, showControlsY);

            /*DrawCircle(showControlsX, showControlsY, 10, RED);
            DrawCircle(showControlsX + showControlsWIDTH, showControlsY + showControlsHEIGHT, 10, RED);*/
        }

    }

    void handle_zoom() {
        if (fluid_attributes.frame_context.leftMouseDown && scene_renderer.getZoomObjectActive()) {
            scene_renderer.dragCamera();
        }
        if (fluid_attributes.frame_context.zooming) {
            float zoom = fluid_attributes.frame_context.zoom_amount;
            fluid_handler.objectSimRadius = fluid_handler.objectRenderRadius / zoom;
        }
    }

    void render_interaction_objects() {
        fluid_handler.render_objects();
        obstacle_handler.render_objects();
    }

    float getCombinedTime() {
        return FillGridTime + miscellaneousTime + CollisionTime + ObstacleCollisionTime + ToGridTime + DensityUpdateTime + ProjectionTime + ToParticlesTime + RenderingTime;
    }

    float getSimStepTime() {
        return SimStepTime;
    }

    float getFillGridTime() {
        return FillGridTime;
    }

    float getMiscellaneousTime() {
        return miscellaneousTime;
    }

    float getCollisionTime() {
        return CollisionTime;
    }

    float getObstacleCollisionTime() {
        return ObstacleCollisionTime;
    }

    float getToGridTime() {
        return ToGridTime;
    }

    float getDensityUpdateTime() {
        return DensityUpdateTime;
    }

    float getProjectionTime() {
        return ProjectionTime;
    }

    float getToParticlesTime() {
        return ToParticlesTime;
    }

    float getRenderingTime() {
        return RenderingTime;
    }

};