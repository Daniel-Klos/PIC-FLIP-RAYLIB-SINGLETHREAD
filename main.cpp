#include <raylib.h>
#include <iostream>

#include "scene_handler.hpp"

int main() 
{

    // SETTINGS
    // ---------------------------------------
    // numparticles and gravity are self explanatory
    // gridNumX is how large the horizontal area is
    // flipRatio is the blend between PIC and FLIP, expects a float from [0.f, 1.f].
    // vorticitystrength is how strong vorticity confinement forces are if you choose to include that in the sim 
    // maxFps sets the deltaTime for the simulation and the framerate for the window

    int numParticles = 30000; // 30000 --- start with a kinda large number so that good density sample is taken at the start of the simulation
    float gravityY = 5500.f; // 5500
    float gravityX = 0.f;
    int gridNumX = 250; // 250
    float flipRatio = 0.9f; // 0.9f
    float vorticityStrength = 0.f;
    int maxFps = 120;

    SetTraceLogLevel(LOG_NONE);

    int WIDTH = 2000;
    int HEIGHT = 1300;

    InitWindow(WIDTH, HEIGHT, "Flip Fluid");

    FluidState fluid_attributes = FluidState(numParticles, gridNumX, vorticityStrength, flipRatio, gravityX, gravityY, WIDTH, HEIGHT);

    SceneHandler scene_handler = SceneHandler(fluid_attributes, maxFps);

    while (!WindowShouldClose())
    {
        scene_handler.simulate();
    }

    return 0;
}