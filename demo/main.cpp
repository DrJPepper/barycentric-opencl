#include "mesh.h"
#include <toml++/toml.hpp>

int main(int argc, char *argv[]) {
    if (argc < 2) {
        cout << "ERROR: Input file name required\n";
        exit(1);
    }
    string inputFileName = argv[1];
    std::ifstream inFile(inputFileName);
    if (!inFile.good()) {
        cerr << "ERROR: input file invalid or doesn't exist\n";
        return 1;
    }

    toml::table setts;
    try {
        setts = toml::parse_file(inputFileName);
    } catch (const toml::parse_error& err) {
        cerr << "Parsing TOML file failed:\n" << err << "\n";
        return 1;
    }

    string modelFile = setts["model"]["model"].value_or("../obj_models/model.obj");
    string inputType = setts["model"]["type"].value_or("obj");

    Mesh mesh = Mesh(modelFile, inputType);
    mesh.scale = setts["model"]["scale"].value_or(1.0);

    if (!inputType.compare("obj")) {
        auto pc = setts["model"]["param_corners"];
        mesh.paramCorners << pc[0].value_or(0), pc[1].value_or(0),
                     pc[2].value_or(0), pc[3].value_or(0);
    } else if (!inputType.compare("raytrace")) {
        auto rb1 = setts["ray_trace"]["ray_start"];
        auto rb2 = setts["ray_trace"]["ray_box_size"];
        mesh.rayBox << rb1[0].value_or(0.0), rb1[1].value_or(0.0), rb2[0].value_or(10.0), rb2[1].value_or(10.0);
        mesh.rtsettings->scale = setts["ray_trace"]["scale"].value_or(1.0);
        auto rot = setts["ray_trace"]["rotate"];
        mesh.rtsettings->rotate << rot[0].value_or(0), rot[1].value_or(0), rot[2].value_or(0);
        mesh.rtsettings->x = setts["ray_trace"]["x"].value_or(20);
        mesh.rtsettings->y = setts["ray_trace"]["y"].value_or(20);
        auto Zv = setts["ray_trace"]["Zv"];
        mesh.rtsettings->Zv << Zv[0].value_or(1), Zv[1].value_or(0), Zv[2].value_or(0);
        auto Vup = setts["ray_trace"]["Vup"];
        mesh.rtsettings->Vup << Vup[0].value_or(0), Vup[1].value_or(1), Vup[2].value_or(0);
    }

    mesh.initializeMesh();
    cout << mesh.barycentricLoopingF(Vector2d{0.5, 0.5}) << endl;
    
    return 0;
}
