#include "../imrt/cerr_instance.h"
#include "../imrt/imrt_fmo.h"
#include <iostream>

int main()
{
    emili::imrt::CerrFmoSource src("instances/CERR_Prostate");
    emili::imrt::ImrtFmoSolver solver(src);
    if (!solver.isReady()) {
        std::cerr << "solver not ready\n";
        return 1;
    }
    auto res = solver.solve({0, 90, 180, 270});
    std::cout.precision(6);
    std::cout << "objective: " << std::fixed << res.objective << std::endl;
    return 0;
}
