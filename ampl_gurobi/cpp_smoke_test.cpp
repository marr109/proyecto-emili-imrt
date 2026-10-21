#include <ampl/ampl.h>
#include <cstdlib>
#include <iostream>

int main()
{
    ampl::Environment env(AMPL_BIN_DIR);
    ampl::AMPL ampl(env);

    ampl.eval("var x >= 0; maximize obj: x; subject to c: x <= 10;");
    // Path to the gurobi AMPL driver, e.g. <venv>/lib/python3.X/site-packages/ampl_module_gurobi/bin/gurobi
    const char* solver = std::getenv("EMILI_GUROBI_BIN");
    if (!solver) { std::cerr << "set EMILI_GUROBI_BIN to the gurobi AMPL driver" << std::endl; return 1; }
    ampl.setOption("solver", solver);
    ampl.solve();

    std::cout << "solve_result: " << ampl.getValue("solve_result").str() << std::endl;
    std::cout << "x = " << ampl.getVariable("x").value() << std::endl;
    return 0;
}
