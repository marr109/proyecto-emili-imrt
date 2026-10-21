
//  Created by Federico Pagnozzi on 28/11/14.
//  Copyright (c) 2014 Federico Pagnozzi. All rights reserved.
//  This file is distributed under the BSD 2-Clause License. See LICENSE.TXT
//  for details.

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cstring>
#include <algorithm>
#include "generalParser.h"
#include "setup.h"
#include "imrt/imrt_builder.h"
#include "imrt/imrt.h"
#include "imrt/imrt_bao.h"
#include "imrt/imrt_report.h"


/*---------------------------------------------------------------------------*
 * Final plan report
 *
 * Registered with atexit() so it also runs when the search is cut short by the
 * wall-clock budget (-it): emili::finalise() handles SIGALRM and calls exit(0)
 * directly, which skips the tail of main() where this used to live. The guard
 * makes the normal path and the atexit path idempotent -- whichever runs first
 * emits, the other is a no-op.
 *---------------------------------------------------------------------------*/
static emili::LocalSearch* g_ls = nullptr;
static bool g_report_emitted = false;

static void emitPlanReport()
{
    if (g_report_emitted || g_ls == nullptr) return;
    g_report_emitted = true;

    emili::Solution* solution = g_ls->getBestSoFar();
    if (solution == nullptr) return;
    double solval = solution->getSolutionValue();
    emili::LocalSearch* ls = g_ls;

    // ── Clinical-style plan report (ICRU-83 metrics) ──────────────────
    // Re-implemented on top of IFmoDataSource (see imrt/imrt_report.h)
    // after the old ImrtInstance-based reportPlan was removed with the
    // CORT-format loader.
    emili::Problem* prob = &ls->getInitialSolution().getProblem();
    if (auto* baoProb = dynamic_cast<emili::imrt::BaoProblem*>(prob)) {
        auto* bs = dynamic_cast<emili::imrt::BaoSolution*>(solution);
        if (bs) {
            emili::imrt::reportPlan(baoProb->getSource(), bs->intensities_,
                                     solval, bs->angle_degrees_, std::cout);
        }
    } else if (auto* imrtProb = dynamic_cast<emili::imrt::ImrtProblem*>(prob)) {
        auto* is = dynamic_cast<emili::imrt::ImrtSolution*>(solution);
        if (is) {
            // getActiveAngles() empty means "all angles active" (see
            // ImrtProblem::isAngleActive); map the active indices (or
            // the full catalog) to their degree values for the header.
            const std::vector<int>& active = imrtProb->getActiveAngles();
            const std::vector<int>& catalog = imrtProb->getAngleDegrees();
            std::vector<int> deg;
            if (active.empty()) {
                deg = catalog;
            } else {
                deg.reserve(active.size());
                for (int idx : active) deg.push_back(catalog[idx]);
            }
            emili::imrt::reportPlan(imrtProb->getSource(), is->getIntensities(),
                                     solval, deg, std::cout);
        }
    }

}

int main(int argc, char *argv[])
{
    prs::emili_header();
    srand ( time(0) );
    clock_t time = clock();
    if (argc < 3 )
    {
        prs::info();
        return 1;
    }

    float pls = 0;
    emili::LocalSearch* ls;

    prs::GeneralParserE ps(argv, argc);
    prs::EmBaseBuilder emb(ps, ps.getTokenManager());
    prs::imrt::ImrtBuilder imrtb(ps, ps.getTokenManager());
    ps.addBuilder(&emb);
    ps.addBuilder(&imrtb);

    ls = ps.parseParams();
    g_ls = ls;
    atexit(emitPlanReport);

    if(ls != nullptr)
    {
        pls = ls->getSearchTime();
        emili::Solution* solution;
        std::cout << std::endl << "Searching..." << std::endl;
        if(pls > 0)
        {
            solution = ls->timedSearch(pls);
        }
        else
        {
            solution = ls->search();
        }
        if(!emili::get_print())
        {
            solution = ls->getBestSoFar();
            double time_elapsed = (double)(clock()-time)/CLOCKS_PER_SEC;
            double solval = solution->getSolutionValue();
            std::cout << "time : " << time_elapsed << std::endl;
            std::cout << "iteration counter : " << emili::iteration_counter() << std::endl;
            std::cerr << solution->getSolutionValue() << std::endl;
            std::cout << "Objective function value: " << std::fixed << solval << std::endl;
            std::cerr << std::fixed << solval << std::endl;
            std::cout << "Found solution: ";
            std::cout << solution->getSolutionRepresentation() << std::endl;
            std::cout << std::endl;

            emitPlanReport();
        }
        g_ls = nullptr;
        delete ls;
    }
}
