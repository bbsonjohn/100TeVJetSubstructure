#include <string.h>
#include <map>
#include "processes.h"

using namespace std;

Processes::Processes(void)
{
	Catalog["ttB"] = ttB;
	Catalog["tt"] = tt;
	Catalog["tB"] = tB;
	Catalog["Stops"] = Stops;

    ttB["0-1500"] = 206.01;
    ttB["1500-3000"] = 12.58;
    ttB["3000-5500"] = 1.18;
    ttB["5500-9000"] = 0.092;
    ttB["9000-100000"] = 0.009;

    tt["0-1000"] = 29141.30;
    tt["1000-2000"] = 1777.28;
    tt["2000-3500"] = 185.22;
    tt["3500-5500"] = 18.92;
    tt["5500-8500"] = 2.39;
    tt["8500-100000"] = 0.277;

    tB["0-1000"] = 3399.65;
    tB["1000-2000"] = 165.25;
    tB["2000-3500"] = 15.57;
    tB["3500-6000"] = 1.59;
    tB["6000-9000"] = 0.11;
    tB["9000-100000"] = 0.013;

    Stops["1500"] = 0.535;
    Stops["2000"] = 0.123;
    Stops["2500"] = 0.038;
    Stops["3000"] = 0.0142;
    Stops["3500"] = 0.00593;
    Stops["4000"] = 0.00276;
    Stops["4500"] = 0.00136;
    Stops["5000"] = 0.00071;
    Stops["5500"] = 0.00038;
    Stops["6000"] = 0.00022;
    Stops["6500"] = 0.00013;
    Stops["7000"] = 0.000077;
    Stops["7500"] = 0.000047;
    Stops["8000"] = 0.000029;
    Stops["8500"] = 0.000019;
    Stops["9000"] = 0.000012;
}