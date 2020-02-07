#ifndef ROOT_SOIL_INTERACTION_H_
#define ROOT_SOIL_INTERACTION_H_

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>
#include <dune/pybindxi/numpy.h>
namespace py = pybind11;


#include <memory>
#include <MappedOrganism.h>

/**
 * Describes the exchange of root and soil via sink terms
 */
class ROSI_HybridRichards
{
public:

    ROSI_HybridRichards(std::shared_ptr<CPlantBox::MappedRootSystem> rs) : rs(rs) { }
    virtual ~ROSI_HybridRichards() { };

    /**
     * Pressure of soil surrounding the root segment with index sIdx [cm]
     * (needed by the root part)
     */
    virtual double soilPressure(int sIdx) {
        int cellIdx = rs->seg2cell.at(sIdx);
        return oldSoilX[cellIdx];
    }

    /**
     * Fluxes from root segments into a the soil cell with cell index cIdx [TODO]
     * (needed by the soil part)
     */
    virtual double roots2cell(int cIdx, double soilP)
    {
        if (rs->cell2seg.count(cIdx)>0) {
            auto sIdxs = rs->cell2seg.at(cIdx);
            double flux = 0;
            for (int i : sIdxs) {
                double f = 0.;
                if (i < oldRootX.size()) {
                    double rootP = oldRootX[i];
                    double a = rs->radii[i];
                    auto n1 = rs->nodes[rs->segments[i].x];
                    auto n2 = rs->nodes[rs->segments[i].y];
                    double l = (n2.minus(n1)).length();

                    // f = 2*

                }
                flux += f;
            }

            return flux;
        } else {
            return 0.;
        }
    }

    /**
     * Root radial conductivity [TODO]
     * needed by the soil part
     */
    virtual double kr(int sIdx) {

        double age = rs->getSimTime() - rs->nodeCTs.at(sIdx+1);
        int type = rs->types.at(sIdx);
        return 0.; // TODO

    }

    /**
     * Root asxial conductivity [TODO]
     * needed by both parts
     */
    virtual double kz(int sIdx) {

        double age = rs->getSimTime() - rs->nodeCTs.at(sIdx+1);
        int type = rs->types.at(sIdx);
        return 0.; // TODO

    }

    void setKr() {

    }

    void setKz() {

    }

    // solutions of old time step (TODO) smart pointer?
    std::vector<double> oldSoilX;
    std::vector<double> oldRootX;

    //
    std::vector<double> krAge;
    std::vector<double> krValue;
    std::vector<double> kzAge;
    std::vector<double> kzValue;

    std::shared_ptr<CPlantBox::MappedRootSystem> rs;

};


#endif
