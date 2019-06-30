// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup EmbeddedCoupling
 * \ingroup CCModel
 * \brief Coupling manager for low-dimensional domains embedded in the bulk
 *        domain. Intersection computation relies on the dune grid-glue backend.
 *        Both domain are assumed to be discretized using a cc finite volume scheme.
 */

#ifndef DUMUX_CC_GRIDGLUE_EMBEDDEDCOUPLINGMANAGER_CIRCLEAVERAGE_HH
#define DUMUX_CC_GRIDGLUE_EMBEDDEDCOUPLINGMANAGER_CIRCLEAVERAGE_HH

#include <dumux/multidimension/embeddedcoupling/cellcentered/gridgluecouplingmanager.hh>
#include <dumux/multidimension/embeddedcoupling/circlegridfactory.hh>

namespace Dumux
{

/*!
 * \ingroup EmbeddedCoupling
 * \ingroup CCModel
 * \brief Manages the coupling between bulk elements and lower dimensional elements
 * The pressure for the source term is computed as an integral over the circle describing
 * the surface of the one-dimensional tube. This exact determination of the bulk pressure
 * is necessary for convergence studies.
 */
template<class TypeTag>
class CCGridGlueEmbeddedCouplingManagerCircleAverage : public CCGridGlueEmbeddedCouplingManager<TypeTag>
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using ParentType = Dumux::CCGridGlueEmbeddedCouplingManager<TypeTag>;

    // obtain the type tags of the sub problems
    using BulkProblemTypeTag = typename GET_PROP_TYPE(TypeTag, BulkProblemTypeTag);
    using LowDimProblemTypeTag = typename GET_PROP_TYPE(TypeTag, LowDimProblemTypeTag);

    using BulkGridView = typename GET_PROP_TYPE(BulkProblemTypeTag, GridView);
    using LowDimGridView = typename GET_PROP_TYPE(LowDimProblemTypeTag, GridView);

    using BulkProblem = typename GET_PROP_TYPE(BulkProblemTypeTag, Problem);
    using LowDimProblem = typename GET_PROP_TYPE(LowDimProblemTypeTag, Problem);

    using BulkGrid = typename GET_PROP_TYPE(BulkProblemTypeTag, Grid);
    using LowDimGrid = typename GET_PROP_TYPE(LowDimProblemTypeTag, Grid);

    using BulkPrimaryVariables = typename GET_PROP_TYPE(BulkProblemTypeTag, PrimaryVariables);
    using LowDimPrimaryVariables = typename GET_PROP_TYPE(LowDimProblemTypeTag, PrimaryVariables);

    enum {
        bulkDim = BulkGridView::dimension,
        lowDimDim = LowDimGridView::dimension,
        dimWorld = BulkGridView::dimensionworld
    };

    using PointSourceData = Dumux::PointSourceDataCircleAverage<TypeTag>;

    using BulkElement = typename BulkGridView::template Codim<0>::Entity;
    using LowDimElement = typename LowDimGridView::template Codim<0>::Entity;
    using BulkExtractor = typename Dune::GridGlue::Codim0Extractor<BulkGridView>;
    using LowDimExtractor = typename Dune::GridGlue::Codim0Extractor<LowDimGridView>;
    using GlueType = typename Dune::GridGlue::GridGlue<BulkExtractor, LowDimExtractor>;

public:

    /*!
     * \brief Constructor
     */
    CCGridGlueEmbeddedCouplingManagerCircleAverage(BulkProblem& bulkProblem, LowDimProblem& lowDimProblem)
    : ParentType(bulkProblem, lowDimProblem) {}

    /* \brief Compute integration point point sources and associated data
     *
     * This method uses grid glue to intersect the given grids. Over each intersection
     * we later need to integrate a source term. This method places point sources
     * at each quadrature point and provides the point source with the necessary
     * information to compute integrals (quadrature weight and integration element)
     * \param order The order of the quadrature rule for integration of sources over an intersection
     * \param verbose If the point source computation is verbose
     */
    void computePointSourceData(unsigned int order = 1, bool verbose = false)
    {
        // Initialize the descriptors (all elements)
        auto bulkDescriptor = [](const BulkElement& element, unsigned int subentity) { return true; };
        auto lowDimDescriptor = [](const LowDimElement& element, unsigned int subentity) { return true; };
        // construct element extractors
        BulkExtractor bulkExtractor(this->bulkGridView(), bulkDescriptor);
        LowDimExtractor lowDimExtractor(this->lowDimGridView(), lowDimDescriptor);
        // choose overlapping merger
        Dune::GridGlue::OverlappingMerge<bulkDim, lowDimDim, dimWorld> merger;
        // build the glue object
        GlueType glue(bulkExtractor, lowDimExtractor, &merger);

        // enable fallback brute-force step if the advancing front algorithm didn't find a seed
        merger.enableFallback(true);

        // compute intersections
        glue.build();

        // write output for debug purposes
        if (verbose)
        {
            GridGlueVtkWriter glueVtk;
            glueVtk.write(glue, "glue", verbose);
        }

        // initilize the maps
        // do some logging and profiling
        Dune::Timer watch;
        std::cout << "Initializing the point sources..." << std::endl;

        // clear all internal members like pointsource vectors and stencils
        // initializes the point source id counter
        clear();

        // iterate over all intersections
        // we choose the target side (1, the low dimensional domain)
        const auto isEndIt = glue.template iend<1>();
        for (auto isIt = glue.template ibegin<1>(); isIt != isEndIt; ++isIt)
        {
            // all inside elements are identical...
            const auto& inside = isIt->inside(0);
            // get the intersection geometry for integrating over it
            const auto intersectionGeometry = isIt->geometry();

            // get the Gaussian quadrature rule for the local intersection
            using Quad = Dune::QuadratureRule<Scalar, lowDimDim>;
            const Quad &quad = Dune::QuadratureRules<Scalar, lowDimDim>::rule(intersectionGeometry.type(), order);

            const unsigned int lowDimElementIdx = this->lowDimProblem().elementMapper().index(inside);

            // iterate over all quadrature points
            for (auto&& qp : quad)
            {
                // get global position of quadrature point
                const auto globalPos = intersectionGeometry.global(qp.position());

                //////////////////////////////////////////////////////////
                // get circle average connectivity and interpolation data
                //////////////////////////////////////////////////////////

                // make circle grid
                auto numSegments = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, MultiDimension, NumCircleSegments);
                const auto radius = this->lowDimProblem().spatialParams().radius(lowDimElementIdx);
                const auto normal = intersectionGeometry.corner(1)-intersectionGeometry.corner(0);
                std::shared_ptr<LowDimGrid> circleGrid(Dumux::CircleGridFactory<LowDimGrid>::createGrid(globalPos, normal, radius, numSegments));

                // Debug output
                // Dune::VTKWriter<LowDimGridView> writer(circleGrid->leafGridView());
                // writer.write("circle_" + std::to_string(circleCounter_++));

                // redirect grid-glue output to file because it would spam the console
                std::ofstream out("grid-glue-circle.log");
                auto coutbuf = std::cout.rdbuf(out.rdbuf()); //save and redirect

                // only intersect with bulk elements within range of the integration point
                auto bulkDescriptor = [&globalPos, &radius](const BulkElement& element, unsigned int subentity)
                {
                    const auto geometry = element.geometry();
                    const auto dist = (geometry.center()-globalPos).two_norm();
                    const auto diag = 2*(geometry.center()-geometry.corner(0)).two_norm();

                    return dist < std::max(2*radius, diag) && dist > std::min(radius/2, diag);
                };
                auto lowDimDescriptor = [](const LowDimElement& element, unsigned int subentity) { return true; };
                LowDimExtractor circleExtractor(circleGrid->leafGridView(), lowDimDescriptor);
                BulkExtractor bulkExtractor(this->bulkGridView(), bulkDescriptor);
                Dune::GridGlue::OverlappingMerge<bulkDim, lowDimDim, dimWorld> circleMerger;
                GlueType circleGlue(bulkExtractor, circleExtractor, &circleMerger);

                // calculate circle intersections with bulk
                circleGlue.build();

                //reset to standard output again
                std::cout.rdbuf(coutbuf);

                // loop over local intersections
                std::vector<Scalar> circleIpWeight;
                std::vector<unsigned int> circleStencil;
                const auto circleIsEndIt = circleGlue.template iend<1>();
                for (auto circleIsIt = circleGlue.template ibegin<1>(); circleIsIt != circleIsEndIt; ++circleIsIt)
                {
                    // use trapezoidal rule
                    const auto cirlceIntersectionGeometry = circleIsIt->geometry();
                    const auto length = cirlceIntersectionGeometry.volume();

                    circleStencil.push_back(this->bulkProblem().elementMapper().index(circleIsIt->outside(0)));
                    circleIpWeight.push_back(length);
                }

                // export low dim circle stencil
                lowDimCircleStencils_[lowDimElementIdx].insert(lowDimCircleStencils_[lowDimElementIdx].end(),
                                                               circleStencil.begin(),
                                                               circleStencil.end());

                /////////////////////////////////////////////////////////////////////////////////////////////

                // compute the coupling point sources
                for (unsigned int outsideIdx = 0; outsideIdx < isIt->neighbor(0); ++outsideIdx)
                {
                    const auto& outside = isIt->outside(outsideIdx);
                    const unsigned int bulkElementIdx = this->bulkProblem().elementMapper().index(outside);

                    // each quadrature point will be a point source for the sub problem
                    const auto id = idCounter_++;
                    const auto qpweight = qp.weight();
                    const auto ie = intersectionGeometry.integrationElement(qp.position());
                    this->bulkPointSources().emplace_back(globalPos, id, qpweight, ie, std::vector<unsigned int>({bulkElementIdx}));
                    this->bulkPointSources().back().setEmbeddings(isIt->neighbor(0));
                    this->lowDimPointSources().emplace_back(globalPos, id, qpweight, ie, std::vector<unsigned int>({lowDimElementIdx}));
                    this->lowDimPointSources().back().setEmbeddings(isIt->neighbor(0));

                    // pre compute additional data used for the evaluation of
                    // the actual solution dependent source term
                    PointSourceData psData;
                    psData.addLowDimInterpolation(lowDimElementIdx);
                    psData.addBulkInterpolation(bulkElementIdx);
                    // add data needed to compute integral over the circle
                    psData.addCircleInterpolation(circleIpWeight, circleStencil);

                    // publish point source data in the global vector
                    pointSourceData_.push_back(psData);

                    // compute the coupling stencils
                    this->bulkCouplingStencils()[bulkElementIdx].push_back(lowDimElementIdx);

                    // export bulk circle stencil
                    bulkCircleStencils_[bulkElementIdx].insert(bulkCircleStencils_[bulkElementIdx].end(),
                                                               circleStencil.begin(),
                                                               circleStencil.end());
                }
            }
        }

        // make the circle stencils unique
        for (auto&& stencil : bulkCircleStencils_)
        {
            std::sort(stencil.second.begin(), stencil.second.end());
            stencil.second.erase(std::unique(stencil.second.begin(), stencil.second.end()), stencil.second.end());
        }

        for (auto&& stencil : lowDimCircleStencils_)
        {
            std::sort(stencil.second.begin(), stencil.second.end());
            stencil.second.erase(std::unique(stencil.second.begin(), stencil.second.end()), stencil.second.end());
        }

        // coupling stencils
        for (auto&& stencil : this->bulkCouplingStencils())
        {
            std::sort(stencil.second.begin(), stencil.second.end());
            stencil.second.erase(std::unique(stencil.second.begin(), stencil.second.end()), stencil.second.end());
        }

        // the low dim coupling stencil
        for (auto&& stencil : this->lowDimCouplingStencils())
        {
            // add the circle stencil to the coupling stencil
            stencil.second.insert(stencil.second.end(),
                                  lowDimCircleStencils_.at(stencil.first).begin(),
                                  lowDimCircleStencils_.at(stencil.first).end());
            std::sort(stencil.second.begin(), stencil.second.end());
            stencil.second.erase(std::unique(stencil.second.begin(), stencil.second.end()), stencil.second.end());
        }

        // compute the bulk stencil (like in the fv element geometry)
        for (const auto& element : elements(this->bulkGridView()))
        {
            // the element itself is first
            const unsigned int eIdx = this->bulkProblem().elementMapper().index(element);
            this->bulkStencils().insert({eIdx, {eIdx}});

            for (const auto& intersection : intersections(this->bulkGridView(), element))
                if (intersection.neighbor())
                    this->bulkStencils().at(eIdx).push_back(this->bulkProblem().elementMapper().index(intersection.outside()));
        }

        // ...and afterwards add the bulk circle stencil so that the order doesn't get messed up.
        for (auto&& stencil : bulkCircleStencils_)
        {
            auto&& bulkStencil = this->bulkStencils().at(stencil.first);
            // only insert elements that are not already in the stencil
            for (auto&& eIdx : stencil.second)
                if (std::find(bulkStencil.begin(), bulkStencil.end(), eIdx) == bulkStencil.end())
                    bulkStencil.push_back(eIdx);
        }

        // compute the low dimensional stencil (like in the fv geometry)
        for (const auto& element : elements(this->lowDimGridView()))
        {
            // the element itself is first
            const unsigned int eIdx = this->lowDimProblem().elementMapper().index(element);
            this->lowDimStencils().insert({eIdx, {eIdx}});

            for (const auto& intersection : intersections(this->lowDimGridView(), element))
                if (intersection.neighbor())
                    this->lowDimStencils().at(eIdx).push_back(this->lowDimProblem().elementMapper().index(intersection.outside()));
        }

        std::cout << "took " << watch.elapsed() << " seconds." << std::endl;
    }

    /*!
     * \brief Methods to be accessed by the subproblems
     */
    // \{

    //! Return a reference to the pointSource data
    const PointSourceData& pointSourceData(unsigned int id) const
    {
        return pointSourceData_[id];
    }

    //! Return data for a bulk point source with the identifier id
    const BulkPrimaryVariables& bulkPriVars(unsigned int id) const
    {
        auto& data = pointSourceData_[id];
        data.interpolateBulk(this->bulkProblem().model().curSol());
        return data.bulkPriVars();
    }

    //! Return data for a low dim point source with the identifier id
    const LowDimPrimaryVariables& lowDimPriVars(unsigned int id) const
    {
        auto& data = pointSourceData_[id];
        data.interpolateLowDim(this->lowDimProblem().model().curSol());
        return data.lowDimPriVars();
    }

    //! Return a reference to the bulk problem
    Scalar radius(unsigned int id) const
    {
        const auto& data = pointSourceData_[id];
        return this->lowDimProblem().spatialParams().radius(data.lowDimElementIdx());
    }
    // \}

    //! Clear all internal data members
    void clear()
    {
        ParentType::clear();
        pointSourceData_.clear();
        bulkCircleStencils_.clear();
        lowDimCircleStencils_.clear();
        idCounter_ = 0;
    }

private:
    // id generator for point sources
    unsigned int idCounter_;

    // special point source data for the circle average
    mutable std::vector<PointSourceData> pointSourceData_;

    // circle stencils
    std::unordered_map<unsigned int, std::vector<unsigned int> > bulkCircleStencils_;
    std::unordered_map<unsigned int, std::vector<unsigned int> > lowDimCircleStencils_;
};

} // namespace Dumux

#endif
