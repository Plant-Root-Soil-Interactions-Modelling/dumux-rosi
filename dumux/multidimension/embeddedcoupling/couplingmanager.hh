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
 * \brief Coupling manager for low-dimensional domains embedded in the bulk
          domain. Intersection computation relies on the dune grid-glue backend.
 */

#ifndef DUMUX_EMBEDDEDCOUPLINGMANAGER_HH
#define DUMUX_EMBEDDEDCOUPLINGMANAGER_HH

#include <iostream>
#include <fstream>
#include <string>
#include <utility>

#include <dune/common/timer.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/localfunctions/lagrange/pqkfactory.hh>

#include <dune/grid-glue/gridglue.hh>
#include <dune/grid-glue/extractors/extractorpredicate.hh>
#include <dune/grid-glue/extractors/codim0extractor.hh>
#include <dune/grid-glue/merging/overlappingmerge.hh>

#include <dumux/io/gridgluevtkwriter.hh>
//#include <dumux/io/pointvtkwriter.hh> where is the working version (martins?)
#include <dumux/multidimension/embeddedcoupling/circlegridfactory.hh>
#include <dumux/multidimension/embeddedcoupling/box/boxgridfactory.hh>
#include <dumux/multidimension/embeddedcoupling/pointsourcedata.hh>

namespace Dumux
{

namespace Properties
{
// Property forward declarations
NEW_PROP_TAG(BulkProblemTypeTag);
NEW_PROP_TAG(LowDimProblemTypeTag);
NEW_PROP_TAG(GridView);
} // namespace Properties

/*!
 * \ingroup EmbeddedCoupling
 * \brief Manages the coupling between bulk elements and lower dimensional elements
 */
template<class TypeTag>
class EmbeddedCouplingManager
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;

    // obtain the type tags of the sub problems
    typedef typename GET_PROP_TYPE(TypeTag, BulkProblemTypeTag) BulkProblemTypeTag;
    typedef typename GET_PROP_TYPE(TypeTag, LowDimProblemTypeTag) LowDimProblemTypeTag;

    typedef typename GET_PROP_TYPE(BulkProblemTypeTag, GridView) BulkGridView;
    typedef typename GET_PROP_TYPE(LowDimProblemTypeTag, GridView) LowDimGridView;

    typedef typename GET_PROP_TYPE(BulkProblemTypeTag, Grid) BulkGrid;
    typedef typename GET_PROP_TYPE(LowDimProblemTypeTag, Grid) LowDimGrid;

    typedef typename GET_PROP_TYPE(BulkProblemTypeTag, Problem) BulkProblem;
    typedef typename GET_PROP_TYPE(LowDimProblemTypeTag, Problem) LowDimProblem;

    typedef typename GET_PROP_TYPE(BulkProblemTypeTag, PointSource) BulkPointSource;
    typedef typename GET_PROP_TYPE(LowDimProblemTypeTag, PointSource) LowDimPointSource;

    typedef typename GET_PROP_TYPE(BulkProblemTypeTag, PrimaryVariables) BulkPrimaryVariables;
    typedef typename GET_PROP_TYPE(LowDimProblemTypeTag, PrimaryVariables) LowDimPrimaryVariables;

    enum {
        bulkDim = BulkGridView::dimension,
        lowDimDim = LowDimGridView::dimension,
        dimWorld = BulkGridView::dimensionworld
    };

    enum {
        bulkIsBox = GET_PROP_VALUE(BulkProblemTypeTag, ImplicitIsBox),
        lowDimIsBox = GET_PROP_VALUE(LowDimProblemTypeTag, ImplicitIsBox)
    };

    typedef Dumux::PointSourceData<TypeTag> PointSourceData;

    typedef typename BulkGridView::template Codim<0>::Entity BulkElement;
    typedef typename LowDimGridView::template Codim<0>::Entity LowDimElement;
    typedef typename Dune::GridGlue::Codim0Extractor<BulkGridView> BulkExtractor;
    typedef typename Dune::GridGlue::Codim0Extractor<LowDimGridView> LowDimExtractor;
    typedef typename Dune::GridGlue::GridGlue<BulkExtractor, LowDimExtractor> GlueType;

public:

    /*!
     * \brief Constructor
     */
    EmbeddedCouplingManager(BulkProblem& bulkProblem, LowDimProblem& lowDimProblem)
    : bulkProblem_(bulkProblem), lowDimProblem_(lowDimProblem) {
        circleCounter_ = 0;
    }

    /*!
     * \brief Methods to be accessed by the coupled problem
     */
    // \{

    /*!
     * \brief Called by the coupled problem before
     *        initializing the subproblems / models
     */
    void preInit()
    {
        computePointSourceData(1, false);
    }

    /*!
     * \brief Called by the coupled problem after
     *        initializing the subproblems / models
     */
    void postInit()
    {}

    const std::vector<unsigned int>& couplingStencil(const BulkElement& element) const
    {
        const unsigned int eIdx = bulkProblem().elementMapper().index(element);
        if (bulkCouplingStencils_.count(eIdx))
            return bulkCouplingStencils_.at(eIdx);
        else
            return emptyStencil_;
    }

    const std::vector<unsigned int>& couplingStencil(const LowDimElement& element) const
    {
        const unsigned int eIdx = lowDimProblem().elementMapper().index(element);
        if (lowDimCouplingStencils_.count(eIdx))
            return lowDimCouplingStencils_.at(eIdx);
        else
            return emptyStencil_;
    }

    const std::vector<unsigned int>& stencil(const BulkElement& element) const
    {
        const unsigned int eIdx = bulkProblem().elementMapper().index(element);
        return bulkStencils_.at(eIdx);
    }

    const std::vector<unsigned int>& stencil(const LowDimElement& element) const
    {
        const unsigned int eIdx = lowDimProblem().elementMapper().index(element);
        return lowDimStencils_.at(eIdx);
    }

    // \}

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
        auto bulkExtractor = std::make_shared<BulkExtractor>(this->bulkGridView(), bulkDescriptor);
        auto lowDimExtractor = std::make_shared<LowDimExtractor>(this->lowDimGridView(), lowDimDescriptor);
        // choose overlapping merger
        auto merger = std::make_shared<Dune::GridGlue::OverlappingMerge<bulkDim, lowDimDim, dimWorld>>();
        // build the glue object
        GlueType glue(bulkExtractor, lowDimExtractor, merger);

        // enable fallback brute-force step if the advancing front algorithm didn't find a seed
        merger->enableFallback(true);

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

        // empty the vectors
        clearPointSources_();

        // initialize the id counter
        idCounter_ = 0;

        // iterate over all intersections
        // we choose the target side (1, the low dimensional domain)
        const auto isEndIt = glue.template iend<1>();
        for (auto isIt = glue.template ibegin<1>(); isIt != isEndIt; ++isIt)
        {
            if (bulkIsBox && lowDimIsBox)
                treatIntersectionBox_(isIt, order, verbose);

            else if (!bulkIsBox && !lowDimIsBox)
                treatIntersectionCC_(isIt, order, verbose);

            else
                DUNE_THROW(Dune::NotImplemented, "Mixed discretization.");
        }

        std::cout << "took " << watch.elapsed() << " seconds." << std::endl;
    }

    void computeStencils()
    {
        // coupling stencils
        for (auto&& stencil : bulkCouplingStencils_)
        {
            std::sort(stencil.second.begin(), stencil.second.end());
            stencil.second.erase(std::unique(stencil.second.begin(), stencil.second.end()), stencil.second.end());
        }
        // the low dim coupling stencil
        for (auto&& stencil : lowDimCouplingStencils_)
        {
            std::sort(stencil.second.begin(), stencil.second.end());
            stencil.second.erase(std::unique(stencil.second.begin(), stencil.second.end()), stencil.second.end());
        }

        // compute the bulk stencil (like in the fv element geometry)
        for (const auto& element : elements(bulkGridView()))
        {
            if (bulkIsBox)
            {
                const unsigned int eIdx =  bulkProblem().elementMapper().index(element);
                for (unsigned int i = 0; i < element.subEntities(bulkDim); ++i)
                {
                    auto bulkVertexIdx = bulkProblem().vertexMapper().subIndex(element, i, bulkDim);
                    bulkStencils_[eIdx].push_back(bulkVertexIdx);
                }
            }
            else
            {
                // the element itself is first
                const unsigned int eIdx =  bulkProblem().elementMapper().index(element);
                bulkStencils_.insert({eIdx, {eIdx}});

                for (const auto& intersection : intersections(bulkGridView(), element))
                    if (intersection.neighbor())
                        bulkStencils_.at(eIdx).push_back(bulkProblem().elementMapper().index(intersection.outside()));
            }
        }

        // compute the low dimensional stencil (like in the fv geometry)
        for (const auto& element : elements(lowDimGridView()))
        {
            if (lowDimIsBox)
            {
                const unsigned int eIdx =  lowDimProblem().elementMapper().index(element);
                for (unsigned int i = 0; i < element.subEntities(lowDimDim); ++i)
                {
                    auto lowDimVertexIdx = lowDimProblem().vertexMapper().subIndex(element, i, lowDimDim);
                    lowDimStencils_[eIdx].push_back(lowDimVertexIdx);
                }
            }
            else
            {
                // the element itself is first
                const unsigned int eIdx =  lowDimProblem().elementMapper().index(element);
                lowDimStencils_.insert({eIdx, {eIdx}});

                for (const auto& intersection : intersections(lowDimGridView(), element))
                    if (intersection.neighbor())
                        lowDimStencils_.at(eIdx).push_back(lowDimProblem().elementMapper().index(intersection.outside()));
            }
        }
    }

    /*!
     * \brief Methods to be accessed by the subproblems
     */
    // \{

    //! Return a reference to the bulk problem
    Scalar radius(unsigned int id) const
    {
        const auto& data = pointSourceData_[id];
        return lowDimProblem().spatialParams().radius(data.lowDimElementIdx());
    }

    //! Return a reference to the pointSource data
    const PointSourceData& pointSourceData(unsigned int id) const
    {
        return pointSourceData_[id];
    }

    //! Return a reference to the bulk problem
    const BulkProblem& bulkProblem() const
    {
        return bulkProblem_;
    }

    //! Return a reference to the low dimensional problem
    const LowDimProblem& lowDimProblem() const
    {
        return lowDimProblem_;
    }

    //! Return a reference to the bulk problem finite element cache
    const Dune::PQkLocalFiniteElementCache<Scalar, Scalar, bulkDim, 1>&
    bulkFeCache() const
    {
        return bulkFeCache_;
    }

    //! Return a reference to the low dimensional finite element cache
    const Dune::PQkLocalFiniteElementCache<Scalar, Scalar, lowDimDim, 1>&
    lowDimFeCache() const
    {
        return lowDimFeCache_;
    }

    //! Return a reference to the bulk problem
    const BulkGridView& bulkGridView() const
    {
        return bulkProblem().gridView();
    }

    //! Return a reference to the low dimensional problem
    const LowDimGridView& lowDimGridView() const
    {
        return lowDimProblem().gridView();
    }

    //! Return a reference to the point sources
    const std::vector<BulkPointSource>& bulkPointSources() const
    {
        return bulkPointSources_;
    }

    //! Return a reference to the low dimensional problem
    const std::vector<LowDimPointSource>& lowDimPointSources() const
    {
        return lowDimPointSources_;
    }

    //! Return data for a bulk point source with the identifier id
    const BulkPrimaryVariables& bulkPriVars(unsigned int id)
    {
        auto& data = pointSourceData_[id];
        data.interpolateBulk(bulkProblem().model().curSol());
        return data.bulkPriVars();
    }

    //! Return data for a low dim point source with the identifier id
    const LowDimPrimaryVariables& lowDimPriVars(unsigned int id)
    {
        auto& data = pointSourceData_[id];
        data.interpolateLowDim(lowDimProblem().model().curSol());
        return data.lowDimPriVars();
    }
    // \}

private:

    template <class IntersectionIterator>
    void treatIntersectionBox_(const IntersectionIterator& isIt,
                               unsigned int order = 1,
                               bool verbose = false)
    {
        // all inside elements are identical...
        const auto& inside = isIt->inside(0);
        const auto insideGeometry = inside.geometry();
        const auto& lowDimLocalFiniteElement = lowDimFeCache_.get(insideGeometry.type());

        // we need only one bulk element as interpolations are only dependent on the corner values
        const auto& outside = isIt->outside(0);
        const auto outsideGeometry = outside.geometry();
        const auto& bulkLocalFiniteElement = bulkFeCache_.get(outsideGeometry.type());

        // build the box grids
        // all inside elements are identical...
        std::shared_ptr<LowDimGrid> boxLowDimGrid(Dumux::LocalBoxGridFactory<LowDimGrid, lowDimDim>::createGrid(insideGeometry));
        // we only need to consider one outside element as they are conforming
        std::shared_ptr<BulkGrid> boxBulkGrid(Dumux::LocalBoxGridFactory<BulkGrid, bulkDim>::createGrid(outsideGeometry));

        // redirect grid-glue output to file
        std::ofstream out("grid-glue.log");
        auto coutbuf = std::cout.rdbuf(out.rdbuf()); //save and redirect

        auto bulkDescriptor = [](const BulkElement& element, unsigned int subentity) -> bool { return true; };
        auto lowDimDescriptor = [](const LowDimElement& element, unsigned int subentity) -> bool { return true; };
        auto localLowDimExtractor = std::make_shared<LowDimExtractor>(boxLowDimGrid->leafGridView(), lowDimDescriptor);
        auto localBulkExtractor = std::make_shared<BulkExtractor>(boxBulkGrid->leafGridView(), bulkDescriptor);
        auto localMerger = std::make_shared<Dune::GridGlue::OverlappingMerge<bulkDim, lowDimDim, dimWorld>>();
        GlueType localGlue(localBulkExtractor, localLowDimExtractor, localMerger);
        // calculate box intersections
        localGlue.build();

        //reset to standard output again
        std::cout.rdbuf(coutbuf);

        // loop over local intersections
        const auto localIsEndIt = localGlue.template iend<1>();
        for (auto localIsIt = localGlue.template ibegin<1>(); localIsIt != localIsEndIt; ++localIsIt)
        {
            const auto intersectionGeometry = localIsIt->geometry();

            // get the Gaussian quadrature rule for the local intersection
            typedef Dune::QuadratureRule<Scalar, lowDimDim> Quad;
            const Quad &quad = Dune::QuadratureRules<Scalar, lowDimDim>::rule(intersectionGeometry.type(), order);

            // iterate over all quadrature points
            for (auto&& qp : quad)
            {
                // each quadrature point will be a point source for the sub problem
                const auto globalPos = intersectionGeometry.global(qp.position());
                const auto id = idCounter_++;
                const Scalar value = qp.weight()*intersectionGeometry.integrationElement(qp.position());
                bulkPointSources_.emplace_back(globalPos, value, id);
                lowDimPointSources_.emplace_back(globalPos, value, id);

                // pre compute additional data used for the evaluation of
                // the actual solution dependent source term
                PointSourceData psData;

                std::vector<Dune::FieldVector<Scalar, 1> > lowDimShapeValues, bulkShapeValues;
                std::vector<unsigned int> lowDimCornerIndices, bulkCornerIndices;

                lowDimLocalFiniteElement.localBasis().evaluateFunction(insideGeometry.local(globalPos), lowDimShapeValues);
                bulkLocalFiniteElement.localBasis().evaluateFunction(outsideGeometry.local(globalPos), bulkShapeValues);

                for (unsigned int cIdx = 0; cIdx < insideGeometry.corners(); ++cIdx)
                    lowDimCornerIndices.push_back(lowDimGridView().indexSet().index(inside.template subEntity<lowDimDim>(cIdx)));
                for (unsigned int cIdx = 0; cIdx < outsideGeometry.corners(); ++cIdx)
                    bulkCornerIndices.push_back(bulkGridView().indexSet().index(outside.template subEntity<bulkDim>(cIdx)));

                const unsigned int lowDimElementIdx = lowDimProblem().elementMapper().index(inside);
                const unsigned int bulkElementIdx = bulkProblem().elementMapper().index(outside);

                psData.addLowDimInterpolation(lowDimShapeValues, lowDimCornerIndices, lowDimElementIdx);
                psData.addBulkInterpolation(bulkShapeValues, bulkCornerIndices, bulkElementIdx);

                // publish point source data in the global vector
                pointSourceData_.push_back(psData);
            }
        }

        auto insideIdx = lowDimProblem().elementMapper().index(inside);
        // compute the coupling stencils
        for (unsigned int outsideIdx = 0; outsideIdx < isIt->neighbor(0); ++outsideIdx)
        {
            const unsigned int outsideBulkElementIdx = bulkProblem().elementMapper().index(isIt->outside(outsideIdx));
            for (unsigned int i = 0; i < inside.subEntities(lowDimDim); ++i)
            {
                const auto lowDimVertexIdx = lowDimProblem().vertexMapper().subIndex(inside, i, lowDimDim);
                bulkCouplingStencils_[outsideBulkElementIdx].push_back(lowDimVertexIdx);
            }

            auto bulkElement = isIt->outside(outsideIdx);
            for (unsigned int i = 0; i < bulkElement.subEntities(bulkDim); ++i)
            {
                const auto bulkVertexIdx = bulkProblem().vertexMapper().subIndex(bulkElement, i, bulkDim);
                lowDimCouplingStencils_[insideIdx].push_back(bulkVertexIdx);
            }
        }
    }

    template <class IntersectionIterator>
    void treatIntersectionCC_(IntersectionIterator& isIt,
                              unsigned int order = 1,
                              bool verbose = false)
    {

        // all inside elements are identical...
        const auto& inside = isIt->inside(0);
        // we need only one bulk element as interpolations are only dependent on the corner values
        const auto& outside = isIt->outside(0);
        // get the intersection geometry for integrating over it
        const auto intersectionGeometry = isIt->geometry();

        // get the Gaussian quadrature rule for the local intersection
        typedef Dune::QuadratureRule<Scalar, lowDimDim> Quad;
        const Quad &quad = Dune::QuadratureRules<Scalar, lowDimDim>::rule(intersectionGeometry.type(), order);

        const unsigned int lowDimElementIdx = lowDimProblem().elementMapper().index(inside);
        const unsigned int bulkElementIdx = bulkProblem().elementMapper().index(outside);

        // iterate over all quadrature points
        for (auto&& qp : quad)
        {
            // each quadrature point will be a point source for the sub problem
            const auto globalPos = intersectionGeometry.global(qp.position());
            const auto id = idCounter_++;
            const Scalar value = qp.weight()*intersectionGeometry.integrationElement(qp.position());
            bulkPointSources_.push_back(BulkPointSource(globalPos, value, id));
            lowDimPointSources_.push_back(LowDimPointSource(globalPos, value, id));

            // pre compute additional data used for the evaluation of
            // the actual solution dependent source term
            PointSourceData psData;
            psData.addLowDimInterpolation(lowDimElementIdx);
            psData.addBulkInterpolation(bulkElementIdx);

            // publish point source data in the global vector
            pointSourceData_.push_back(psData);
        }

        // compute the coupling stencils
        for (unsigned int outsideIdx = 0; outsideIdx < isIt->neighbor(0); ++outsideIdx)
        {
            const unsigned int outsideBulkElementIdx = bulkProblem().elementMapper().index(isIt->outside(outsideIdx));
            bulkCouplingStencils_[outsideBulkElementIdx].push_back(lowDimElementIdx);
            lowDimCouplingStencils_[lowDimElementIdx].push_back(outsideBulkElementIdx);
        }
    }

    void clearPointSources_()
    {
        bulkPointSources_.clear();
        lowDimPointSources_.clear();
        pointSourceData_.clear();
    }

    BulkProblem& bulkProblem_;
    LowDimProblem& lowDimProblem_;

    std::vector<BulkPointSource> bulkPointSources_;
    std::vector<LowDimPointSource> lowDimPointSources_;

    std::vector<PointSourceData> pointSourceData_;

    std::unordered_map<unsigned int, std::vector<unsigned int> > bulkCouplingStencils_;
    std::unordered_map<unsigned int, std::vector<unsigned int> > lowDimCouplingStencils_;
    std::unordered_map<unsigned int, std::vector<unsigned int> > bulkStencils_;
    std::unordered_map<unsigned int, std::vector<unsigned int> > lowDimStencils_;
    std::unordered_map<unsigned int, std::vector<unsigned int> > bulkCircleStencils_;
    std::unordered_map<unsigned int, std::vector<unsigned int> > lowDimCircleStencils_;
    std::vector<unsigned int> emptyStencil_;

    // id generator for point sources
    unsigned int idCounter_;

    unsigned int circleCounter_;

    // finite element caches for ansatz functions (only box)
    const Dune::PQkLocalFiniteElementCache<Scalar, Scalar, bulkDim, 1> bulkFeCache_;
    const Dune::PQkLocalFiniteElementCache<Scalar, Scalar, lowDimDim, 1> lowDimFeCache_;
};

} // namespace Dumux

#endif
