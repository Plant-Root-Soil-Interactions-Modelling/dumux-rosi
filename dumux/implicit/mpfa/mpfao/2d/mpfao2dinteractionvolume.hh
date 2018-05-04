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
#ifndef DUMUX_MPFAO2DINTERACTIONVOLUME_HH
#define DUMUX_MPFAO2DINTERACTIONVOLUME_HH

#include <eigen3/Eigen/Dense>

/**
 * @file
 * @brief  Class including the information of an interaction volume of a MPFA O-method that does not change with time.
 */

namespace Dumux
{

//!
/*! \brief Class including the information of an interaction volume of a MPFA O-method that does not change with time.
 *
 * Includes information needed to calculate the transmissibility matrix of an O-interaction-volume.
 *
 */
template<class TypeTag>
class MpfaO2DInteractionVolume
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP(TypeTag, MpfaMethods) MpfaMethods;

    enum
    {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;

    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> WorldMatrix;

public:

    typedef Dune::FieldVector<Scalar, dimWorld> DimVector;
    typedef Eigen::MatrixXd DynamicMatrix;
    typedef Eigen::VectorXd DynamicVector;

    struct Properties
    {
        // In an interior 2D Mpfa-O interaction volume
        // for quadrilaterals there are 9 elements in stencil
        // This value is used for a primary resizing of involved
        // vectors.
        static const int StandardStencilSize = 9;
    };

    struct FaceTypes
    {
        static const int InteriorFace = 0;

        static const int DirichletFace = 1;

        static const int NeumannFace = 2;

        static const int InternalDirichletFace = 3;

        static const int InternalNoFlowFace = 4;

        static const int InternalFluxFace = 5;
    };

    struct SubVolumeFace
    {
        DimVector normal;
        DimVector center;
        Scalar area;
        int positiveSubVolIndex;
        int negativeSubVolIndex;
        int localIdxOnPositive;
        int localIdxOnNegative;
        bool onBoundary;
        int faceType;
        BoundaryTypes boundaryTypes;

        SubVolumeFace()
        {
            normal = 0;
            area = 0;
            positiveSubVolIndex = -1;
            negativeSubVolIndex = -1;
            localIdxOnPositive = -1;
            localIdxOnNegative = -1;
            onBoundary = false;
            faceType = -1;
            boundaryTypes.reset();
        }
    };

    struct SubVolume
    {
        Element element;
        int globalIdx;
        Dune::FieldVector<int, dim> localIsIndex;
        Dune::FieldVector<DimVector, dim> nu;
        Dune::FieldVector<DimVector, dim> contiPoints;
        Scalar T;

        SubVolume()
        {
            element = Element();
            globalIdx = -1;
            nu = DimVector(0);
            contiPoints = DimVector(0);
            for(int i = 0; i < dim; i++)
                localIsIndex[i] = -1;
            T = 0;
        }
    };

    //! Constructs an MpfaO2DInteractionVolume object
    // initialization of the variables according to the
    // case of an interior interaction volume in 2D
    // for quadrilaterals
    MpfaO2DInteractionVolume():
        stored_(false),
        interiorVolume_(true),
        noOfFaces_(0),
        noOfSubVols_(0),
        interactionVolumeType_(MpfaMethods::oMethod)
    {
        subVolumes_.reserve(4);
        subVolumeFaces_.reserve(4);
    }

    int interactionVolumeType() const
    {
        return interactionVolumeType_;
    }

    //! Mark storage as completed
    void setStored()
    {
        stored_ = true;
    }

    // set the coordinates of the central vertex of the interaction volume
    void setCentralVertex(DimVector centralVertex)
    {
        centralVertex_ = centralVertex;
    }

    // set the global index of the central vertex of the interaction volume
    void setCentralVertexIndex(int globalIndex)
    {
        centralVertexIndex_ = globalIndex;
    }

    DimVector getCentralVertex() const
    {
        return centralVertex_;
    }

    int getCentralVertexGlobalIndex() const
    {
        return centralVertexIndex_;
    }

    //! Returns true if information has already been stored
    bool isStored() const
    {
        return stored_;
    }

    void setSubVolToFaceMap(int subVolIdx, int localFaceIdx, int faceIndex)
    {
        std::pair<int, int> key;
        key = std::make_pair(subVolIdx, localFaceIdx);

        facesOnSubVolumes_.insert(std::pair< std::pair<int,int>, int> (key, faceIndex));
    }

    int getFaceIndexFromVolume(int subVolIdx, int localFaceIdx)
    {
        std::pair<int, int> key;
        key = std::make_pair(subVolIdx, localFaceIdx);

        std::map< std::pair<int, int>, int >::iterator it = facesOnSubVolumes_.find(key);
        if (it != facesOnSubVolumes_.end())
            return it->second;
        else
            DUNE_THROW(Dune::InvalidStateException, "No face index found corresponding to local subvolume face index!");
    }

    int getFaceIndexFromVolume(int subVolIdx, int localFaceIdx) const
    {
        std::pair<int, int> key;
        key = std::make_pair(subVolIdx, localFaceIdx);

        std::map< std::pair<int, int>, int >::const_iterator it = facesOnSubVolumes_.find(key);
        if (it != facesOnSubVolumes_.end())
            return it->second;
        else
            DUNE_THROW(Dune::InvalidStateException, "No face index found corresponding to local subvolume face index!");
    }

    void setElementToSubVolMap(int globalElIdx, int subVolIdx)
    {
        std::pair<int, int> indexPair (globalElIdx, subVolIdx);
        globalToSubVolIdx_.insert(indexPair);
    }

    int findSubVolumeIndex(int globalElementIndex) const
    {
        std::map<int, int>::const_iterator it = globalToSubVolIdx_.find(globalElementIndex);
        if (it != globalToSubVolIdx_.end())
            return it->second;
        else
            DUNE_THROW(Dune::InvalidStateException, "No sub volume index found corresponding to global element index!");
    }

    int getFaceIndexInRegion(const Problem &problem, const Element &element, int facetIdxOnElement, bool &switchSign, int &regionFaceIdx) const
    {
        int elementGlobalIdx = problem.elementMapper().index(element);

        // find out sub volume index of the element at hand
        int subVolIdx;
        std::map <int, int>::const_iterator it = globalToSubVolIdx_.find(elementGlobalIdx);
        if (it != globalToSubVolIdx_.end())
            subVolIdx = it->second;
        else
            DUNE_THROW(Dune::InvalidStateException, "No sub volume index found corresponding to global element index!");

        // now find out the local sub volume face index of the element facet
        int facetIdxOnSubVolume = getFacetIndexOnSubVolume(subVolIdx, facetIdxOnElement);

        // now get the corresponding face index within the interaction region
        int faceIdxInRegion = getFaceIndexFromVolume(subVolIdx, facetIdxOnSubVolume);
        // pass the index to the provided container
        regionFaceIdx = faceIdxInRegion;

        // if normal vector is defined pointing inside the element, set boolean switchSign to true, otherwise to false
        if (subVolumeFaces_[faceIdxInRegion].positiveSubVolIndex == subVolIdx)
            switchSign = false;
        else if (subVolumeFaces_[faceIdxInRegion].negativeSubVolIndex == subVolIdx)
            switchSign = true;
        else
            DUNE_THROW(Dune::InvalidStateException, "Check implementation!!!");

        return faceIdxInRegion;
    }

    int getFacetIndexOnSubVolume(int subVolIdx, int facetIdxOnElement) const
    {
        if (subVolumes_[subVolIdx].localIsIndex[0] == facetIdxOnElement)
            return 0;
        else if (subVolumes_[subVolIdx].localIsIndex[1] == facetIdxOnElement)
            return 1;
        else
            DUNE_THROW(Dune::InvalidStateException, "No sub volume face index found corresponding to element facet index!");
    }

    int getInteriorFaceIndex(int faceIndex)
    {
        std::set<int>::iterator it = interiorFaces_.begin();
        for(int i = 0; i < interiorFaces_.size(); i++)
        {
            if(*it == faceIndex)
                return i;

            ++it;
        }
        // if we are here we haven't found the face which should not happen!
        DUNE_THROW(Dune::InvalidStateException, "Face index not found in set of interior faces!");
    }

    int getInteriorFaceIndex(int faceIndex) const
    {
        std::set<int>::const_iterator it = interiorFaces_.begin();
        for(int i = 0; i < interiorFaces_.size(); i++)
        {
            if(*it == faceIndex)
                return i;

            ++it;
        }
        // if we are here we haven't found the face which should not happen!
        DUNE_THROW(Dune::InvalidStateException, "Face index not found in set of interior faces!");
    }

    int getDirichletFaceIndex(int faceIndex)
    {
        std::set<int>::iterator it = dirichletFaces_.begin();
        for(int i = 0; i < dirichletFaces_.size(); i++)
        {
            if(*it == faceIndex)
                return i;

            ++it;
        }
        // if we are here we haven't found the face which should not happen!
        DUNE_THROW(Dune::InvalidStateException, "Face index not found in set of boundary faces!");
    }

    int getDirichletFaceIndex(int faceIndex) const
    {
        std::set<int>::iterator it = dirichletFaces_.begin();
        for(int i = 0; i < dirichletFaces_.size(); i++)
        {
            if(*it == faceIndex)
                return i;

            ++it;
        }
        // if we are here we haven't found the face which should not happen!
        DUNE_THROW(Dune::InvalidStateException, "Face index not found in set of boundary faces!");
    }

    std::set<int>& getDirichletFaceIndexSet()
    {
        return dirichletFaces_;
    }
    const std::set<int>& getDirichletFaceIndexSet() const
    {
        return dirichletFaces_;
    }

    std::set<int>& getInteriorFaceIndexSet()
    {
        return interiorFaces_;
    }
    const std::set<int>& getInteriorFaceIndexSet() const
    {
        return interiorFaces_;
    }

    void setNumberOfSubVols(int input)
    {
        noOfSubVols_ = input;
    }

    int getNumberOfSubVols() const
    {
        return noOfSubVols_;
    }

    void setNumberOfSubFaces(int input)
    {
        noOfFaces_ = input;
    }

    int getNumberOfSubFaces() const
    {
        return noOfFaces_;
    }

    int getNumberOfInteriorFaces()
    {
        return interiorFaces_.size();
    }
    int getNumberOfInteriorFaces() const
    {
        return interiorFaces_.size();
    }

    void addSubVolume(SubVolume &input)
    {
        subVolumes_.push_back(input);
    }

    SubVolume& getSubVolume(int subVolIdx)
    {
        return subVolumes_[subVolIdx];
    }
    const SubVolume& getSubVolume(int subVolIdx) const
    {
        return subVolumes_[subVolIdx];
    }

    void addSubVolumeFace(SubVolumeFace &input)
    {
        subVolumeFaces_.push_back(input);
    }

    SubVolumeFace& getSubVolumeFace(int subVolFaceIdx)
    {
        return subVolumeFaces_[subVolFaceIdx];
    }
    const SubVolumeFace& getSubVolumeFace(int subVolFaceIdx) const
    {
        return subVolumeFaces_[subVolFaceIdx];
    }

    void addInteriorFace(int subFaceIndex)
    {
        interiorFaces_.insert(subFaceIndex);
    }

    void addDirichletFace(int subFaceIndex)
    {
        dirichletFaces_.insert(subFaceIndex);
    }

    void setInteriorVolume(bool input)
    {
        interiorVolume_ = input;
    }

    bool isInteriorVolume() const
    {
        return interiorVolume_;
    }

    int getNumberOfBoundaryFaces() const
    {
        return dirichletFaces_.size();
    }

    DynamicMatrix& getTransMatrix()
    {
        return T_;
    }

    const DynamicMatrix& getTransMatrix() const
    {
        return T_;
    }

    // Ainverse * B can be used later to evaluate gradU (U = primary variable)
    // in order to interpolate on the sub control volume
    DynamicMatrix& getAInvB()
    {
        return AInverseB_;
    }
    const DynamicMatrix& getAInvB() const
    {
        return AInverseB_;
    }

    // C * Ainverse is multiplied on the neumann/interior fluxes for the application
    // of the neumann or interior flux boundary conditions
    DynamicMatrix& getCAInverse()
    {
        return CAInverse_;
    }
    const DynamicMatrix& getCAInverse() const
    {
        return CAInverse_;
    }

    DynamicMatrix& getAInverse()
    {
        return AInverse_;
    }
    const DynamicMatrix& getAInverse() const
    {
        return AInverse_;
    }

    Scalar getW(const Problem &problem, int faceIdx, int subVolIdx, int localDirection) const
    {
        WorldMatrix K = problem.model().getDiffusionCoefficient(subVolumes_[subVolIdx].element);

        DimVector tmp(0);
        DimVector normal = subVolumeFaces_[faceIdx].normal;
        DimVector nu = subVolumes_[subVolIdx].nu[localDirection];

        K.mv(nu, tmp);
        tmp /= subVolumes_[subVolIdx].T;
        tmp *= subVolumeFaces_[faceIdx].area;
        tmp *= -1;
        return normal * tmp;
    }

    void calculateTransmissivityMatrix(const Problem &problem)
    {
        const int noInteriorFaces = interiorFaces_.size();
        const int noDirichletFaces = dirichletFaces_.size();
        const int noOfPotentials = noOfSubVols_ + noDirichletFaces;

        T_.resize(noOfFaces_, noOfPotentials);
        T_ = Eigen::MatrixXd::Zero(noOfFaces_, noOfPotentials);

        AInverseB_.resize(noInteriorFaces, noOfPotentials);
        AInverseB_ = Eigen::MatrixXd::Zero(noInteriorFaces, noOfPotentials);

        CAInverse_.resize(noOfFaces_, noInteriorFaces);
        CAInverse_ = Eigen::MatrixXd::Zero(noOfFaces_, noInteriorFaces);

        if (noInteriorFaces == 0)
        {
            // Loop over all the faces, in this case these are all dirichlet boundaries
            for(int face = 0; face < subVolumeFaces_.size(); face++)
            {
                int posSubVol = subVolumeFaces_[face].positiveSubVolIndex;
                if (posSubVol == -1)
                    DUNE_THROW(Dune::NotImplemented, "posSubVolIdx = -1!!");

                Scalar WiPosj[dim];
                for (int i = 0; i < dim; i++)
                    WiPosj[i] = getW(problem, face, posSubVol, i);

                for (int localDir = 0; localDir < dim; localDir++)
                {
                    int globalFaceIndex = getFaceIndexFromVolume(posSubVol,localDir);
                    int dirichletFaceIdx = getDirichletFaceIndex(globalFaceIndex);
                    T_(face, noOfSubVols_ + dirichletFaceIdx) += WiPosj[localDir];
                    T_(face, posSubVol) -= WiPosj[localDir];
                }
                // if face is an interior dirichlet face, we have to add the entries coming from the right hand side
                if (subVolumeFaces_[face].faceType == FaceTypes::InternalDirichletFace)
                {
                    int negSubVol = subVolumeFaces_[face].negativeSubVolIndex;
                    if (negSubVol == -1)
                        DUNE_THROW(Dune::NotImplemented, "negSubVolIdx = -1!!");

                    Scalar WiNegj[dim];
                    for (int i = 0; i < dim; i++)
                        WiNegj[i] = getW(problem, face, negSubVol, i);

                    for (int localDir = 0; localDir < dim; localDir++)
                    {
                        int globalFaceIndex = getFaceIndexFromVolume(negSubVol,localDir);
                        int dirichletFaceIdx = getDirichletFaceIndex(globalFaceIndex);
                        T_(face, noOfSubVols_ + dirichletFaceIdx) -= WiNegj[localDir];
                        T_(face, negSubVol) += WiNegj[localDir];
                    }
                }
            }
        }
        else
        {
            DynamicMatrix A = Eigen::MatrixXd::Zero(noInteriorFaces, noInteriorFaces);
            DynamicMatrix B = Eigen::MatrixXd::Zero(noInteriorFaces, noOfPotentials);
            DynamicMatrix C = Eigen::MatrixXd::Zero(noOfFaces_, noInteriorFaces);
            DynamicMatrix D = Eigen::MatrixXd::Zero(noOfFaces_, noOfPotentials);

            Scalar xi = GET_PROP_VALUE(TypeTag, XiFactor);

            // Loop over all the faces, assemble the matrices
            for(int face = 0; face < subVolumeFaces_.size(); face++)
            {
                int posSubVol = subVolumeFaces_[face].positiveSubVolIndex;
                if (posSubVol == -1)
                    DUNE_THROW(Dune::NotImplemented, "posSubVolIdx = -1!!");

                bool hasFaceUnknown = false;
                if ( subVolumeFaces_[face].faceType != FaceTypes::DirichletFace
                            && subVolumeFaces_[face].faceType != FaceTypes::InternalDirichletFace )
                    hasFaceUnknown = true;

                // Index of face within the interior faces
                int interiorIndexFace = -1;
                if (hasFaceUnknown)
                    interiorIndexFace = getInteriorFaceIndex(face);

                Scalar WiPosj[dim];
                for (int i = 0; i < dim; i++)
                    WiPosj[i] = getW(problem, face, posSubVol, i);

                // Check the local directions of positive sub volume
                for (int localDir = 0; localDir < dim; localDir++)
                {
                    int globalFaceIndex = getFaceIndexFromVolume(posSubVol,localDir);
                    if ( subVolumeFaces_[globalFaceIndex].faceType != FaceTypes::DirichletFace
                            && subVolumeFaces_[globalFaceIndex].faceType != FaceTypes::InternalDirichletFace)
                    {
                        // get index of face within interior face indices
                        int interiorIndexNeighbourFace = getInteriorFaceIndex(globalFaceIndex);
                        if (subVolumeFaces_[face].faceType != FaceTypes::InternalNoFlowFace)
                            C(face, interiorIndexNeighbourFace) += WiPosj[localDir];
                        if (hasFaceUnknown)
                        {
                            if (subVolumeFaces_[face].faceType == FaceTypes::InternalFluxFace)
                            {
                                A(interiorIndexFace, interiorIndexNeighbourFace) += xi*WiPosj[localDir];
                            }
                            else
                                A(interiorIndexFace, interiorIndexNeighbourFace) += WiPosj[localDir];
                        }
                        // if we are on an internal flux face and facet coupling is turned on,
                        // add entry coming from the fracture domain
                        if (hasFaceUnknown && interiorIndexNeighbourFace == interiorIndexFace
                            && subVolumeFaces_[face].faceType == FaceTypes::InternalFluxFace
                            && GET_PROP_VALUE(TypeTag, FacetCoupling)
                            && GET_PROP_VALUE(TypeTag, CouplingStrategy) == GET_PROP_TYPE(TypeTag, CouplingStrategies)::fluxCoupling)
                            {

                                Scalar facetPerm = problem.getLowDimPermeability(subVolumes_[posSubVol].element, subVolumes_[posSubVol].localIsIndex[localDir]);
                                Scalar aperture = problem.getLowDimDomainAperture(subVolumes_[posSubVol].element, subVolumes_[posSubVol].localIsIndex[localDir]);

                                A(interiorIndexFace, interiorIndexNeighbourFace) -= 2*facetPerm*subVolumeFaces_[face].area/aperture;
                            }
                    }
                    else
                    // if we are here that means that a face is a dirichlet boundary and therefore
                    // creates an entry in matrix D & eventually B
                    {
                        int dirichletFaceIdx = getDirichletFaceIndex(globalFaceIndex);
                        if (subVolumeFaces_[face].faceType != FaceTypes::InternalNoFlowFace)
                            D(face, noOfSubVols_ + dirichletFaceIdx) += WiPosj[localDir];
                        if (hasFaceUnknown)
                            B(interiorIndexFace, noOfSubVols_ + dirichletFaceIdx) -= WiPosj[localDir];
                    }
                    // add entries to matrix D
                    D(face, posSubVol) -= WiPosj[localDir];
                    if (hasFaceUnknown)
                    {
                        if (subVolumeFaces_[face].faceType == FaceTypes::InternalFluxFace)
                            B(interiorIndexFace, posSubVol) += xi*WiPosj[localDir];
                        else
                            B(interiorIndexFace, posSubVol) += WiPosj[localDir];
                    }
                }

                // add remaining entries to matrices A & B for the negative subVolIndex
                // in case we are on an interior face
                if (!subVolumeFaces_[face].onBoundary
                            && subVolumeFaces_[face].faceType != FaceTypes::InternalDirichletFace)
                {
                    int negSubVol = subVolumeFaces_[face].negativeSubVolIndex;
                    if (negSubVol == -1)
                        DUNE_THROW(Dune::NotImplemented, "negSubVolIdx = -1!!");

                    Scalar WiNegj[dim];
                    for (int i = 0; i < dim; i++)
                        WiNegj[i] = getW(problem, face, negSubVol, i);

                    // Check local directions of negative sub volume
                    for (int localDir = 0; localDir < dim; localDir++)
                    {
                        int globalFaceIndexNeg = getFaceIndexFromVolume(negSubVol,localDir);

                        if ( subVolumeFaces_[globalFaceIndexNeg].faceType != FaceTypes::DirichletFace
                                && subVolumeFaces_[globalFaceIndexNeg].faceType != FaceTypes::InternalDirichletFace)
                        {
                            // get index of face within interior face indices
                            int interiorIndexNeighbourFace = getInteriorFaceIndex(globalFaceIndexNeg);
                            if (subVolumeFaces_[face].faceType == FaceTypes::InternalFluxFace)
                                A(interiorIndexFace, interiorIndexNeighbourFace) -= (1-xi)*WiNegj[localDir];
                            else
                                A(interiorIndexFace, interiorIndexNeighbourFace) -= WiNegj[localDir];
                        }
                        else
                        // if we are here that means that a face is a dirichlet boundary and therefore
                        // creates an entry in matrix B
                            B(interiorIndexFace, noOfSubVols_ + getDirichletFaceIndex(globalFaceIndexNeg)) += WiNegj[localDir];

                        // add entries to matrix B
                        if (subVolumeFaces_[face].faceType == FaceTypes::InternalFluxFace)
                            B(interiorIndexFace, negSubVol) -= (1-xi)*WiNegj[localDir];
                        else
                            B(interiorIndexFace, negSubVol) -= WiNegj[localDir];
                    }
                }
            }
            //printMatrices(A, B, C, D);
            // Matrices are assembled, now we calculate the transmissivity matrix
            //Eigen::MatrixXd Ainverse(noInteriorFaces, noInteriorFaces);
            AInverse_ = A.inverse();
            AInverseB_ = AInverse_ * B;
            CAInverse_ = C * AInverse_;

            T_ = C * AInverseB_;
            T_ = T_ + D;
        }
        //printTransmissivityMatrix();
    }

    void passElementsInRegion(std::vector<Element> &storage)
    {
        storage.resize(subVolumes_.size());
        for (int i = 0; i < subVolumes_.size(); i++)
            storage[i] = subVolumes_[i].element;
    }

    void printMatrices(Eigen::MatrixXd &A, Eigen::MatrixXd &B, Eigen::MatrixXd &C, Eigen::MatrixXd &D)
    {
        std::cout << "Printing matrices A, B, C, D of the interaction volume..." << std::endl;

        std::cout << "Matrix A: " << std::endl;
        for(int i = 0; i < A.rows();i++)
        {
            for(int j = 0; j < A.cols(); j++)
                std::cout << A(i,j) << ", ";
            std::cout << std::endl;
        }
        std::cout << "Matrix B: " << std::endl;
        for(int i = 0; i < B.rows();i++)
        {
            for(int j = 0; j < B.cols(); j++)
                std::cout << B(i,j) << ", ";
            std::cout << std::endl;
        }
        std::cout << "Matrix C: " << std::endl;
        for(int i = 0; i < C.rows();i++)
        {
            for(int j = 0; j < C.cols(); j++)
                std::cout << C(i,j) << ", ";
            std::cout << std::endl;
        }
        std::cout << "Matrix D: " << std::endl;
        for(int i = 0; i < D.rows();i++)
        {
            for(int j = 0; j < D.cols(); j++)
                std::cout << D(i,j) << ", ";
            std::cout << std::endl;
        }
    }

    void printTransmissivityMatrix()
    {
        std::cout << "Printing transmissivity matrix: " << std::endl;
        for(int i = 0; i < T_.rows();i++)
        {
            for(int j = 0; j < T_.cols(); j++)
                std::cout << T_(i, j) << ", ";
            std::cout << std::endl;
        }
        std::cout << "Printing CA-1 matrix: " << std::endl;
        for(int i = 0; i < CAInverse_.rows();i++)
        {
            for(int j = 0; j < CAInverse_.cols(); j++)
                std::cout << CAInverse_(i, j) << ", ";
            std::cout << std::endl;
        }
    }

    void printTransmissivityMatrix() const
    {
        std::cout << "Printing transmissivity matrix: " << std::endl;
        for(int i = 0; i < T_.rows();i++)
        {
            for(int j = 0; j < T_.cols(); j++)
                std::cout << T_(i, j) << ", ";
            std::cout << std::endl;
        }
        std::cout << "Printing CA-1 matrix: " << std::endl;
        for(int i = 0; i < CAInverse_.rows();i++)
        {
            for(int j = 0; j < CAInverse_.cols(); j++)
                std::cout << CAInverse_(i, j) << ", ";
            std::cout << std::endl;
        }
    }

    void printInteractionVolume()
    {

        std::cout << std::endl << "Printing interaction volume" << std::endl;
        std::cout << "Number of elements in interaction volume:" << subVolumes_.size() << ", or noOfSubVols: " << noOfSubVols_ << std::endl
                    << "The number of faces is: " << noOfFaces_ << ", or size of FaceContainer: " << subVolumeFaces_.size() << std::endl
                    << "...of which " << interiorFaces_.size() << " are interior/neumann and " << dirichletFaces_.size() << " are dirichlet faces." << std::endl
                    << "The boundary face indices are: ";
                    for(std::set<int>::iterator it = dirichletFaces_.begin() ; it != dirichletFaces_.end();++it)
                        std::cout << *it << ", ";
                    std::cout << std::endl << "The interior face indices are: ";
                    for(std::set<int>::iterator it = interiorFaces_.begin(); it != interiorFaces_.end();++it)
                        std::cout << *it << ", ";
                    std::cout << std::endl;

        for(int i = 0; i < subVolumes_.size();i++)
        {
            int subVolIdx = i;
            std::cout << "Printing Info related to sub volume " << i << std::endl;
            std::cout << "ElementCenter: " << subVolumes_[i].element->geometry().center() << std::endl;
            std::cout << "ContiPoint" << subVolIdx+1 << "-1: " << subVolumes_[i].contiPoints[0] << std::endl
                        << "ContiPoint" << subVolIdx+1 << "-2: " << subVolumes_[i].contiPoints[1] << std::endl
                        << "Nu" << subVolIdx+1 << "-1: " << subVolumes_[i].nu[0] << std::endl
                        << "Nu" << subVolIdx+1 << "-2: " << subVolumes_[i].nu[1] << std::endl
                        << "FaceIndex of local face 1: " << getFaceIndexFromVolume(i,0) << std::endl
                        << "FaceIndex of local face 2: " << getFaceIndexFromVolume(i,1) << std::endl
                        << "FaceArea" << subVolIdx+1 << "-1: " << subVolumeFaces_[getFaceIndexFromVolume(i,0)].area << std::endl
                        << "FaceArea" << subVolIdx+1 << "-2: " << subVolumeFaces_[getFaceIndexFromVolume(i,1)].area << std::endl
                        << "FaceNormal" << subVolIdx+1 << "-1: " << subVolumeFaces_[getFaceIndexFromVolume(i,0)].normal << std::endl
                        << "FaceNormal" << subVolIdx+1 << "-2: " << subVolumeFaces_[getFaceIndexFromVolume(i,1)].normal << std::endl
                        << "T" << subVolIdx+1 << ": " << subVolumes_[i].T << std::endl;
        }
        std::cout << "Now we print the faceInfos: " << std::endl;
        for(int i = 0; i < subVolumeFaces_.size(); i++)
        {
            std::cout << "Printing info related to sub face " << i << std::endl;
            std::cout << "Positive sub volume Index: " << subVolumeFaces_[i].positiveSubVolIndex << std::endl
                        << "Negative sub volume Index: " << subVolumeFaces_[i].negativeSubVolIndex << std::endl
                        << "Is on boundary: " << subVolumeFaces_[i].onBoundary << std::endl;
        }
    }

    void updateSubVolumeData(int subVolIdx, DimVector newContiPoint0, int contiDirection, DimVector newNu1, int nuDirection, Scalar newT)
    {
        if (newContiPoint0 != subVolumes_[subVolIdx].contiPoints[contiDirection])
        {
            std::cout << "ContiPoints different! " << std::endl;
            std::cout << "ContiPoint before: " << subVolumes_[subVolIdx].contiPoints[contiDirection] << std::endl;
            std::cout << "ContiPoint after: " << newContiPoint0 << std::endl;
        }
        subVolumes_[subVolIdx].contiPoints[contiDirection] = newContiPoint0;

        if (newNu1 != subVolumes_[subVolIdx].nu[nuDirection])
        {
            std::cout << "Nu different! " << std::endl;
            std::cout << "nu before: " << subVolumes_[subVolIdx].nu[nuDirection] << std::endl;
            std::cout << "nu after: " << newNu1 << std::endl;
        }
        subVolumes_[subVolIdx].nu[nuDirection] = newNu1;

        if (subVolumes_[subVolIdx].T != newT)
        {
            std::cout << "T different!!"  << std::endl;
            std::cout << "T before: " << subVolumes_[subVolIdx].T << std::endl;
            std::cout << "T after: " << subVolumes_[subVolIdx].T << std::endl;
        }
        subVolumes_[subVolIdx].T = newT;

    }

private:
    DynamicMatrix T_;
    DynamicMatrix CAInverse_;
    DynamicMatrix AInverseB_;
    DynamicMatrix AInverse_;

    DimVector centralVertex_;
    int centralVertexIndex_;

    std::vector<SubVolumeFace> subVolumeFaces_;
    std::vector<SubVolume> subVolumes_;

    bool stored_;
    bool interiorVolume_;
    int noOfFaces_;
    int noOfSubVols_;
    const int interactionVolumeType_;

    std::set<int> interiorFaces_;
    std::set<int> dirichletFaces_;
    std::map < std::pair<int, int>, int> facesOnSubVolumes_;
    std::map <int, int> globalToSubVolIdx_;
};
}
#endif
