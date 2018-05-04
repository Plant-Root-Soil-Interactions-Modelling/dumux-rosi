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
#ifndef DUMUX_MPFAL2DINTERACTIONVOLUME_HH
#define DUMUX_MPFAL2DINTERACTIONVOLUME_HH

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
class MpfaL2DInteractionVolume
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
    typedef typename GridView::template Codim<0>::EntityPointer ElementPointer;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;

    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> WorldMatrix;

public:

    typedef Dune::FieldVector<Scalar, dimWorld> DimVector;
    typedef Dune::FieldMatrix<Scalar, 2, 3> TransmissivityMatrix;

    struct Properties
    {
        // In an interior 2D Mpfa-L interaction volume
        // for quadrilaterals there are 7 elements in stencil
        // This value is used for a primary resizing of involved
        // vectors.
        static const int StandardStencilSize = 7;
    };

    //! Constructs an MpfaO2DInteractionVolume object
    // initialization of the variables according to the
    // case of an interior interaction volume in 2D
    // for quadrilaterals
    MpfaL2DInteractionVolume(bool secondTriangle):
        stored_(false),
        noOfFaces_(2),
        noOfSubVols_(3),
        interiorVolume_(true),
        interactionVolumeType_(MpfaMethods::lMethod),
        secondTriangle_(secondTriangle),
        fluxFaceIndex_(-1)
    {
        R_ = 0;
        R_[0][1] = 1;
        R_[1][0] = -1;
    }

    bool secondTriangle() const
    {
        return secondTriangle_;
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

    //! Returns true if information has already been stored
    bool isStored() const
    {
        return stored_;
    }

    // At the moment the following mappings are not used.
    // TODO: To be erased after testing and verification that not necessary
    void setFaceToFacetMaps(int subVolIdx, int facetIdxOnElement, int faceIndex)
    {
        std::pair<int, int> key(subVolIdx, facetIdxOnElement);
        facesOnSubVolumes_.insert(std::pair< std::pair<int,int>, int> (key, faceIndex) );

        std::pair<int, int> key2(faceIndex, subVolIdx);
        faceToFacetIndex_.insert(std::pair< std::pair<int, int>, int> (key2, facetIdxOnElement) );
    }

    int getFaceIndexFromElement(const Element &element, int facetIdxOnElement)
    {
        int subVolIdx = getSubVolIdx(element);

        std::pair<int, int> key(subVolIdx, facetIdxOnElement);

        std::map< std::pair<int, int>, int >::iterator it = facesOnSubVolumes_.find(key);
        if (it != facesOnSubVolumes_.end())
            return it->second;
        else
            DUNE_THROW(Dune::InvalidStateException, "No face index found corresponding to local subvolume face index!");
    }

    int getFaceIndexFromSubVolume(int subVolIdx, int facetIdxOnElement)
    {
        std::pair<int, int> key(subVolIdx, facetIdxOnElement);

        std::map< std::pair<int, int>, int >::iterator it = facesOnSubVolumes_.find(key);
        if (it != facesOnSubVolumes_.end())
            return it->second;
        else
            DUNE_THROW(Dune::InvalidStateException, "No face index found corresponding to local subvolume face index!");
    }

    int getSubVolIdx(const Element &element) const
    {
        int subVolIdx = -1;
        for(int i = 0; i < 3; i++)
            if(subVolumes_[i] == element)
            {
                subVolIdx = i; break;
            }
        if(subVolIdx != -1)
            return subVolIdx;
        else
            DUNE_THROW(Dune::InvalidStateException, "Element not found in the sub volumes!");
    }

    void setFluxFaceIndex(int index)
    {
        fluxFaceIndex_ = index;
    }

    int getFluxFaceIdx() const
    {
        if (fluxFaceIndex_ == -1)
            DUNE_THROW(Dune::InvalidStateException, "Flux face index not stored correctly!");
        else
            return fluxFaceIndex_;
    }

    int getNumberOfSubVols() const
    {
        return noOfSubVols_;
    }

    int getNumberOfSubFaces() const
    {
        return noOfFaces_;
    }

    bool isInteriorVolume() const
    {
        return interiorVolume_;
    }

    void setNu(DimVector nu, int index)
    {
        nu_[index] = nu;
    }

    DimVector getNu(int index) const
    {
        return nu_[index];
    }

    void setT(Scalar T, int subVolIdx)
    {
        dF_[subVolIdx] = T;
    }

    void setNormal(DimVector normal, int index)
    {
        normal_[index] = normal;
    }

    void setSubVolume(Element element, int subVolIdx)
    {
        subVolumes_[subVolIdx] = element;
    }

    Element getSubVolume(int subVolIdx) const
    {
        return subVolumes_[subVolIdx];
    }

    TransmissivityMatrix& getTransMatrix()
    {
        return T_;
    }

    const TransmissivityMatrix& getTransMatrix() const
    {
        return T_;
    }

    Scalar getW(const Problem &problem, int faceIdx, int subVolIdx, int nuIdx) const
    {
        WorldMatrix K = problem.model().getDiffusionCoefficient(*subVolumes_[subVolIdx]);

        DimVector tmp(0);
        DimVector normal = normal_[faceIdx];
        DimVector nu = nu_[nuIdx];

        K.mv(nu, tmp);
        tmp /= dF_[subVolIdx];
        return normal * tmp;
    }

    Scalar getXi(int nu1Idx, int subVolIdx, int nu2Idx)
    {
        DimVector Rnu2(0);
        R_.umv(nu_[nu2Idx], Rnu2);
        Scalar tmp = Rnu2 * nu_[nu1Idx];
        tmp /= dF_[subVolIdx];
        return tmp;
    }

    void calculateTransmissivityMatrix(const Problem &problem)
    {
        Dune::FieldMatrix<Scalar, 2, 2> A(0);
        Dune::FieldMatrix<Scalar, 2, 3> B(0);
        Dune::FieldMatrix<Scalar, 2, 2> C(0);
        Dune::FieldMatrix<Scalar, 2, 3> D(0);

        Scalar w111 = getW(problem, 0, 0, 0);
        Scalar w112 = getW(problem, 0, 0, 1);
        Scalar w211 = getW(problem, 1, 0, 0);
        Scalar w212 = getW(problem, 1, 0, 1);

        Scalar w124 = getW(problem, 0, 1, 3);
        Scalar w123 = getW(problem, 0, 1, 2);
        Scalar w236 = getW(problem, 1, 2, 5);
        Scalar w235 = getW(problem, 1, 2, 4);

        Scalar X711 = getXi(6, 0, 0);
        Scalar X712 = getXi(6, 0, 1);

        C[0][0] = -w111; C[0][1] = -w112;
        C[1][0] = -w211; C[1][1] = -w212;

        D[0][0] = w111 + w112;
        D[1][0] = w211 + w212;

        A[0][0] = w111 - w124 - w123*X711;
        A[0][1] = w112 - w123*X712;
        A[1][0] = w211 - w236*X711;
        A[1][1] = w212 - w235 - w236*X712;

        B[0][0] = w111 + w112 + w123*(1 - X711 - X712);
        B[0][1] = -w123 - w124;
        B[1][0] = w211 + w212 + w236*(1 - X711 - X712);
        B[1][2] = -w235 - w236;

        //printMatrices(A, B, C, D);

        // Set up Transmissivity Matrix
        A.invert();
        //T_ = C.rightmultiply(B.leftmultiply(A));
        Dune::FieldMatrix<Scalar, 2, 3> temp = B.leftmultiply(A);
        Dune::FieldMatrix<Scalar, 2, 3> temp2(0);
        for(int i = 0; i < 2; i++)
        {
                temp2[i][0] = C[i][0]*temp[0][0] + C[i][1]*temp[1][0];
                temp2[i][1] = C[i][0]*temp[0][1] + C[i][1]*temp[1][1];
                temp2[i][2] = C[i][0]*temp[0][2] + C[i][1]*temp[1][2];

        }
        temp2 += D;
        T_ = temp2;
        //T_ += D;

        //printTransmissivityMatrix();
    }

    void printMatrices(Dune::FieldMatrix<Scalar,2, 2> &A, Dune::FieldMatrix<Scalar,2, 3> &B,
                                        Dune::FieldMatrix<Scalar,2, 2> &C, Dune::FieldMatrix<Scalar,2, 3> &D)
    {
        std::cout << "Printing matrices A, B, C, D of the interaction volume..." << std::endl;

        std::cout << "Matrix A: " << std::endl;
        for(int i = 0; i < A.N();i++)
        {
            for(int j = 0; j < A.M(); j++)
                std::cout << A[i][j] << ", ";
            std::cout << std::endl;
        }
        std::cout << "Matrix B: " << std::endl;
        for(int i = 0; i < B.N();i++)
        {
            for(int j = 0; j < B.M(); j++)
                std::cout << B[i][j] << ", ";
            std::cout << std::endl;
        }
        std::cout << "Matrix C: " << std::endl;
        for(int i = 0; i < C.N();i++)
        {
            for(int j = 0; j < C.M(); j++)
                std::cout << C[i][j] << ", ";
            std::cout << std::endl;
        }
        std::cout << "Matrix D: " << std::endl;
        for(int i = 0; i < D.N();i++)
        {
            for(int j = 0; j < D.M(); j++)
                std::cout << D[i][j] << ", ";
            std::cout << std::endl;
        }
    }

    void printTransmissivityMatrix()
    {
        std::cout << "Printing transmissivity matrix: " << std::endl;
        for(int i = 0; i < 2;i++)
        {
            for(int j = 0; j < 3; j++)
                std::cout << T_[i][j] << ", ";
            std::cout << std::endl;
        }
    }

    void printTransmissivityMatrix() const
    {
        std::cout << "Printing transmissivity matrix: " << std::endl;
        for(int i = 0; i < 2;i++)
        {
            for(int j = 0; j < 3; j++)
                std::cout << T_[i][j] << ", ";
            std::cout << std::endl;
        }
    }

    void printInteractionVolume()
    {
        std::cout << "Printing interaction volume..." << std::endl;
        std::cout << "Sub volume centers..." << std::endl;
        for(int i = 0; i < 3; i++)
        std::cout << "Sub volume " << i+1 << " center: " << subVolumes_[i]->geometry().center() << std::endl;

        std::cout << "Nu vectors..." << std::endl;
        for(int i = 0; i < 7; i++)
            std::cout << "Nu " << i+1 << ": " << nu_[i] << std::endl;

        std::cout << "Normals: " << std::endl;
        for(int i = 0; i < 2; i++)
            std::cout << "normal " << i+1 << ": " << normal_[i] << std::endl;

        std::cout << "dFs... " << std::endl;
        for (int i = 0; i < 3; i++)
            std::cout << "dF-" << i+1 << ": " << dF_[i] << std::endl;
    }

    void passElementsInRegion(std::vector<Element> &storage) const
    {
        storage.resize(3);
        for (int i = 0; i < subVolumes_.size(); i++)
            storage[i] = subVolumes_[i];
    }

private:
    TransmissivityMatrix T_;

    Dune::FieldVector<Element, 3> subVolumes_;
    Dune::FieldVector<DimVector, 7> nu_;
    Dune::FieldVector<DimVector, 2> normal_;
    Dune::FieldVector<Scalar, 3> dF_;

    Dune::FieldMatrix<Scalar, 2, 2> R_;

    bool stored_;
    bool secondTriangle_;
    const bool interiorVolume_;

    int fluxFaceIndex_;
    const int noOfFaces_;
    const int noOfSubVols_;
    const int interactionVolumeType_;

    std::map< std::pair<int, int>, int> facesOnSubVolumes_;
    std::map< std::pair<int, int>, int> faceToFacetIndex_;
};
}
#endif
