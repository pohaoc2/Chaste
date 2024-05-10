/*

Copyright (c) 2005-2024, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef ZULUETAFORCE_HPP_
#define ZULUETAFORCE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "Exception.hpp"

#include "AbstractForce.hpp"
#include "VertexBasedCellPopulation.hpp"

#include <iostream>

/**
 * A force class for use in Vertex-based simulations. This force is based on the
 * Energy function proposed by Farhadifar et al in  Curr. Biol., 2007, 17, 2095-2104.
 */


template<unsigned DIM>
class Zuluetaforce : public AbstractForce<DIM>
{
friend class TestForces;

private:

    friend class boost::serialization::access;
    /**
     * Boost Serialization method for archiving/checkpointing.
     * Archives the object and its member variables.
     *
     * @param archive  The boost archive.
     * @param version  The current version of this class.
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<DIM> >(*this);
        archive & mAreaElasticityParameter;
        archive & mTargetAreaParameter;
    }

protected:

    /**
     * The strength of the area term in the model. Corresponds to K_alpha in Farhadifar's paper.
     */
    double mAreaElasticityParameter;

    /**
     * The target area. This term corresponds to A_0 in Farhadifar's paper.
     *
     * Note: If a subclass of AbstractTargetAreaModifier is present in the
     * cell-based simulation, then we use the CellData item "target area" stored
     * in each Cell in the force calculation instead of this parameter value.
     */
    double mTargetAreaParameter;
    double mbeta;
    double mmu;
    double mgamma;
    double mmc;
    double mAmax;
    double mK;
    double mkplusConstant;
    double mHillCoefficient;
    double mk1;
    double mk2;
    double mDt;

    std::vector<double> mRestLengths;
    std::vector<double> mMyosinLevels;
    std::vector<double> mLineTension;

public:

    /**
     * Constructor.
     */
    Zuluetaforce();

    /**
     * Destructor.
     */
    virtual ~Zuluetaforce();

    /**
     * Overridden AddForceContribution() method.
     *
     * Calculates the force on each node in the vertex-based cell population based on the energy function
     * Farhadifar's model.
     *
     * @param rCellPopulation reference to the cell population
     */
    virtual void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation);

    /**
     * @return mAreaElasticityParameter
     */
    double GetAreaElasticityParameter();

    /**
     * @return mTargetAreaParameter
     */
    double GetTargetAreaParameter();

    /**
     * Set mAreaElasticityParameter.
     *
     * @param areaElasticityParameter the new value of mAreaElasticityParameter
     */
    void SetAreaElasticityParameter(double areaElasticityParameter);

    /**
     * Set mTargetAreaParameter.
     *
     * @param targetAreaParameter the new value of mTargetAreaParameter
     */


    void SetTargetAreaParameter(double targetAreaParameter);
    void InitializeRestLengths(VertexBasedCellPopulation<DIM>* p_cell_population);
    void InitializeMyosinLevel(VertexBasedCellPopulation<DIM>* p_cell_population);
    void InitializeLineTension(VertexBasedCellPopulation<DIM>* p_cell_population);

    void UpdateMyosinLevelLineTension(VertexBasedCellPopulation<DIM>* p_cell_population);
    void CalculateLineTension(unsigned edgeLocalIndex, double edgeLength, std::set<unsigned> shared_elements);
    void CalculateMyosinLevel(unsigned edgeLocalIndex, double edgeLength);
    double GetBoundaryLineTensionParameter(unsigned edgeIndex, double edgeLength);
    double GetNonBoundaryLineTensionParameter(unsigned edgeIndex, double edgeLength);
    std::set<unsigned> GetSharedElements(Node<DIM>* pNodeA, Node<DIM>* pNodeB, VertexBasedCellPopulation<DIM>& rVertexCellPopulation);
    unsigned GetEdgeLocalIndex(Node<DIM>* pNodeA, Node<DIM>* pNodeB, VertexBasedCellPopulation<DIM>& rVertexCellPopulation, std::set<unsigned> shared_elements);
    double GetLineTensionParameter(Node<DIM>* pNodeA, Node<DIM>* pNodeB, VertexBasedCellPopulation<DIM>& rVertexCellPopulation);
    double CalculateMyosinContribution(unsigned edgeIndex, double edgeLength, double deformation);
    
    void SetTensionParameter(double beta, double mu, double gamma);
    void Setmc(double mc);
    void Setkplus(double Amax, double K, double kplusConstant, double HillCoefficient);
    void Setkminus(double k1, double k2);
    void SetDt(double dt);
    

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputForceParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(Zuluetaforce)

#endif /*ZULUETAFORCE_HPP_*/
