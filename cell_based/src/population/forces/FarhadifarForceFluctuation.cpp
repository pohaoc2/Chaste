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
#include <iostream>
#include <random>
#include "FarhadifarForceFluctuation.hpp"

template<unsigned DIM>
FarhadifarForceFluctuation<DIM>::FarhadifarForceFluctuation()
   : AbstractForce<DIM>(),
     mAreaElasticityParameter(1.0), // These parameters are Case I in Farhadifar's paper
     mPerimeterContractilityParameter(0.04),
     mLineTensionParameter(0.12),
     mPreviousLineTensionParameter(0.12),
     mBoundaryLineTensionParameter(0.12), // This parameter as such does not exist in Farhadifar's model
     mPreviousBoundaryLineTensionParameter(0.12),
     mTargetAreaParameter(1.0),
     mDt(0.01),
     mTau(0.1),
     mSigma(0.1)
{
}

template<unsigned DIM>
FarhadifarForceFluctuation<DIM>::~FarhadifarForceFluctuation()
{
}

template<unsigned DIM>
void FarhadifarForceFluctuation<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{
    // Throw an exception message if not using a VertexBasedCellPopulation
    ///\todo: check whether this line influences profiling tests - if so, we should remove it.
    if (dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation) == nullptr)
    {
        EXCEPTION("FarhadifarForce is to be used with a VertexBasedCellPopulation only");
    }
    if (mline_tension_map.size() == 0 || mline_tension_map.size() != static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation)->rGetMesh().GetNumEdges())
    {
        InitializeLineTensionMap(static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation));
    }
        
    //std::cout << "Line tension of edge 100 = " << mline_tension_map[100] << std::endl;

    VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);
    UpdateLineTensionMap(p_cell_population);

    

    // Define some helper variables

    unsigned num_nodes = p_cell_population->GetNumNodes();
    unsigned num_elements = p_cell_population->GetNumElements();

    /*
     * Check if a subclass of AbstractTargetAreaModifier is being used by
     * interrogating the first cell (assuming either all cells, or no cells, in
     * the population have the CellData item "target area").
     */
    bool using_target_area_modifier = p_cell_population->Begin()->GetCellData()->HasItem("target area");

    // Begin by computing the area and perimeter of each element in the mesh, to avoid having to do this multiple times
    std::vector<double> element_areas(num_elements);
    std::vector<double> element_perimeters(num_elements);
    std::vector<double> target_areas(num_elements);
    for (typename VertexMesh<DIM,DIM>::VertexElementIterator elem_iter = p_cell_population->rGetMesh().GetElementIteratorBegin();
         elem_iter != p_cell_population->rGetMesh().GetElementIteratorEnd();
         ++elem_iter)
    { 
        unsigned elem_index = elem_iter->GetIndex();
        element_areas[elem_index] = p_cell_population->rGetMesh().GetVolumeOfElement(elem_index);
        element_perimeters[elem_index] = p_cell_population->rGetMesh().GetSurfaceAreaOfElement(elem_index);

        if (using_target_area_modifier)
        {
            target_areas[elem_index] = p_cell_population->GetCellUsingLocationIndex(elem_index)->GetCellData()->GetItem("target area");
        }
        else
        {
            target_areas[elem_index] = mTargetAreaParameter;
        }
    }

    // Iterate over vertices in the cell population
    for (unsigned node_index=0; node_index<num_nodes; node_index++)
    {
        Node<DIM>* p_this_node = p_cell_population->GetNode(node_index);

        /*
         * The force on this Node is given by the gradient of the total free
         * energy of the CellPopulation, evaluated at the position of the vertex. This
         * free energy is the sum of the free energies of all CellPtrs in
         * the cell population. The free energy of each CellPtr is comprised of three
         * terms - an area deformation energy, a perimeter deformation energy
         * and line tension energy.
         *
         * Note that since the movement of this Node only affects the free energy
         * of the CellPtrs containing it, we can just consider the contributions
         * to the free energy gradient from each of these CellPtrs.
         */
        c_vector<double, DIM> area_elasticity_contribution = zero_vector<double>(DIM);
        c_vector<double, DIM> perimeter_contractility_contribution = zero_vector<double>(DIM);
        c_vector<double, DIM> line_tension_contribution = zero_vector<double>(DIM);

        // Find the indices of the elements owned by this node
        std::set<unsigned> containing_elem_indices = p_cell_population->GetNode(node_index)->rGetContainingElementIndices();

        // Iterate over these elements
        for (std::set<unsigned>::iterator iter = containing_elem_indices.begin();
             iter != containing_elem_indices.end();
             ++iter)
        {
            // Get this element, its index and its number of nodes
            VertexElement<DIM, DIM>* p_element = p_cell_population->GetElement(*iter);
            unsigned elem_index = p_element->GetIndex();
            unsigned num_nodes_elem = p_element->GetNumNodes();

            // Find the local index of this node in this element
            unsigned local_index = p_element->GetNodeLocalIndex(node_index);

            // Add the force contribution from this cell's area elasticity (note the minus sign)
            c_vector<double, DIM> element_area_gradient =
                    p_cell_population->rGetMesh().GetAreaGradientOfElementAtNode(p_element, local_index);
            area_elasticity_contribution -= GetAreaElasticityParameter()*(element_areas[elem_index] -
                    target_areas[elem_index])*element_area_gradient;

            // Get the previous and next nodes in this element
            unsigned previous_node_local_index = (num_nodes_elem+local_index-1)%num_nodes_elem;
            Node<DIM>* p_previous_node = p_element->GetNode(previous_node_local_index);

            unsigned next_node_local_index = (local_index+1)%num_nodes_elem;
            Node<DIM>* p_next_node = p_element->GetNode(next_node_local_index);

            // Compute the line tension parameter for each of these edges - be aware that this is half of the actual
            // value for internal edges since we are looping over each of the internal edges twice

            double previous_edge_line_tension_parameter = GetLineTensionParameter(p_previous_node, p_this_node, *p_cell_population);
            double next_edge_line_tension_parameter = GetLineTensionParameter(p_this_node, p_next_node, *p_cell_population);

            // Compute the gradient of each these edges, computed at the present node
            c_vector<double, DIM> previous_edge_gradient =
                    -p_cell_population->rGetMesh().GetNextEdgeGradientOfElementAtNode(p_element, previous_node_local_index);
            c_vector<double, DIM> next_edge_gradient = p_cell_population->rGetMesh().GetNextEdgeGradientOfElementAtNode(p_element, local_index);

            // Add the force contribution from cell-cell and cell-boundary line tension (note the minus sign)
            line_tension_contribution -= previous_edge_line_tension_parameter*previous_edge_gradient +
                    next_edge_line_tension_parameter*next_edge_gradient;

            // Add the force contribution from this cell's perimeter contractility (note the minus sign)
            c_vector<double, DIM> element_perimeter_gradient;
            element_perimeter_gradient = previous_edge_gradient + next_edge_gradient;
            perimeter_contractility_contribution -= GetPerimeterContractilityParameter()* element_perimeters[elem_index]*
                                                                                                     element_perimeter_gradient;
        }

        c_vector<double, DIM> force_on_node = area_elasticity_contribution + perimeter_contractility_contribution + line_tension_contribution;
        p_cell_population->GetNode(node_index)->AddAppliedForceContribution(force_on_node);
    }
}

template<unsigned DIM>
void FarhadifarForceFluctuation<DIM>::InitializeLineTensionMap(VertexBasedCellPopulation<DIM>* p_cell_population)
{
    std::vector<double> line_tension_map;
    unsigned numEdges = p_cell_population->rGetMesh().GetNumEdges();
    unsigned num_nodes = p_cell_population->GetNumNodes();
    line_tension_map.resize(numEdges);

    for (unsigned node_index=0; node_index<num_nodes; node_index++)
    {
        Node<DIM>* p_this_node = p_cell_population->GetNode(node_index);
        // Find the indices of the elements owned by this node
        std::set<unsigned> containing_elem_indices = p_cell_population->GetNode(node_index)->rGetContainingElementIndices();

        // Iterate over these elements
        for (std::set<unsigned>::iterator iter = containing_elem_indices.begin();
             iter != containing_elem_indices.end();
             ++iter)
        {
            // Get this element, its index and its number of nodes
            VertexElement<DIM, DIM>* p_element = p_cell_population->GetElement(*iter);
            unsigned num_nodes_elem = p_element->GetNumNodes();
            // Find the local index of this node in this element
            unsigned local_index = p_element->GetNodeLocalIndex(node_index);
            // Get the previous and next nodes in this element
            unsigned previous_node_local_index = (num_nodes_elem+local_index-1)%num_nodes_elem;
            Node<DIM>* p_previous_node = p_element->GetNode(previous_node_local_index);
            std::set<unsigned> shared_elements = GetSharedElements(p_previous_node, p_this_node, *p_cell_population);
            unsigned edgeLocalIndex = GetEdgeLocalIndex(p_previous_node, p_this_node, *p_cell_population, shared_elements);
            if (shared_elements.size() == 1){
                line_tension_map[edgeLocalIndex] = GetBoundaryLineTensionParameter();
            }
            else if (shared_elements.size() == 2){
                line_tension_map[edgeLocalIndex] = GetLineTensionParameter();
            }
        }
    }
    mline_tension_map = line_tension_map;
}


template<unsigned DIM>
std::set<unsigned> FarhadifarForceFluctuation<DIM>::GetSharedElements(Node<DIM>* pNodeA, Node<DIM>* pNodeB, VertexBasedCellPopulation<DIM>& rVertexCellPopulation)
{
    // Find the indices of the elements owned by each node
    std::set<unsigned> elements_containing_nodeA = pNodeA->rGetContainingElementIndices();
    std::set<unsigned> elements_containing_nodeB = pNodeB->rGetContainingElementIndices();

    // Find common elements
    std::set<unsigned> shared_elements;
    std::set_intersection(elements_containing_nodeA.begin(),
                          elements_containing_nodeA.end(),
                          elements_containing_nodeB.begin(),
                          elements_containing_nodeB.end(),
                          std::inserter(shared_elements, shared_elements.begin()));
    assert(!shared_elements.empty());
    return shared_elements;
}

template<unsigned DIM>
unsigned FarhadifarForceFluctuation<DIM>::GetEdgeLocalIndex(Node<DIM>* pNodeA, Node<DIM>* pNodeB, VertexBasedCellPopulation<DIM>& rVertexCellPopulation, std::set<unsigned> shared_elements)
{
    std::set<unsigned>::iterator iter = shared_elements.begin();
    // Get this element, its index and its number of nodes
    VertexElement<DIM, DIM>* p_element = rVertexCellPopulation.GetElement(*iter);
    unsigned num_edges_elem = p_element->GetNumEdges();
    unsigned edgeLocalIndex = 0;
    for (unsigned i = 0; i < num_edges_elem; i++)
    {
        if (p_element->GetEdge(i)->ContainsNode(pNodeA) && p_element->GetEdge(i)->ContainsNode(pNodeB))
        {
            edgeLocalIndex = p_element->GetEdge(i)->GetIndex();
            break;
        }
    }
    return edgeLocalIndex;
}

template<unsigned DIM>
void FarhadifarForceFluctuation<DIM>::CalculatePerturbedLineTension(Node<DIM>* pNodeA, Node<DIM>* pNodeB, VertexBasedCellPopulation<DIM>& rVertexCellPopulation)
{
    std::set<unsigned> shared_elements = GetSharedElements(pNodeA, pNodeB, rVertexCellPopulation);
    unsigned edgeLocalIndex = GetEdgeLocalIndex(pNodeA, pNodeB, rVertexCellPopulation, shared_elements);


    double line_tension_parameter_in_calculation = mline_tension_map[edgeLocalIndex];
    // Add the fluctuation term: truncated normal distribution with mean 0 and variance 1
    std::random_device rd; // Create a random device
    std::default_random_engine generator(rd());
    std::normal_distribution<double> dist(0, mSigma);
    double normal_random = dist(generator);
    double fluctuation = mSigma * std::sqrt(2 * mDt / mTau) * normal_random;
    line_tension_parameter_in_calculation += fluctuation;
    line_tension_parameter_in_calculation = std::max(0.0, line_tension_parameter_in_calculation);
    if (shared_elements.size() == 2)
    {
        line_tension_parameter_in_calculation -= mDt/mTau*(line_tension_parameter_in_calculation - GetLineTensionParameter());
        mline_tension_map[edgeLocalIndex] = line_tension_parameter_in_calculation;
    }
    
    else if (shared_elements.size() == 1)
    {
        line_tension_parameter_in_calculation -= mDt/mTau*(line_tension_parameter_in_calculation - GetBoundaryLineTensionParameter());
        // Add the fluctuation term: truncated normal distribution with mean 0 and variance 1
        mline_tension_map[edgeLocalIndex] = line_tension_parameter_in_calculation;
    }

}

template<unsigned DIM>
void FarhadifarForceFluctuation<DIM>::UpdateLineTensionMap(VertexBasedCellPopulation<DIM>* p_cell_population)
{
    unsigned num_nodes = p_cell_population->GetNumNodes();

    for (unsigned node_index=0; node_index<num_nodes; node_index++)
    {
        Node<DIM>* p_this_node = p_cell_population->GetNode(node_index);
        // Find the indices of the elements owned by this node
        std::set<unsigned> containing_elem_indices = p_cell_population->GetNode(node_index)->rGetContainingElementIndices();

        // Iterate over these elements
        for (std::set<unsigned>::iterator iter = containing_elem_indices.begin();
             iter != containing_elem_indices.end();
             ++iter)
        {
            // Get this element, its index and its number of nodes
            VertexElement<DIM, DIM>* p_element = p_cell_population->GetElement(*iter);
            unsigned num_nodes_elem = p_element->GetNumNodes();
            // Find the local index of this node in this element
            unsigned local_index = p_element->GetNodeLocalIndex(node_index);
            // Get the previous and next nodes in this element
            unsigned previous_node_local_index = (num_nodes_elem+local_index-1)%num_nodes_elem;
            Node<DIM>* p_previous_node = p_element->GetNode(previous_node_local_index);
            CalculatePerturbedLineTension(p_previous_node, p_this_node, *p_cell_population);
        }
    }
}

template<unsigned DIM>
double FarhadifarForceFluctuation<DIM>::GetLineTensionParameter(Node<DIM>* pNodeA, Node<DIM>* pNodeB, VertexBasedCellPopulation<DIM>& rVertexCellPopulation)
{
    std::set<unsigned> shared_elements = GetSharedElements(pNodeA, pNodeB, rVertexCellPopulation);
    unsigned edgeLocalIndex = GetEdgeLocalIndex(pNodeA, pNodeB, rVertexCellPopulation, shared_elements);
    // Since each internal edge is visited twice in the loop above, we have to use half the line tension parameter
    // for each visit.
    double line_tension_parameter_in_calculation = mline_tension_map[edgeLocalIndex];
        if (shared_elements.size() == 2)
    {
        line_tension_parameter_in_calculation = 0.5*line_tension_parameter_in_calculation;
    }


    return line_tension_parameter_in_calculation;
}

template<unsigned DIM>
double FarhadifarForceFluctuation<DIM>::GetAreaElasticityParameter()
{
    return mAreaElasticityParameter;
}

template<unsigned DIM>
double FarhadifarForceFluctuation<DIM>::GetPerimeterContractilityParameter()
{
    return mPerimeterContractilityParameter;
}

template<unsigned DIM>
double FarhadifarForceFluctuation<DIM>::GetLineTensionParameter()
{
    return mLineTensionParameter;
}

template<unsigned DIM>
double FarhadifarForceFluctuation<DIM>::GetBoundaryLineTensionParameter()
{
    return mBoundaryLineTensionParameter;
}

template<unsigned DIM>
double FarhadifarForceFluctuation<DIM>::GetTargetAreaParameter()
{
    return mTargetAreaParameter;
}

template<unsigned DIM>
void FarhadifarForceFluctuation<DIM>::SetAreaElasticityParameter(double areaElasticityParameter)
{
    mAreaElasticityParameter = areaElasticityParameter;
}

template<unsigned DIM>
void FarhadifarForceFluctuation<DIM>::SetPerimeterContractilityParameter(double perimeterContractilityParameter)
{
    mPerimeterContractilityParameter = perimeterContractilityParameter;
}

template<unsigned DIM>
void FarhadifarForceFluctuation<DIM>::SetLineTensionParameter(double lineTensionParameter)
{
    mLineTensionParameter = lineTensionParameter;
    mPreviousLineTensionParameter = lineTensionParameter;
}

template<unsigned DIM>
void FarhadifarForceFluctuation<DIM>::SetBoundaryLineTensionParameter(double boundaryLineTensionParameter)
{
    mBoundaryLineTensionParameter = boundaryLineTensionParameter;
    mPreviousBoundaryLineTensionParameter = boundaryLineTensionParameter;
}

template<unsigned DIM>
void FarhadifarForceFluctuation<DIM>::SetTargetAreaParameter(double targetAreaParameter)
{
    mTargetAreaParameter = targetAreaParameter;
}

template<unsigned DIM>
void FarhadifarForceFluctuation<DIM>::SetDt(double Dt)
{
    mDt = Dt;
}

template<unsigned DIM>
void FarhadifarForceFluctuation<DIM>::SetTau(double tau)
{
    mTau = tau;
}

template<unsigned DIM>
void FarhadifarForceFluctuation<DIM>::SetSigma(double sigma)
{
    mSigma = sigma;
}


template<unsigned DIM>
void FarhadifarForceFluctuation<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<AreaElasticityParameter>" << mAreaElasticityParameter << "</AreaElasticityParameter>\n";
    *rParamsFile << "\t\t\t<PerimeterContractilityParameter>" << mPerimeterContractilityParameter << "</PerimeterContractilityParameter>\n";
    *rParamsFile << "\t\t\t<LineTensionParameter>" << mLineTensionParameter << "</LineTensionParameter>\n";
    *rParamsFile << "\t\t\t<BoundaryLineTensionParameter>" << mBoundaryLineTensionParameter << "</BoundaryLineTensionParameter>\n";
    *rParamsFile << "\t\t\t<TargetAreaParameter>" << mTargetAreaParameter << "</TargetAreaParameter>\n";

    // Call method on direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class FarhadifarForceFluctuation<1>;
template class FarhadifarForceFluctuation<2>;
template class FarhadifarForceFluctuation<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(FarhadifarForceFluctuation)
